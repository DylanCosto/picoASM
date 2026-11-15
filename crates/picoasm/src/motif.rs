//! Bases that one strand's reads misread: modified bases, learned from the reads.
//!
//! A strain modifies the bases of a few motifs (methylation: AGCT in some *K. pneumoniae*, GATC
//! in most enterobacteria), and ONT reads of the strand that carries the modification often
//! call another base there. Where the reads of that strand outnumber the others, polishing and
//! the plain vote follow them. The motifs differ from strain to strain, so they are learned
//! from the reads themselves: at every base the reads of one strand call cleanly (at least
//! [`CLEAN_SHARE`] of them as the contig has it), the bases round it, in the other strand's
//! reads' own orientation, are counted, and so are the times those reads call another base
//! there ([`NOISY_SHARE`]). A context of a few bases where that strand's reads misread the base
//! far more often than the genome's average ([`ENRICHMENT`]) is a learned motif. Then where the
//! strands disagree, one strand calling the contig's base and most of the other's reads
//! another, and the other base would make the first strand's context a motif it misreads that
//! way, the other base is taken, unless the other strand misreading the contig's base would
//! explain the split about as well (by the rates the motifs are misread at).

use crate::vote::{Column, code};

/// Contexts round a base: so many bases before it and after it, in a read's own orientation.
/// The longest pins a motif; the shorter ones gather the sites of motifs with fewer fixed
/// bases (GATC's) or a short genome's.
const SHAPES: [(usize, usize); 5] = [(2, 2), (1, 2), (2, 1), (1, 3), (3, 1)];

/// No shape reaches further than this from its base.
const REACH: usize = 3;

/// Each strand has at least this many reads at a base learned from or voted on.
const MIN_DEPTH: u32 = 5;

/// A strand calls a base cleanly when at least this share of its reads (in percent) do.
const CLEAN_SHARE: u32 = 80;

/// A strand misreads a base when at least this share of its reads (in percent) call one
/// other base.
const NOISY_SHARE: u32 = 25;

/// A base is corrected when the strand misreading it in a learned motif explains a split at
/// least this many times as often (by the motifs' rates) as the other strand misreading the
/// contig's base would.
const RATE_RATIO: f32 = 2.0;

/// A motif is learned from at least this many sites...
const MIN_SITES: u32 = 30;

/// ...where its misreading is at least this many times as common as the shape's average.
const ENRICHMENT: u64 = 3;

/// Whether the base vote learns modified-base motifs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Modified {
    /// Learned and corrected (see the module docs).
    Learn,
    /// Left to the plain vote.
    Ignore,
}

/// The motifs learned: for each shape, each context's rate of misreading as each base `x` of
/// `ACGT` (the share of its sites where a strand's reads do; 0 where none is learned).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Motifs {
    misread: Vec<Vec<[f32; 4]>>,
}

/// Each strand's reads at a column: how many, and how many call each base or none.
fn depths(column: &Column) -> [u32; 2] {
    [0, 1].map(|s| column.at[s].iter().map(|&n| u32::from(n)).sum())
}

/// The 2-bit codes of the context of shape `(before, after)` round `p` of `seq`, with `base`
/// at `p`, in the orientation of the reads of `strand` (0 forward, 1 reverse); `None` if a
/// base there is not `ACGT`.
fn context(
    seq: &[u8],
    p: usize,
    strand: usize,
    (before, after): (usize, usize),
    base: usize,
) -> Option<usize> {
    let mut k = 0;
    for i in 0..=before + after {
        let b = if i == before {
            base
        } else if strand == 0 {
            code(seq[p + i - before])?
        } else {
            3 - code(seq[p + before - i])?
        };
        let b = if i == before && strand == 1 { 3 - b } else { b };
        k = k << 2 | b;
    }
    Some(k)
}

impl Motifs {
    /// Learns the motifs from the reads' calls at the bases of `seqs` (see the module docs);
    /// `columns[i]` are the reads' calls on `seqs[i]`.
    pub(crate) fn learn(seqs: &[&[u8]], columns: &[Vec<Column>]) -> Self {
        let mut occurs: Vec<Vec<u32>> = SHAPES
            .iter()
            .map(|&(a, b)| vec![0; 1 << (2 * (a + b + 1))])
            .collect();
        let mut sites: Vec<Vec<[u32; 4]>> = SHAPES
            .iter()
            .map(|&(a, b)| vec![[0; 4]; 1 << (2 * (a + b + 1))])
            .collect();
        for (seq, columns) in seqs.iter().zip(columns) {
            for p in REACH..seq.len().saturating_sub(REACH) {
                let Some(b) = code(seq[p]) else { continue };
                let column = &columns[p];
                let depth = depths(column);
                if depth.iter().any(|&d| d < MIN_DEPTH) {
                    continue;
                }
                for s in 0..2 {
                    let o = 1 - s;
                    let clean = u32::from(column.at[o][b]) * 100 >= depth[o] * CLEAN_SHARE;
                    let misread = (0..4)
                        .filter(|&x| x != b)
                        .max_by_key(|&x| column.at[s][x])
                        .filter(|&x| {
                            clean && u32::from(column.at[s][x]) * 100 >= depth[s] * NOISY_SHARE
                        });
                    for (i, &shape) in SHAPES.iter().enumerate() {
                        let Some(k) = context(seq, p, s, shape, b) else {
                            continue;
                        };
                        occurs[i][k] += 1;
                        if let Some(x) = misread {
                            // In the reads' own orientation.
                            sites[i][k][if s == 0 { x } else { 3 - x }] += 1;
                        }
                    }
                }
            }
        }
        let misread = occurs
            .iter()
            .zip(&sites)
            .map(|(occurs, sites)| {
                let all: u64 = occurs.iter().map(|&n| u64::from(n)).sum();
                let misread: u64 = sites.iter().flatten().map(|&n| u64::from(n)).sum();
                occurs
                    .iter()
                    .zip(sites)
                    .map(|(&occurs, sites)| {
                        [0, 1, 2, 3].map(|x| {
                            let n = u64::from(sites[x]);
                            // n / occurs >= ENRICHMENT * misread / all, without dividing.
                            let learned = sites[x] >= MIN_SITES
                                && n * all >= ENRICHMENT * misread * u64::from(occurs);
                            #[expect(
                                clippy::cast_precision_loss,
                                reason = "a rate, not an exact count"
                            )]
                            let rate = n as f32 / occurs.max(1) as f32;
                            if learned { rate } else { 0.0 }
                        })
                    })
                    .collect()
            })
            .collect();
        Self { misread }
    }

    /// How many contexts are learned to be misread (as one base or more).
    pub(crate) fn len(&self) -> usize {
        self.misread
            .iter()
            .flatten()
            .filter(|rates| rates.iter().any(|&r| r > 0.0))
            .count()
    }

    /// Returns the base (`ACGT` code) that most reads of one strand at `p` of `seq` call and
    /// those of the other misread as `seq[p]` in a learned motif; `None` if there is none.
    pub(crate) fn corrected(&self, seq: &[u8], p: usize, column: &Column) -> Option<usize> {
        if self.misread.is_empty() || p < REACH || p + REACH >= seq.len() {
            return None;
        }
        let have = code(seq[p])?;
        let depth = depths(column);
        if depth.iter().any(|&d| d < MIN_DEPTH) {
            return None;
        }
        for s in 0..2 {
            let o = 1 - s;
            // Strand `s` calls the contig's base most; strand `o` another, more than that one.
            if (0..5).any(|k| column.at[s][k] > column.at[s][have]) {
                continue;
            }
            let Some(y) = (0..4).max_by_key(|&y| (column.at[o][y], std::cmp::Reverse(y))) else {
                continue;
            };
            // Most of strand `o`'s reads: a learned misread is against a clean strand.
            if y == have || 2 * u32::from(column.at[o][y]) < depth[o] {
                continue;
            }
            // Strand `s` misreading `y` as the contig's base explains the split; strand `o`
            // misreading the contig's base as `y` may explain it too (a motif on both strands,
            // such as C. jejuni's CCTGA and TCAGG, is learned either way round). The likelier
            // explanation, by the learned rates, decides: by [`RATE_RATIO`] or the contig keeps
            // its base.
            let rate = |strand: usize, truth: usize, read_as: usize| {
                let x = if strand == 0 { read_as } else { 3 - read_as };
                SHAPES
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &shape)| {
                        context(seq, p, strand, shape, truth).map(|k| self.misread[i][k][x])
                    })
                    .fold(0.0_f32, f32::max)
            };
            let (change, keep) = (rate(s, y, have), rate(o, have, y));
            if change > 0.0 && change >= RATE_RATIO * keep {
                return Some(y);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn column(fwd: [u16; 5], rev: [u16; 5]) -> Column {
        Column {
            at: [fwd, rev],
            after: [[0; 4]; 2],
        }
    }

    #[test]
    fn contexts_are_read_in_each_strands_orientation() {
        // AGCTC round the C at 2: forward reads see AG[C]TC; reverse reads, GA[G]CT.
        let seq = b"AGCTC";
        let code_of = |s: &[u8]| s.iter().fold(0, |k, &b| k << 2 | code(b).expect("ACGT"));
        assert_eq!(context(seq, 2, 0, (2, 2), 1), Some(code_of(b"AGCTC")));
        assert_eq!(context(seq, 2, 1, (2, 2), 1), Some(code_of(b"GAGCT")));
        // With T put at 2, and a shorter shape.
        assert_eq!(context(seq, 2, 0, (1, 2), 3), Some(code_of(b"GTTC")));
        assert_eq!(context(seq, 2, 1, (1, 2), 3), Some(code_of(b"AACT")));
    }

    #[test]
    fn a_motif_one_strand_misreads_is_learned_and_its_misreads_corrected() {
        // A sequence with many AGCTs, forward reads of their C calling T a third of the time,
        // and elsewhere clean reads.
        let unit = b"AAGCTGACCATTGCAGTTGCA";
        let seq: Vec<u8> = unit.repeat(60);
        let mut columns = vec![column([0; 5], [0; 5]); seq.len()];
        for (p, &b) in seq.iter().enumerate() {
            let b = code(b).expect("ACGT");
            columns[p].at[0][b] = 20;
            columns[p].at[1][b] = 20;
            if p % unit.len() == 3 {
                // The C of AGCT.
                columns[p].at[0] = [0, 13, 0, 7, 0];
            }
        }
        let motifs = Motifs::learn(&[&seq], &[columns.clone()]);
        assert_ne!(motifs.len(), 0);
        // A contig with T where the reverse reads have C and forward ones mostly T.
        let mut wrong = seq.clone();
        let p = 10 * unit.len() + 3;
        wrong[p] = b'T';
        let split = column([0, 7, 0, 12, 0], [0, 15, 0, 3, 0]);
        assert_eq!(motifs.corrected(&wrong, p, &split), Some(1));
        // Elsewhere (a T of ATTG) the same split is left to the plain vote.
        let q = 10 * unit.len() + 10;
        assert_eq!(
            motifs.corrected(&seq, q, &column([0, 7, 0, 12, 0], [0, 15, 0, 3, 0])),
            None
        );
    }

    #[test]
    fn a_split_two_motifs_explain_goes_to_the_likelier() {
        // The contig has T at the third base of AAGCTG..., forward reads say T and reverse
        // ones C. Either the forward strand misreads C as T there (the contig is wrong), or
        // the reverse strand misreads T as C (in its own orientation, A as G; the contig is
        // right). Each explanation has a learned rate; the likelier one decides.
        let seq = b"AAGCTGACCATTGCAGTTGCA".repeat(3);
        let p = 3;
        let split = column([0, 7, 0, 12, 0], [0, 15, 0, 3, 0]);
        let with = |forward: f32, reverse: f32| {
            let mut motifs = Motifs {
                misread: SHAPES
                    .iter()
                    .map(|&(a, b)| vec![[0.0; 4]; 1 << (2 * (a + b + 1))])
                    .collect(),
            };
            let (i, shape) = (0, SHAPES[0]);
            // Forward reads of AG[C]TG misread the C as T.
            let k = context(&seq, p, 0, shape, 1).expect("ACGT");
            motifs.misread[i][k][3] = forward;
            // Reverse reads of the T version misread (in their orientation, A as G).
            let mut other = seq.clone();
            other[p] = b'T';
            let k = context(&other, p, 1, shape, 3).expect("ACGT");
            motifs.misread[i][k][2] = reverse;
            motifs
        };
        // Forward misreading C as T is far likelier: C is taken.
        let mut wrong = seq.clone();
        wrong[p] = b'T';
        assert_eq!(with(0.4, 0.05).corrected(&wrong, p, &split), Some(1));
        // About as likely the other way: the contig keeps its base.
        assert_eq!(with(0.4, 0.3).corrected(&wrong, p, &split), None);
    }
}
