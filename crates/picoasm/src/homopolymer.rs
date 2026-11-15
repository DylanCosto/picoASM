//! Homopolymer lengths called again from the reads, after polishing.
//!
//! Most of the errors polishing leaves in an Oxford Nanopore assembly are homopolymers a base
//! or more short: reads under-call the length of a long run, the more so the longer it is, so
//! the length most reads give is short, and a consensus of them is too. Here the reads are
//! mapped onto the polished contigs once more, and each read's length is measured at every run
//! of at least [`MIN_RUN`] bases (the longest stretch of the run's base in the read's bases
//! aligned from the base before the run to the one after it, both included), each strand's
//! reads apart. How reads give runs of each length, by the run's base as they read it
//! ([`Lengths`]), is learned from the contigs' own runs: reads that see a run as A or G give
//! its length far more often than those that see it as T or C (a run of 8 G's: 74% of G reads
//! give 8, 56% of C reads), so a run is learned from only where its A or G reads mostly give
//! its length, and a run called short does not teach that runs of its length read longer. Then
//! each run is given the length that makes its reads likeliest, if that beats its own by
//! [`MIN_LENGTH_GAIN`]; the learning and calling are done again from the lengths called
//! ([`ROUNDS`]). The changes are applied with those of `vote.rs`, from the same mapping.

use picoasm_overlap::{EditOp, Overlap, ReadIndex, Reads};
use picoasm_seq::{Strand, reverse_complement_into};
use rayon::prelude::*;

use crate::polish::Draft;
use crate::vote::Change;

/// Runs shorter than this are left as polishing called them.
pub(crate) const MIN_RUN: usize = 4;

/// A run's reads are measured only if there are at least this many.
const MIN_READS: usize = 5;

/// Lengths are counted up to this; longer ones count as this.
const MAX_LEN: usize = 63;

/// A run of one base in a contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    /// The contig's index.
    pub contig: usize,
    /// Its first base.
    pub start: u32,
    /// Its length.
    pub len: u32,
    /// Its base.
    pub base: u8,
}

/// The reads' lengths at one run (for development: see [`crate::BactParams::keep_layout`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Evidence {
    /// The run, as polishing left it.
    pub run: Run,
    /// `counts[k]`: the reads that give the run `k` bases.
    pub counts: Vec<u32>,
    /// The same, for the reads of each strand (forward, reverse) apart.
    pub by_strand: [Vec<u32>; 2],
}

/// Measures the reads' lengths at every run of the drafts (see the module docs), from the
/// reads' `mappings` onto them ([`picoasm_overlap::map_reads`]).
pub(crate) fn evidence(
    drafts: &[Draft],
    reads: &Reads,
    mappings: &[Option<Overlap>],
) -> Vec<Evidence> {
    let runs: Vec<Vec<Run>> = drafts
        .iter()
        .enumerate()
        .map(|(c, d)| runs_of(c, &d.seq))
        .collect();
    let mut counts: Vec<Vec<[Vec<u32>; 2]>> = runs
        .iter()
        .map(|r| vec![[vec![0; MAX_LEN + 1], vec![0; MAX_LEN + 1]]; r.len()])
        .collect();
    // Each read's lengths are measured in parallel, a batch of reads at a time, and counted
    // (in any order: they are sums) as each batch comes back.
    let mapped: Vec<(ReadIndex, &Overlap)> = reads
        .indices()
        .zip(mappings)
        .filter_map(|(read, m)| m.as_ref().map(|m| (read, m)))
        .collect();
    for batch in mapped.chunks(EVIDENCE_BATCH) {
        let lengths: Vec<Vec<(usize, usize)>> = batch
            .par_iter()
            .map_init(Vec::new, |oriented, &(read, m)| {
                lengths_of(&runs[m.target.index()], reads, read, m, oriented)
            })
            .collect();
        for (&(_, m), lengths) in batch.iter().zip(lengths) {
            let strand = usize::from(m.strand.is_reverse());
            for (i, len) in lengths {
                counts[m.target.index()][i][strand][len] += 1;
            }
        }
    }
    runs.into_iter()
        .zip(counts)
        .flat_map(|(runs, counts)| {
            runs.into_iter()
                .zip(counts)
                .map(|(run, by_strand)| Evidence {
                    run,
                    counts: by_strand[0]
                        .iter()
                        .zip(&by_strand[1])
                        .map(|(a, b)| a + b)
                        .collect(),
                    by_strand,
                })
                .filter(|e| e.counts.iter().sum::<u32>() as usize >= MIN_READS)
        })
        .collect()
}

/// Mapped reads whose run lengths are measured at once in [`evidence`].
const EVIDENCE_BATCH: usize = 8_192;

/// Returns the length `read` (mapped by `m`, onto the contig whose runs are `contig_runs`)
/// gives each run it spans, as (the run's index, the length); `oriented` is a buffer.
fn lengths_of(
    contig_runs: &[Run],
    reads: &Reads,
    read: ReadIndex,
    m: &Overlap,
    oriented: &mut Vec<u8>,
) -> Vec<(usize, usize)> {
    let first = contig_runs.partition_point(|r| r.start < m.target_start + 1);
    let last = contig_runs.partition_point(|r| r.start + r.len + 1 < m.target_end);
    if first >= last {
        return Vec::new();
    }
    let bases = reads.get(read);
    let oriented_start = match m.strand {
        Strand::Forward => {
            oriented.clear();
            oriented.extend_from_slice(bases);
            m.query_start
        }
        Strand::Reverse => {
            reverse_complement_into(bases, oriented);
            reads.read_len(read) - m.query_end
        }
    };
    let at = query_positions(m, oriented_start);
    let mut lengths = Vec::new();
    for (i, run) in contig_runs.iter().enumerate().take(last).skip(first) {
        let (Some(from), Some(to)) = (at(run.start - 1), at(run.start + run.len + 1)) else {
            continue;
        };
        let segment = &oriented[pos(from)..pos(to).max(pos(from))];
        lengths.push((i, longest_stretch(segment, run.base).min(MAX_LEN)));
    }
    lengths
}

/// How the reads give runs of each length: `rates[b][len][k]`, the share of the reads that
/// read a run of `len` of base `b` (in their own orientation, ACGT) as `k` bases.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Lengths {
    rates: Vec<Vec<Vec<f64>>>,
}

/// A length's rates weigh its own reads and, as this many reads, the rates of one shorter
/// shifted by a base.
const LENGTH_SMOOTHING: f64 = 50.0;

/// No length is taken to be read rarer than this.
const MIN_LENGTH_RATE: f64 = 1e-4;

/// A run is given another length when that makes its reads at least this many times as likely
/// (a natural log). On six training read sets, 5, 10 and 20 left 22, 23 and 26 errors.
const MIN_LENGTH_GAIN: f64 = 5.0;

/// The rounds of learning and calling.
const ROUNDS: usize = 3;

/// A run is learned from only if at least this many of its reads see it as A or G.
const MIN_PURINE_READS: u32 = 3;

/// A run's length is weighed up to this many bases either side of its length.
const MAX_LENGTH_MOVE: usize = 3;

/// The base of `run`, in the orientation of the reads of `strand`.
fn oriented_base(base: u8, strand: usize) -> Option<usize> {
    let b = crate::vote::code(base)?;
    Some(if strand == 0 { b } else { 3 - b })
}

impl Lengths {
    /// Learns how the reads give runs of each length from `evidence`, the runs having
    /// `lengths` bases.
    fn learn(evidence: &[Evidence], lengths: &[usize]) -> Self {
        let mut counts = vec![vec![vec![0.0; MAX_LEN + 1]; MAX_LEN + 1]; 4];
        for (e, &len) in evidence.iter().zip(lengths) {
            let len = len.min(MAX_LEN);
            // Only runs whose reads that see them as A or G mostly give them that length.
            let Some(purine) =
                (0..2).find(|&s| matches!(oriented_base(e.run.base, s), Some(0 | 2)))
            else {
                continue;
            };
            let reads = &e.by_strand[purine];
            let mode = (0..reads.len()).max_by_key(|&k| (reads[k], std::cmp::Reverse(k)));
            if reads.iter().sum::<u32>() < MIN_PURINE_READS || mode != Some(len) {
                continue;
            }
            for s in 0..2 {
                let Some(b) = oriented_base(e.run.base, s) else {
                    continue;
                };
                for (k, &n) in e.by_strand[s].iter().enumerate() {
                    counts[b][len][k] += f64::from(n);
                }
            }
        }
        let rates = counts
            .into_iter()
            .map(|by_len| {
                let mut rates: Vec<Vec<f64>> = Vec::with_capacity(MAX_LEN + 1);
                for c in by_len {
                    let n: f64 = c.iter().sum();
                    // The rates of one shorter, a base longer; below the shortest, flat.
                    let prior: Vec<f64> = match rates.last() {
                        Some(shorter) => (0..=MAX_LEN)
                            .map(|k| if k == 0 { shorter[0] } else { shorter[k - 1] })
                            .collect(),
                        None => vec![1.0 / 64.0; MAX_LEN + 1],
                    };
                    let r: Vec<f64> = (0..=MAX_LEN)
                        .map(|k| (c[k] + LENGTH_SMOOTHING * prior[k]) / (n + LENGTH_SMOOTHING))
                        .collect();
                    rates.push(r);
                }
                rates
                    .into_iter()
                    .map(|r| r.into_iter().map(|x| x.max(MIN_LENGTH_RATE).ln()).collect())
                    .collect()
            })
            .collect();
        Self { rates }
    }

    /// The log likelihood of a run's reads (`e`) under its having `len` bases.
    fn likelihood(&self, e: &Evidence, len: usize) -> f64 {
        let mut sum = 0.0;
        for s in 0..2 {
            let Some(b) = oriented_base(e.run.base, s) else {
                continue;
            };
            for (k, &n) in e.by_strand[s].iter().enumerate() {
                if n > 0 {
                    sum += f64::from(n) * self.rates[b][len][k];
                }
            }
        }
        sum
    }

    /// The length the reads of `e` favour, by at least [`MIN_LENGTH_GAIN`] over its own
    /// (runs longer than [`MAX_LEN`] keep theirs).
    fn call(&self, e: &Evidence) -> usize {
        let current = pos(e.run.len);
        if current > MAX_LEN {
            return current;
        }
        let own = self.likelihood(e, current);
        (current.saturating_sub(MAX_LENGTH_MOVE).max(MIN_RUN)
            ..=(current + MAX_LENGTH_MOVE).min(MAX_LEN))
            .map(|len| (self.likelihood(e, len) - own, len))
            .filter(|&(gain, _)| gain >= MIN_LENGTH_GAIN)
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .map_or(current, |(_, len)| len)
    }
}

/// Returns, per draft, the runs the reads (`evidence`) give another length (see the module
/// docs), as changes to its bases.
pub(crate) fn changes(drafts: &[Draft], evidence: &[Evidence]) -> Vec<Vec<Change>> {
    let mut lengths: Vec<usize> = evidence.iter().map(|e| pos(e.run.len)).collect();
    for _ in 0..ROUNDS {
        let model = Lengths::learn(evidence, &lengths);
        lengths = evidence.iter().map(|e| model.call(e)).collect();
    }
    let mut out = vec![Vec::new(); drafts.len()];
    for (e, &len) in evidence.iter().zip(&lengths) {
        if len != pos(e.run.len) {
            out[e.run.contig].push(Change {
                start: e.run.start,
                len: e.run.len,
                bases: vec![e.run.base; len],
            });
        }
    }
    out
}

/// Returns the runs of at least [`MIN_RUN`] bases of `seq`, with a base either side.
fn runs_of(contig: usize, seq: &[u8]) -> Vec<Run> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < seq.len() {
        let mut j = i + 1;
        while j < seq.len() && seq[j] == seq[i] {
            j += 1;
        }
        if j - i >= MIN_RUN && i > 0 && j < seq.len() {
            out.push(Run {
                contig,
                start: u32::try_from(i).expect("contigs are shorter than 4 Gb"),
                len: u32::try_from(j - i).expect("contigs are shorter than 4 Gb"),
                base: seq[i],
            });
        }
        i = j;
    }
    out
}

/// Returns a function giving, for a target position the mapping `m` spans, the position in
/// the oriented read aligned to it (as `polish::cut_points` places window ends).
pub(crate) fn query_positions(
    m: &Overlap,
    oriented_start: u32,
) -> impl Fn(u32) -> Option<u32> + '_ {
    // Runs of target positions with their query offsets, one per edit.
    let mut steps: Vec<(u32, u32, EditOp, u32)> = Vec::with_capacity(m.edits.len());
    let (mut t, mut q) = (m.target_start, oriented_start);
    for e in &m.edits {
        steps.push((t, q, e.op(), e.len()));
        match e.op() {
            EditOp::Match | EditOp::Mismatch => {
                t += e.len();
                q += e.len();
            }
            EditOp::Insertion => q += e.len(),
            EditOp::Deletion => t += e.len(),
        }
    }
    let end = (t, q);
    move |target: u32| {
        if target < m.target_start || target > end.0 {
            return None;
        }
        if target == end.0 {
            return Some(end.1);
        }
        // The last step starting at or before `target` that covers target bases.
        let i = steps.partition_point(|s| s.0 <= target).checked_sub(1)?;
        let (mut i, mut found) = (i, None);
        loop {
            let (t0, q0, op, len) = steps[i];
            match op {
                EditOp::Match | EditOp::Mismatch if target < t0 + len => {
                    found = Some(q0 + (target - t0));
                }
                EditOp::Deletion if target < t0 + len => found = Some(q0),
                _ => {}
            }
            if found.is_some() || i == 0 {
                break;
            }
            i -= 1;
        }
        found
    }
}

/// Returns the length of the longest stretch of `base` in `seq`.
fn longest_stretch(seq: &[u8], base: u8) -> usize {
    seq.split(|&b| b != base)
        .map(<[u8]>::len)
        .max()
        .unwrap_or(0)
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_are_found_with_a_base_either_side() {
        let runs = runs_of(0, b"AAAACGGGGGTCCCCA");
        let found: Vec<(u32, u32, u8)> = runs.iter().map(|r| (r.start, r.len, r.base)).collect();
        assert_eq!(found, [(5, 5, b'G'), (11, 4, b'C')]);
    }

    /// A run of `len` of `base` at `start`, its forward reads giving `fwd` and its reverse
    /// reads `rev` (pairs of length and reads).
    fn evidence(
        start: u32,
        len: u32,
        base: u8,
        fwd: &[(usize, u32)],
        rev: &[(usize, u32)],
    ) -> Evidence {
        let mut by_strand = [vec![0; MAX_LEN + 1], vec![0; MAX_LEN + 1]];
        for (s, counts) in [fwd, rev].into_iter().enumerate() {
            for &(k, n) in counts {
                by_strand[s][k] = n;
            }
        }
        Evidence {
            run: Run {
                contig: 0,
                start,
                len,
                base,
            },
            counts: by_strand[0]
                .iter()
                .zip(&by_strand[1])
                .map(|(a, b)| a + b)
                .collect(),
            by_strand,
        }
    }

    #[test]
    fn a_run_is_given_the_length_runs_whose_reads_read_alike_have() {
        // Runs of 6, 7 and 8 G's: G reads give them right, C reads one short or two.
        let mut all = Vec::new();
        for i in 0..200 {
            all.push(evidence(
                i,
                6,
                b'G',
                &[(6, 19), (5, 1)],
                &[(6, 15), (5, 4), (4, 1)],
            ));
        }
        for i in 0..200 {
            all.push(evidence(
                i,
                7,
                b'G',
                &[(7, 18), (6, 2)],
                &[(7, 12), (6, 6), (5, 2)],
            ));
            all.push(evidence(
                i,
                8,
                b'G',
                &[(8, 17), (7, 3)],
                &[(8, 9), (7, 8), (6, 3)],
            ));
        }
        // A run polished two short, its reads like those of 8.
        all.push(evidence(
            5000,
            6,
            b'G',
            &[(8, 16), (7, 4)],
            &[(8, 8), (7, 9), (6, 3)],
        ));
        // A run of 7 whose C reads give 6 more often than most.
        all.push(evidence(
            6000,
            7,
            b'G',
            &[(7, 18), (6, 2)],
            &[(7, 9), (6, 9), (5, 2)],
        ));
        let mut lengths: Vec<usize> = all.iter().map(|e| pos(e.run.len)).collect();
        for _ in 0..ROUNDS {
            let model = Lengths::learn(&all, &lengths);
            lengths = all.iter().map(|e| model.call(e)).collect();
        }
        assert_eq!(lengths[600], 8);
        assert_eq!(lengths[601], 7);
        assert!(
            lengths[..600]
                .iter()
                .zip(&all)
                .all(|(&l, e)| l == pos(e.run.len))
        );
    }

    #[test]
    fn runs_called_another_length_become_changes() {
        let drafts = vec![Draft::new(b"ACGTAAAACGTGGGGGCA".to_vec(), false)];
        // Runs of 4 and 5 A's and G's read right, and one of 4 A's whose reads give 5.
        let mut evidence = Vec::new();
        for (len, base) in [(4, b'A'), (5, b'A'), (4, b'G'), (5, b'G')] {
            for i in 0..100 {
                evidence.push(evidence_at(100 + i, len, base));
            }
        }
        evidence.push(evidence_at(4, 4, b'A'));
        let n = evidence.len();
        evidence[n - 1] = evidence_like(4, 4, b'A', 5);
        evidence.push(evidence_like(11, 5, b'G', 5));
        let c = changes(&drafts, &evidence);
        assert_eq!(
            c[0],
            [Change {
                start: 4,
                len: 4,
                bases: b"AAAAA".to_vec()
            }]
        );
    }

    /// A run of `len` at `start` read right by 10 reads of each strand.
    fn evidence_at(start: u32, len: u32, base: u8) -> Evidence {
        evidence_like(start, len, base, pos(len))
    }

    /// A run of `len` at `start` read as `read` bases by 10 reads of each strand.
    fn evidence_like(start: u32, len: u32, base: u8, read: usize) -> Evidence {
        evidence(start, len, base, &[(read, 10)], &[(read, 10)])
    }
}
