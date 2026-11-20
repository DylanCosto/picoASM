//! Bases chosen by how the reads misread each context, learned from the genome itself.
//!
//! ONT reads misread some contexts far more than others, in the same way on every copy: a
//! modified base (AGCT in some *K. pneumoniae*, read as AACT or ACT by both strands' reads)
//! can leave the right base with as little as a quarter of the reads, and a plain vote, which
//! asks only which base most reads have, then follows the misreads. Here the reads at every
//! base of the contigs are counted by what they have there, in their own orientation and by the
//! bases round it ([`BEFORE`] and [`AFTER`]), and so are the bases they insert between two
//! bases ([`GAP_SIDE`]). Nearly every base of a polished contig is right, so these counts are
//! how the reads read each context: AGCT's G as G 40%, A 35%, nothing 15%.
//!
//! Then at each base where the reads disagree with the contig ([`DOUBT`]), each other version
//! (another base, none, one inserted) is weighed against the contig's: each read's piece round
//! the base is aligned to each version by a pair hidden Markov model whose rates are those
//! counts at each base's context (so a read is weighed over all its alignments, not the one its
//! mapping chose), and the version that makes the reads likeliest is taken when it beats the
//! contig's by `min_gain` (a natural log). A contig that has AACT where the strain has AGCT
//! loses: reads of AACT almost never call G there.
//!
//! Weighed only by the reads, a context the genome rarely has can win where the reads misread
//! a common one: in *H. pylori* DSM1 the reads garble the methylated TCTTGA on both strands,
//! the contexts learned from the sites the contigs have right do not foresee the worst-read
//! ones, and TCTCGA (15 sites in the strain) took 233 of its 790 TCTTGAs. So a substitution
//! that would make the genome's 6-mers round it rarer ([`Composition`]) must gain that much
//! more from the reads, if one of the 6-mers it makes is one the reads alone would make at
//! many sites, a good share of what the contigs have ([`MIN_FLOOD`]): weighed at every substitution, the
//! cost kept right bases of rare contexts from being put back in genomes with no such motif
//! (*C. jejuni* 19 errors to 48, *Streptomyces* 5,972 to 6,701). Only rarer: rewarding the common ones as well turned right bases of
//! rare contexts into common ones in genomes the reads read cleanly (*S. flexneri* 6 errors
//! to 120, *P. rettgeri* 2 to 171).
//!
//! A run of one base is lengthened or shortened here too, by a base, the reads aligned to each
//! length: up to [`MAX_INDEL_RUN`], beyond which homopolymer calling's length stands.

use crate::overlap::{Overlap, Reads};
use crate::seq::{Strand, reverse_complement, reverse_complement_into};
use rayon::prelude::*;

use crate::homopolymer::query_positions;
use crate::polish::Draft;
use crate::vote::{BASES, Change, Column, code};

/// The bases before and after a base that make its context, in a read's orientation.
const BEFORE: usize = 3;
const AFTER: usize = 3;

/// The bases either side of a gap between two bases that make its context.
const GAP_SIDE: usize = 3;

/// The shorter contexts that smooth the counts of contexts with few sites.
const SHORT_SIDE: usize = 2;

/// How many reads of a context's own count the shorter context's rates weigh as.
const SMOOTHING: f64 = 20.0;

/// No outcome is taken to be rarer than this.
const MIN_RATE: f64 = 1e-4;

/// A strand's reads count only where at least this many cover a base.
const MIN_STRAND_DEPTH: u32 = 3;

/// A base is weighed again where, on either strand, less than this share of the reads have
/// it, or at least [`1 - DOUBT`] of them insert a base after it.
const DOUBT: f64 = 0.8;

/// Another version is weighed only if at least this share of a strand's reads have it.
const MIN_VERSION_SHARE: f64 = 0.1;

/// Another version is taken when it makes the reads at least this many times as likely
/// (as a natural log) as the contig's version does (the default
/// [`crate::BactParams::context_gain`]). On two read sets of each of the five training genomes,
/// 10 put two false bases into S. flexneri (gains of 10.0 and 10.9); 12, 16 and 20 left 35, 37
/// and 38 errors of the 71 before.
pub(crate) const MIN_GAIN: f64 = 12.0;

/// How many times the fall in the natural logs of the 6-mers' counts round a substitution it
/// costs, against the reads' gain (the default [`crate::BactParams::composition`]). Weighed both
/// ways on *H. pylori* DSM1 (two cultures; errors 641 and 735 with none), 10 left 479 and 469,
/// 20 left 411 and 398, 25 left 423 and 415, 40 left 667; as a cost only, 20 left 415 and 410
/// and changed none of the core set's clean genomes.
pub(crate) const COMPOSITION_WEIGHT: f64 = 20.0;

/// The length of the k-mers the composition counts.
const COMPOSITION_K: usize = 6;

/// The composition weighs a substitution only if one of the k-mers it makes is one the reads
/// alone would make at at least this many sites, and at more than [`MIN_FLOOD_SHARE`] of the
/// number the contigs have: a context the reads misread the same way everywhere (*H. pylori*'s
/// TCTCGA, made at 63 sites over the 9 it had), not a rare one a correction happens to make.
/// 3, 5 and 10 left the same errors.
const MIN_FLOOD: u32 = 10;

/// See [`MIN_FLOOD`]. On *H. pylori* DSM1 (errors 415 and 410 weighing every substitution),
/// 1 left 416 and 458, 0.25 left 414 and 410, 1/16 left 407 and 406 but put 8 errors into
/// *C. difficile*; at 0.25 no 6-mer floods in *C. jejuni*, *Streptomyces*, *C. difficile* or
/// *C. lari*.
const MIN_FLOOD_SHARE: f64 = 0.25;

/// The most rounds of learning and weighing: each learns from the contigs the last corrected.
pub(crate) const ROUNDS: usize = 2;

/// The bases either side of a version that each read's piece spans.
const REACH: usize = 12;

/// A read's piece is weighed only if its length is within this factor of the contig's.
const MAX_LENGTH_RATIO: usize = 2;

/// Insertions and deletions are not weighed in runs of this many of one base, nor ones that
/// make such a run: on six training read sets, weighing them in runs of up to 6 (or up to 10)
/// left 24 errors of 31; longer runs are homopolymer calling's alone.
const MAX_INDEL_RUN: usize = 7;

/// What the reads have at a base: A, C, G, T, nothing.
const BASE_OUTCOMES: usize = 5;

/// What the reads have between two bases: nothing, or an inserted A, C, G, T.
const GAP_OUTCOMES: usize = 5;

/// The bases of contexts beyond a template's ends.
const PAD: usize = if BEFORE > GAP_SIDE { BEFORE } else { GAP_SIDE };

/// How the reads read each context: rates of each outcome, by context code, in the reads'
/// orientation.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Model {
    base: Vec<[f64; BASE_OUTCOMES]>,
    gap: Vec<[f64; GAP_OUTCOMES]>,
}

/// What each strand's reads have at the base of `column`.
fn at_base(column: &Column) -> [[u32; BASE_OUTCOMES]; 2] {
    column.at.map(|s| s.map(u32::from))
}

/// What each strand's reads have between the base of `column` and the next: an inserted base,
/// and nothing for the rest of those covering the base.
fn after_base(column: &Column) -> [[u32; GAP_OUTCOMES]; 2] {
    [0, 1].map(|s| {
        let depth: u32 = column.at[s].iter().map(|&n| u32::from(n)).sum();
        let inserted = column.after[s].map(u32::from);
        let mut out = [0; GAP_OUTCOMES];
        out[0] = depth.saturating_sub(inserted.iter().sum());
        out[1..].copy_from_slice(&inserted);
        out
    })
}

/// The 2-bit code of `seq[from..from + len]` read on `strand` (1 the reverse complement);
/// `None` if a base is not ACGT or the range leaves `seq`.
fn kmer(seq: &[u8], from: usize, len: usize, strand: usize) -> Option<usize> {
    let bases = seq.get(from..from + len)?;
    let mut k = 0;
    if strand == 0 {
        for &b in bases {
            k = k * 4 + code(b)?;
        }
    } else {
        for &b in bases.iter().rev() {
            k = k * 4 + 3 - code(b)?;
        }
    }
    Some(k)
}

/// The context code of the base at `p` of `seq`, read on `strand`, `side` bases either side.
fn base_context(seq: &[u8], p: usize, strand: usize, side: usize) -> Option<usize> {
    kmer(seq, p.checked_sub(side)?, 2 * side + 1, strand)
}

/// The context code of the gap after `p` of `seq`, read on `strand`, `side` bases either side.
fn gap_context(seq: &[u8], p: usize, strand: usize, side: usize) -> Option<usize> {
    kmer(seq, (p + 1).checked_sub(side)?, 2 * side, strand)
}

/// Outcome `k` at a base (or, with `gap`, between two bases), as the reads of `strand` have it
/// in their own orientation.
fn oriented(k: usize, strand: usize, gap: bool) -> usize {
    match (strand, gap) {
        (0, _) => k,
        (_, false) if k < 4 => 3 - k,
        (_, false) => k,
        (_, true) if k == 0 => 0,
        (_, true) => 5 - k,
    }
}

/// Rates from `counts`, each context smoothed towards its shorter context's rates (`short`
/// maps a context to its shorter one).
fn rates<const N: usize>(
    counts: &[[f64; N]],
    short_counts: &[[f64; N]],
    short: impl Fn(usize) -> usize,
) -> Vec<[f64; N]> {
    let mut all = [0.0; N];
    for c in short_counts {
        for (a, n) in all.iter_mut().zip(c) {
            *a += n;
        }
    }
    let total: f64 = all.iter().sum();
    #[expect(clippy::cast_precision_loss, reason = "a small count")]
    let overall = all.map(|n| (n + 1.0) / (total + N as f64));
    let short_rates: Vec<[f64; N]> = short_counts
        .iter()
        .map(|c| {
            let n: f64 = c.iter().sum();
            std::array::from_fn(|k| (c[k] + SMOOTHING * overall[k]) / (n + SMOOTHING))
        })
        .collect();
    counts
        .iter()
        .enumerate()
        .map(|(ctx, c)| {
            let prior = &short_rates[short(ctx)];
            let n: f64 = c.iter().sum();
            std::array::from_fn(|k| ((c[k] + SMOOTHING * prior[k]) / (n + SMOOTHING)).max(MIN_RATE))
        })
        .collect()
}

/// How often each [`COMPOSITION_K`]-mer appears in the contigs, either strand.
#[derive(Clone, Debug, PartialEq)]
struct Composition {
    counts: Vec<u32>,
}

impl Composition {
    /// Counts the k-mers of `seqs` and their reverse complements.
    fn count(seqs: &[&[u8]]) -> Self {
        let mut counts = vec![0; 1 << (2 * COMPOSITION_K)];
        for seq in seqs {
            for i in 0..seq.len().saturating_sub(COMPOSITION_K - 1) {
                for strand in 0..2 {
                    if let Some(k) = kmer(seq, i, COMPOSITION_K, strand) {
                        counts[k] += 1;
                    }
                }
            }
        }
        Self { counts }
    }

    /// The change in the sum of ln(count + 1) over the k-mers that cover `p` of `seq` when it
    /// is `base`; 0 near the ends.
    fn gain(&self, seq: &[u8], p: usize, base: u8) -> f64 {
        let Some((window, changed)) = windows(seq, p, base) else {
            return 0.0;
        };
        let score = |w: &[u8]| -> f64 {
            (0..COMPOSITION_K)
                .filter_map(|i| kmer(w, i, COMPOSITION_K, 0))
                .map(|k| f64::from(self.counts[k] + 1).ln())
                .sum()
        };
        score(&changed) - score(window)
    }

    /// Which k-mers the substitutions in `made` (counts of the k-mers they make, either
    /// strand) would flood the contigs with: made at least [`MIN_FLOOD`] times, and more than
    /// [`MIN_FLOOD_SHARE`] of the times the contigs have them.
    fn flooded(&self, made: &[u32]) -> Vec<bool> {
        made.iter()
            .zip(&self.counts)
            .map(|(&m, &c)| m >= MIN_FLOOD && f64::from(m) > MIN_FLOOD_SHARE * f64::from(c))
            .collect()
    }
}

/// The k-mers that cover `p` of `seq` once it is `base`, read on `strand`.
fn made(seq: &[u8], p: usize, base: u8, strand: usize) -> Vec<usize> {
    windows(seq, p, base).map_or_else(Vec::new, |(_, changed)| {
        (0..COMPOSITION_K)
            .filter_map(|i| kmer(&changed, i, COMPOSITION_K, strand))
            .collect()
    })
}

/// The bases of `seq` whose [`COMPOSITION_K`]-mers cover `p`, as they are and with `base` at
/// `p`; `None` near the ends.
fn windows(seq: &[u8], p: usize, base: u8) -> Option<(&[u8], Vec<u8>)> {
    let from = (p + 1).checked_sub(COMPOSITION_K)?;
    let window = seq.get(from..p + COMPOSITION_K)?;
    let mut changed = window.to_vec();
    changed[COMPOSITION_K - 1] = base;
    Some((window, changed))
}

impl Model {
    /// Learns how the reads read each context from the `columns` of `seqs`.
    pub(crate) fn learn(seqs: &[&[u8]], columns: &[Vec<Column>]) -> Self {
        let mut base = vec![[0.0; BASE_OUTCOMES]; 1 << (2 * (BEFORE + AFTER + 1))];
        let mut base_short = vec![[0.0; BASE_OUTCOMES]; 1 << (2 * (2 * SHORT_SIDE + 1))];
        let mut gap = vec![[0.0; GAP_OUTCOMES]; 1 << (4 * GAP_SIDE)];
        let mut gap_short = vec![[0.0; GAP_OUTCOMES]; 1 << (4 * SHORT_SIDE)];
        for (seq, columns) in seqs.iter().zip(columns) {
            for (p, column) in columns.iter().enumerate() {
                let (at, after) = (at_base(column), after_base(column));
                for s in 0..2 {
                    if at[s].iter().sum::<u32>() < MIN_STRAND_DEPTH {
                        continue;
                    }
                    if let (Some(c), Some(short)) = (
                        base_context(seq, p, s, BEFORE),
                        base_context(seq, p, s, SHORT_SIDE),
                    ) {
                        for (k, &n) in at[s].iter().enumerate() {
                            let o = oriented(k, s, false);
                            base[c][o] += f64::from(n);
                            base_short[short][o] += f64::from(n);
                        }
                    }
                    if let (Some(c), Some(short)) = (
                        gap_context(seq, p, s, GAP_SIDE),
                        gap_context(seq, p, s, SHORT_SIDE),
                    ) {
                        for (k, &n) in after[s].iter().enumerate() {
                            let o = oriented(k, s, true);
                            gap[c][o] += f64::from(n);
                            gap_short[short][o] += f64::from(n);
                        }
                    }
                }
            }
        }
        // A base context's shorter one drops a base each side; a gap's, one each side.
        let base_short_of = |c: usize| (c >> 2) & ((1 << (2 * (2 * SHORT_SIDE + 1))) - 1);
        let gap_short_of = |c: usize| (c >> 2) & ((1 << (4 * SHORT_SIDE)) - 1);
        Self {
            base: rates(&base, &base_short, base_short_of),
            gap: rates(&gap, &gap_short, gap_short_of),
        }
    }

    /// The log likelihood of `read` being read off `template[pad..len - pad]` (the bases
    /// either side give the contexts), summed over all alignments: each base of the template
    /// is read as a base or not at all, and then any bases are inserted after it, at the rates
    /// of its context and the gap's.
    fn likelihood(&self, template: &[u8], pad: usize, read: &[u8]) -> f64 {
        let m = template.len() - 2 * pad;
        let n = read.len();
        // `before[i]`: the read's first `i` bases read off the template's bases so far.
        let mut before = vec![0.0; n + 1];
        before[0] = 1.0;
        let mut at = vec![0.0; n + 1];
        let codes: Vec<Option<usize>> = read.iter().map(|&b| code(b)).collect();
        for j in 0..m {
            let t = pad + j;
            let base =
                base_context(template, t, 0, BEFORE).map_or([0.2; BASE_OUTCOMES], |c| self.base[c]);
            let gap = gap_context(template, t, 0, GAP_SIDE)
                .map_or([0.96, 0.01, 0.01, 0.01, 0.01], |c| self.gap[c]);
            for i in 0..=n {
                // The template base read as the read's base `i - 1`, or not at all.
                let read_as = if i > 0 {
                    before[i - 1] * codes[i - 1].map_or(MIN_RATE, |k| base[k])
                } else {
                    0.0
                };
                at[i] = read_as + before[i] * base[4];
                // Then the read's base `i - 1` inserted after it.
                if i > 0 {
                    at[i] += at[i - 1] * codes[i - 1].map_or(MIN_RATE, |k| gap[1 + k]);
                }
            }
            for i in 0..=n {
                before[i] = at[i] * gap[0];
            }
        }
        before[n].max(f64::MIN_POSITIVE).ln()
    }
}

/// A version of a contig at one base.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Version {
    /// The contig's own.
    Contig,
    /// The base is another.
    Base(u8),
    /// The base is not there.
    Deleted,
    /// A base is inserted after it.
    Inserted(u8),
    /// It and the next base change places.
    Swapped,
}

impl Version {
    /// The change to `seq` that makes this version at `p`.
    fn change(self, seq: &[u8], p: usize) -> Option<Change> {
        let start = u32::try_from(p).expect("contigs are shorter than 4 Gb");
        match self {
            Self::Contig => None,
            Self::Base(b) => Some(Change {
                start,
                len: 1,
                bases: vec![b],
            }),
            Self::Deleted => Some(Change {
                start,
                len: 1,
                bases: Vec::new(),
            }),
            Self::Inserted(b) => Some(Change {
                start: start + 1,
                len: 0,
                bases: vec![b],
            }),
            Self::Swapped => Some(Change {
                start,
                len: 2,
                bases: vec![seq[p + 1], seq[p]],
            }),
        }
    }
}

/// Returns the length of the run of one base that holds `seq[p]`.
fn run_len(seq: &[u8], p: usize) -> usize {
    let base = seq[p];
    let before = seq[..p].iter().rev().take_while(|&&b| b == base).count();
    let after = seq[p + 1..].iter().take_while(|&&b| b == base).count();
    before + 1 + after
}

/// Returns each strand's share of its reads that `count` counts, for strands with at least
/// [`MIN_STRAND_DEPTH`] reads.
fn shares(column: &Column, count: impl Fn(usize) -> u32) -> impl Iterator<Item = f64> {
    (0..2).filter_map(move |s| {
        let depth: u32 = column.at[s].iter().map(|&n| u32::from(n)).sum();
        (depth >= MIN_STRAND_DEPTH).then(|| f64::from(count(s)) / f64::from(depth))
    })
}

/// The versions other than the contig's worth weighing at `p` of `seq` (see [`DOUBT`],
/// [`MIN_VERSION_SHARE`]).
fn versions(seq: &[u8], p: usize, column: &Column) -> Vec<Version> {
    let Some(current) = code(seq[p]) else {
        return Vec::new();
    };
    let doubtful = shares(column, |s| u32::from(column.at[s][current])).any(|x| x < DOUBT)
        || shares(column, |s| {
            column.after[s].iter().map(|&n| u32::from(n)).sum()
        })
        .any(|x| x >= 1.0 - DOUBT);
    if !doubtful {
        return Vec::new();
    }
    let common =
        |count: &dyn Fn(usize) -> u32| shares(column, count).any(|x| x >= MIN_VERSION_SHARE);
    let mut out = Vec::new();
    for (k, &base) in BASES.iter().enumerate() {
        if k != current && common(&|s| u32::from(column.at[s][k])) {
            out.push(Version::Base(base));
        }
    }
    if run_len(seq, p) < MAX_INDEL_RUN && common(&|s| u32::from(column.at[s][4])) {
        out.push(Version::Deleted);
    }
    for (k, &base) in BASES.iter().enumerate() {
        if common(&|s| u32::from(column.after[s][k])) {
            // The run the inserted base joins: its neighbours' runs and itself.
            let left = if seq[p] == base { run_len(seq, p) } else { 0 };
            let right = match seq.get(p + 1) {
                Some(&b) if b == base => run_len(seq, p + 1),
                _ => 0,
            };
            if left + right + 1 < MAX_INDEL_RUN {
                out.push(Version::Inserted(base));
            }
        }
    }
    // The base and the next the other way round (a polished AGCTTTTT, whose C most reads
    // miss, as AGTCTTTT).
    if seq.get(p + 1).is_some_and(|&b| b != seq[p]) {
        out.push(Version::Swapped);
    }
    out
}

/// A base weighed again: where, its other versions, and the reads' pieces round it (each in
/// the contig's orientation, with the strand it was read on).
#[derive(Clone, Debug)]
struct Site {
    p: usize,
    versions: Vec<Version>,
    pieces: Pieces,
}

/// Reads' pieces, in the order added, their bases one after another: held as a vector each,
/// the pieces of a noisy run's sites took a third of its memory at the peak (*P. kielensis*,
/// 2.9 Gb of ONT reads 6.7% from the contigs).
#[derive(Clone, Debug, Default)]
struct Pieces {
    bases: Vec<u8>,
    /// Each piece's strand, and where its bases end in `bases`.
    ends: Vec<(Strand, u32)>,
}

impl Pieces {
    fn push(&mut self, strand: Strand, piece: &[u8]) {
        self.bases.extend_from_slice(piece);
        let end = u32::try_from(self.bases.len()).expect("a site's pieces are short");
        self.ends.push((strand, end));
    }

    fn len(&self) -> usize {
        self.ends.len()
    }

    /// Each piece with its strand.
    fn iter(&self) -> impl Iterator<Item = (Strand, &[u8])> {
        let starts = std::iter::once(0).chain(self.ends.iter().map(|&(_, end)| end));
        self.ends
            .iter()
            .zip(starts)
            .map(|(&(strand, end), start)| (strand, &self.bases[start as usize..end as usize]))
    }
}

impl<P: AsRef<[u8]>> FromIterator<(Strand, P)> for Pieces {
    fn from_iter<I: IntoIterator<Item = (Strand, P)>>(iter: I) -> Self {
        let mut pieces = Self::default();
        for (strand, piece) in iter {
            pieces.push(strand, piece.as_ref());
        }
        pieces
    }
}

/// The template of `version` at `p` of `seq`: the bases from `p - REACH` to `p + REACH`, with
/// [`PAD`] more either side.
fn template(seq: &[u8], p: usize, version: Version) -> Vec<u8> {
    let mut t = seq[p - REACH - PAD..=p + REACH + PAD].to_vec();
    let at = REACH + PAD;
    match version {
        Version::Contig => {}
        Version::Base(b) => t[at] = b,
        Version::Deleted => {
            t.remove(at);
        }
        Version::Inserted(b) => t.insert(at + 1, b),
        Version::Swapped => t.swap(at, at + 1),
    }
    t
}

/// The log likelihood of the reads' `pieces` under `version` at `p` of `seq`.
fn weigh(model: &Model, seq: &[u8], p: usize, version: Version, pieces: &Pieces) -> f64 {
    let forward = template(seq, p, version);
    let reverse = reverse_complement(&forward);
    let mut read = Vec::new();
    pieces
        .iter()
        .map(|(strand, piece)| match strand {
            Strand::Forward => model.likelihood(&forward, PAD, piece),
            Strand::Reverse => {
                reverse_complement_into(piece, &mut read);
                model.likelihood(&reverse, PAD, &read)
            }
        })
        .sum()
}

/// Returns, per draft, the sites worth weighing, with the pieces of the reads (`mappings`)
/// that span them.
fn sites(
    drafts: &[Draft],
    columns: &[Vec<Column>],
    reads: &Reads,
    mappings: &[Option<Overlap>],
) -> Vec<Vec<Site>> {
    let mut sites: Vec<Vec<Site>> = drafts
        .iter()
        .zip(columns)
        .map(|(d, columns)| {
            let seq = &d.seq;
            (REACH + PAD..seq.len().saturating_sub(REACH + PAD + 1))
                .filter_map(|p| {
                    let versions = versions(seq, p, &columns[p]);
                    (!versions.is_empty()).then(|| Site {
                        p,
                        versions,
                        pieces: Pieces::default(),
                    })
                })
                .collect()
        })
        .collect();
    let mut oriented = Vec::new();
    let span = 2 * REACH + 1;
    for (read, mapping) in reads.indices().zip(mappings) {
        let Some(m) = mapping else { continue };
        let contig = &mut sites[m.target.index()];
        let first = contig.partition_point(|s| s.p < pos(m.target_start) + REACH);
        let last = contig.partition_point(|s| s.p + REACH < pos(m.target_end));
        if first >= last {
            continue;
        }
        let bases = reads.get(read);
        let oriented_start = match m.strand {
            Strand::Forward => {
                oriented.clear();
                oriented.extend_from_slice(bases);
                m.query_start
            }
            Strand::Reverse => {
                reverse_complement_into(bases, &mut oriented);
                reads.read_len(read) - m.query_end
            }
        };
        let at = query_positions(m, oriented_start);
        for site in &mut contig[first..last] {
            let from = u32::try_from(site.p - REACH).expect("contigs are shorter than 4 Gb");
            let to = u32::try_from(site.p + REACH + 1).expect("contigs are shorter than 4 Gb");
            let (Some(a), Some(b)) = (at(from), at(to)) else {
                continue;
            };
            let (a, b) = (pos(a), pos(b));
            if b < a || (b - a) * MAX_LENGTH_RATIO < span || b - a > span * MAX_LENGTH_RATIO {
                continue;
            }
            site.pieces.push(m.strand, &oriented[a..b]);
        }
    }
    sites
}

/// Returns, per draft, the changes the reads (`mappings`, and their `columns`) favour by at
/// least `min_gain` (see the module docs), sorted by start, none within twice [`REACH`] bases
/// of another (whose reads it was weighed with); the model is learned from all the drafts.
pub(crate) fn changes(
    drafts: &[Draft],
    columns: &[Vec<Column>],
    reads: &Reads,
    mappings: &[Option<Overlap>],
    min_gain: f64,
    composition_weight: f64,
) -> Vec<Vec<Change>> {
    let seqs: Vec<&[u8]> = drafts.iter().map(|d| d.seq.as_slice()).collect();
    let model = Model::learn(&seqs, columns);
    let sites = sites(drafts, columns, reads, mappings);
    // What the reads alone gain by each version of each site.
    let by_reads: Vec<Vec<Vec<(f64, Version)>>> = sites
        .iter()
        .zip(&seqs)
        .map(|(sites, seq)| {
            sites
                .par_iter()
                .map(|site| {
                    let current = weigh(&model, seq, site.p, Version::Contig, &site.pieces);
                    site.versions
                        .iter()
                        .map(|&v| (weigh(&model, seq, site.p, v, &site.pieces) - current, v))
                        .collect()
                })
                .collect()
        })
        .collect();
    let composition = Composition::count(&seqs);
    let mut made_by_reads = vec![0; 1 << (2 * COMPOSITION_K)];
    for ((sites, gains), seq) in sites.iter().zip(&by_reads).zip(&seqs) {
        for (site, gains) in sites.iter().zip(gains) {
            if let Some(&(_, Version::Base(b))) = best(gains).filter(|g| g.0 >= min_gain) {
                for strand in 0..2 {
                    for k in made(seq, site.p, b, strand) {
                        made_by_reads[k] += 1;
                    }
                }
            }
        }
    }
    let flooded = composition.flooded(&made_by_reads);
    log::debug!(
        "bacteria: context: {} {COMPOSITION_K}-mers the reads would flood the contigs with",
        flooded.iter().filter(|&&f| f).count()
    );
    sites
        .into_iter()
        .zip(by_reads)
        .zip(&seqs)
        .map(|((sites, by_reads), seq)| {
            let mut found: Vec<(f64, Change)> = sites
                .par_iter()
                .zip(by_reads)
                .filter_map(|(site, by_reads)| {
                    let gains: Vec<(f64, Version)> = by_reads
                        .into_iter()
                        .map(|(gain, v)| {
                            let prior = match v {
                                Version::Base(b)
                                    if made(seq, site.p, b, 0).iter().any(|&k| flooded[k]) =>
                                {
                                    (composition_weight * composition.gain(seq, site.p, b)).min(0.0)
                                }
                                _ => 0.0,
                            };
                            (gain + prior, v)
                        })
                        .collect();
                    log::trace!(
                        "bacteria: context at {}: {}|{}|{}, {} reads, gains {gains:?}",
                        site.p,
                        String::from_utf8_lossy(&seq[site.p - REACH..site.p]),
                        char::from(seq[site.p]),
                        String::from_utf8_lossy(&seq[site.p + 1..=site.p + REACH]),
                        site.pieces.len()
                    );
                    let &(gain, version) = best(&gains).filter(|b| b.0 >= min_gain)?;
                    Some((gain, version.change(seq, site.p)?))
                })
                .collect();
            // The changes that gain most first; one weighed with reads round another is left.
            found.sort_by(|a, b| b.0.total_cmp(&a.0));
            let mut kept: Vec<Change> = Vec::new();
            for (_, c) in found {
                if kept
                    .iter()
                    .all(|k| pos(k.start.abs_diff(c.start)) > 2 * REACH)
                {
                    kept.push(c);
                }
            }
            kept.sort_by_key(|c| c.start);
            kept
        })
        .collect()
}

/// The version that gains most.
fn best(gains: &[(f64, Version)]) -> Option<&(f64, Version)> {
    gains.iter().max_by(|a, b| a.0.total_cmp(&b.0))
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A pseudo-random sequence of `n` bases.
    fn genome(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                BASES[usize::try_from(x % 4).expect("small")]
            })
            .collect()
    }

    /// A model learned from `seq` read with 5% of each base read as nothing, and AGCT's G read
    /// as G 40%, A 40%, nothing 20% on both strands.
    fn model(seq: &[u8]) -> Model {
        let columns: Vec<Column> = (0..seq.len())
            .map(|p| {
                let mut at = [[0_u16; 5]; 2];
                let k = code(seq[p]).expect("ACGT");
                let modified = p >= 1 && p + 2 < seq.len() && &seq[p - 1..p + 3] == b"AGCT";
                let modified_rev = p >= 2 && p + 1 < seq.len() && &seq[p - 2..p + 2] == b"AGCT";
                for (s, m) in [(0, modified), (1, modified_rev)] {
                    if m {
                        at[s][k] = 8;
                        at[s][if s == 0 { 0 } else { 3 }] = 8;
                        at[s][4] = 4;
                    } else {
                        at[s][k] = 19;
                        at[s][4] = 1;
                    }
                }
                Column {
                    at,
                    after: [[0; 4]; 2],
                }
            })
            .collect();
        Model::learn(&[seq], &[columns])
    }

    #[test]
    fn a_base_that_makes_rare_kmers_loses_composition() {
        // TCTTGA many times, TCTCGA once: making a TCTTGA a TCTCGA loses, the reverse gains as
        // much.
        let mut seq = Vec::new();
        for i in 0..100 {
            seq.extend_from_slice(if i == 50 {
                b"AATCTCGAAT"
            } else {
                b"AATCTTGAAT"
            });
        }
        let composition = Composition::count(&[&seq]);
        let p = 20 * 10 + 5;
        assert_eq!(&seq[p - 3..p + 3], b"TCTTGA");
        let to_c = composition.gain(&seq, p, b'C');
        assert!(to_c < -10.0, "TCTTGA to TCTCGA: {to_c}");
        let q = 50 * 10 + 5;
        let to_t = composition.gain(&seq, q, b'T');
        assert!((to_t + to_c).abs() < 1e-9, "{to_t} vs {to_c}");
        assert!(
            composition.gain(&seq, 2, b'C').abs() < f64::EPSILON,
            "near the start"
        );
    }

    #[test]
    fn only_kmers_the_reads_would_make_everywhere_are_flooded() {
        // TCTCGA at 2 sites (and TCGAGA, its reverse complement, on the other strand); made at
        // 10 more it floods, at 3 it does not, nor does TCTTGA made at 10 over its 98.
        let mut seq = Vec::new();
        for i in 0..100 {
            seq.extend_from_slice(if i < 2 { b"AATCTCGAAT" } else { b"AATCTTGAAT" });
        }
        let composition = Composition::count(&[&seq]);
        let rare = kmer(b"TCTCGA", 0, COMPOSITION_K, 0).expect("ACGT");
        let common = kmer(b"TCTTGA", 0, COMPOSITION_K, 0).expect("ACGT");
        let mut by_reads = vec![0; 1 << (2 * COMPOSITION_K)];
        by_reads[rare] = MIN_FLOOD;
        by_reads[common] = MIN_FLOOD;
        let flooded = composition.flooded(&by_reads);
        assert!(flooded[rare]);
        assert!(!flooded[common]);
        by_reads[rare] = 3;
        assert!(!composition.flooded(&by_reads)[rare]);
        // Making a TCTTGA a TCTCGA makes TCTCGA.
        assert!(made(&seq, 20 * 10 + 5, b'C', 0).contains(&rare));
    }

    #[test]
    fn pieces_come_back_as_they_were_added() {
        let added = [
            (Strand::Forward, b"ACGTA".to_vec()),
            (Strand::Reverse, b"GG".to_vec()),
            (Strand::Forward, b"TTTACG".to_vec()),
        ];
        let pieces: Pieces = added.iter().cloned().collect();
        assert_eq!(pieces.len(), 3);
        let back: Vec<(Strand, Vec<u8>)> = pieces.iter().map(|(s, p)| (s, p.to_vec())).collect();
        assert_eq!(back, added);
    }

    #[test]
    fn reads_of_a_modified_base_favour_it_over_the_base_they_misread_it_as() {
        let mut truth = genome(20_000, 7);
        for i in (100..19_900).step_by(97) {
            truth[i..i + 4].copy_from_slice(b"AGCT");
        }
        let model = model(&truth);
        let p = 1071;
        assert_eq!(&truth[p - 1..p + 3], b"AGCT");
        let mut contig = truth.clone();
        contig[p] = b'A';
        // Reads of the truth: 4 read the G as G, 4 as A, 2 skip it.
        let piece = |b: Option<u8>| {
            let mut r = truth[p - REACH..=p + REACH].to_vec();
            match b {
                Some(b) => r[REACH] = b,
                None => {
                    r.remove(REACH);
                }
            }
            r
        };
        let mut pieces = Pieces::default();
        for (b, n) in [(Some(b'G'), 4), (Some(b'A'), 4), (None, 2)] {
            for _ in 0..n {
                pieces.push(Strand::Forward, &piece(b));
            }
        }
        let a = weigh(&model, &contig, p, Version::Contig, &pieces);
        let g = weigh(&model, &contig, p, Version::Base(b'G'), &pieces);
        assert!(g > a + MIN_GAIN, "G {g} vs A {a}");
        // And 20 reads of AACT, all reading A, favour A.
        let pieces: Pieces = vec![(Strand::Forward, piece(Some(b'A'))); 20]
            .into_iter()
            .collect();
        let a = weigh(&model, &contig, p, Version::Contig, &pieces);
        let g = weigh(&model, &contig, p, Version::Base(b'G'), &pieces);
        assert!(a > g + MIN_GAIN, "A {a} vs G {g}");
    }

    #[test]
    fn two_bases_the_contig_has_the_wrong_way_round_are_swapped() {
        let truth = genome(20_000, 17);
        let model = model(&truth);
        let p = (3000..4000)
            .find(|&p| truth[p] != truth[p + 1])
            .expect("two different bases");
        let mut contig = truth.clone();
        contig.swap(p, p + 1);
        let pieces: Pieces = vec![(Strand::Forward, &truth[p - REACH..=p + REACH]); 8]
            .into_iter()
            .collect();
        let current = weigh(&model, &contig, p, Version::Contig, &pieces);
        let swapped = weigh(&model, &contig, p, Version::Swapped, &pieces);
        assert!(swapped > current + MIN_GAIN, "{swapped} vs {current}");
        assert_eq!(
            Version::Swapped.change(&contig, p),
            Some(Change {
                start: u32::try_from(p).expect("small"),
                len: 2,
                bases: truth[p..p + 2].to_vec()
            })
        );
    }

    #[test]
    fn a_missing_base_is_put_back() {
        let truth = genome(20_000, 13);
        let model = model(&truth);
        // A base outside runs, which the contig lacks; every read has it.
        let p = (2000..3000)
            .find(|&p| {
                truth[p - 1] != truth[p] && truth[p + 1] != truth[p] && truth[p - 1] != truth[p + 1]
            })
            .expect("such a base");
        let mut contig = truth.clone();
        contig.remove(p);
        // The site is the base before it; the reads span its reach and the missing base.
        let pieces: Pieces = vec![(Strand::Reverse, &truth[p - 1 - REACH..=p + REACH]); 8]
            .into_iter()
            .collect();
        let current = weigh(&model, &contig, p - 1, Version::Contig, &pieces);
        let inserted = weigh(&model, &contig, p - 1, Version::Inserted(truth[p]), &pieces);
        assert!(inserted > current + MIN_GAIN, "{inserted} vs {current}");
    }
}
