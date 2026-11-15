//! Overlaps between reads: seeded by shared minimizers, chained, and aligned base by base to
//! the reads' ends.
//!
//! Every read is sketched into homopolymer-compressed universe minimizers ([`Sketcher`]) and
//! indexed. A minimizer held more than [`OverlapParams::max_occurrences`] times seeds nothing:
//! inside a repeat whose copies differ at a few bases, the minimizers all copies share are
//! frequent, and those spanning a copy's own bases are about as frequent as the read depth, so
//! what is left seeds overlaps between reads of one copy (hifiasm's high-frequency filter).
//! For each read, the seeds it shares with each read after it on one relative strand are
//! chained (co-linear, near one diagonal); the best chains are aligned piece by piece between
//! their seeds and extended to the reads' ends with [`DiffAligner`]. A chain whose extension
//! fails on either side (an internal match of a shared repeat; the extension gives up soon after
//! the reads' sequences part) or that diverges by more than [`OverlapParams::max_divergence`] is
//! no overlap.

use std::cmp::Reverse;
use std::collections::{BTreeSet, BinaryHeap};

use picoasm_seq::{Anchor, InvalidSketchParamsError, SketchParams, Sketcher, Strand};
use rayon::prelude::*;
use thiserror::Error;

use crate::align::{DiffAligner, Edit, EditOp, Extent, Stall, push};
use crate::reads::{ReadIndex, Reads};

/// Settings for [`find_overlaps`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OverlapParams {
    /// Minimizer length, in homopolymer-compressed bases (at most 31).
    pub kmer_len: u8,
    /// The share of compressed k-mers kept as minimizers.
    pub density: f64,
    /// A minimizer held more times than this in the reads seeds nothing.
    pub max_occurrences: u32,
    /// A chain needs at least this many seeds.
    pub min_seeds: usize,
    /// Each read aligns at most this many chains, the best first.
    pub max_candidates: usize,
    /// An overlap differs at no more than this share of its columns.
    pub max_divergence: f64,
    /// An overlap covers at least this many bases of each read.
    pub min_overlap: u32,
    /// How far an overlap's alignment must reach.
    pub extension: Extension,
    /// A piece of the alignment between two seeds may differ at this share of its columns
    /// (`None`: [`OverlapParams::max_divergence`]); the overlap as a whole must still be within
    /// [`OverlapParams::max_divergence`]. Raw Oxford Nanopore reads have stretches of a few
    /// hundred bases with errors at 20% or more, which would otherwise lose whole overlaps.
    pub max_piece_divergence: Option<f64>,
    /// Two chained seeds may leave their diagonal by this many bases, plus a twentieth of the
    /// gap between them ([`DIAGONAL_SLACK`] by default). Mapping reads onto a draft contig
    /// takes more, so a read chains across a stretch the draft lacks or holds twice and its
    /// alignment shows the indel.
    pub chain_slack: u32,
    /// Chains are taken from one pair's seeds on one strand at most this many times
    /// ([`CHAINS_PER_PAIR`] by default): long sequences that share many copies of a repeat
    /// need more for the overlap to be among them.
    pub chains_per_pair: usize,
    /// [`map_reads`] only: a chain that does not align is split up to this many times at its
    /// widest gap between seeds, and its parts tried (0, the default: none is split).
    pub max_chain_splits: u32,
    /// A diagonal of the aligner that lags the furthest one by more than this many bases
    /// (target and query together) is dropped (`None`, the default: none is, and each piece
    /// has the fewest edits). Between reads that each differ from the genome at 15-20% the
    /// aligner's work grows with the square of a piece's edits; a band of a few hundred bases
    /// cuts it to a few times the piece's length and hardly changes the alignments. An indel
    /// longer than about half the band is not crossed as the fewest edits would cross it.
    pub max_lag: Option<u32>,
}

/// How far an overlap's alignment must reach past its chain of seeds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Extension {
    /// To an end of the reads on both sides: an overlap, or no overlap at all.
    #[default]
    Full,
    /// As far as the extensions get: an extension that cannot be aligned ends at its seed (a
    /// read's end with a burst of errors, or an adapter, is left out of the overlap, and the
    /// caller decides what the overhang means).
    Clipped,
}

impl Default for OverlapParams {
    fn default() -> Self {
        Self {
            kmer_len: 31,
            density: 0.05,
            max_occurrences: 200,
            min_seeds: 3,
            max_candidates: 400,
            max_divergence: 0.05,
            min_overlap: 1_000,
            extension: Extension::Full,
            max_piece_divergence: None,
            chain_slack: DIAGONAL_SLACK,
            chains_per_pair: CHAINS_PER_PAIR,
            max_chain_splits: 0,
            max_lag: None,
        }
    }
}

/// The settings given to [`find_overlaps`] cannot be used.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InvalidOverlapParamsError {
    /// The minimizer length or density is out of range.
    #[error("invalid minimizer settings")]
    Sketch(#[from] InvalidSketchParamsError),
    /// The divergence is not in `[0, 1)`.
    #[error("maximum divergence {0} is not in [0, 1)")]
    Divergence(f64),
    /// A chain needs at least one seed.
    #[error("a chain needs at least one seed")]
    MinSeeds,
}

/// An overlap of two reads: the query's bases (reverse-complemented on [`Strand::Reverse`])
/// aligned to the target's.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Overlap {
    /// The read aligned to the target.
    pub query: ReadIndex,
    /// The read the query is aligned to.
    pub target: ReadIndex,
    /// [`Strand::Reverse`] if the query is reverse-complemented.
    pub strand: Strand,
    /// The query bases covered, on the query's forward strand.
    pub query_start: u32,
    /// One past the last query base covered, on the query's forward strand.
    pub query_end: u32,
    /// The first target base covered.
    pub target_start: u32,
    /// One past the last target base covered.
    pub target_end: u32,
    /// The alignment's columns in target order.
    pub edits: Vec<Edit>,
}

impl Overlap {
    /// Returns the number of alignment columns.
    #[must_use]
    pub fn columns(&self) -> u64 {
        self.edits.iter().map(|e| u64::from(e.len())).sum()
    }

    /// Returns the number of columns that are not matches.
    #[must_use]
    pub fn differences(&self) -> u64 {
        self.edits
            .iter()
            .filter(|e| e.op() != EditOp::Match)
            .map(|e| u64::from(e.len()))
            .sum()
    }

    /// Replaces the alignment with its matches, then its differences as mismatches: one run
    /// each (more past [`Edit::MAX_LEN`]). [`Self::columns`] and [`Self::differences`] stay as
    /// they were; which columns differ, and how, is lost.
    pub fn keep_counts(&mut self) {
        let differences = self.differences();
        let matches = self.columns() - differences;
        let mut edits = Vec::with_capacity(2);
        for (op, mut left) in [(EditOp::Match, matches), (EditOp::Mismatch, differences)] {
            while left > 0 {
                let len = left.min(u64::from(Edit::MAX_LEN));
                edits.push(Edit::new(op, u32::try_from(len).expect("at most MAX_LEN")));
                left -= len;
            }
        }
        self.edits = edits;
    }
}

/// What [`find_overlaps`] and its kin keep of each overlap's alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Alignment {
    /// Every column.
    Whole,
    /// Its matches and differences only ([`Overlap::keep_counts`]).
    Counts,
}

/// Two reads, the earlier first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Pair {
    /// The earlier read.
    pub query: ReadIndex,
    /// The later read.
    pub target: ReadIndex,
}

/// What [`find_overlaps_skipping`] found.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Overlaps {
    /// At most one overlap per pair, in query order and then target order.
    pub overlaps: Vec<Overlap>,
    /// The pairs that had chains, none of which could be aligned to the reads' ends, sorted.
    pub unaligned: Vec<Pair>,
}

/// Returns the overlaps between `reads`, at most one per pair (the query is the earlier read),
/// in query order and then target order.
///
/// The result is the same for any number of threads; work runs on the caller's rayon pool.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
///
/// # Panics
///
/// Panics if a read is longer than 2^31 bases.
pub fn find_overlaps(
    reads: &Reads,
    params: &OverlapParams,
) -> Result<Vec<Overlap>, InvalidOverlapParamsError> {
    Ok(find_overlaps_skipping(reads, params, &[])?.overlaps)
}

/// Returns the overlaps between `reads` as [`find_overlaps`] does, each keeping only the counts
/// of its alignment ([`Overlap::keep_counts`]), dropped as each is aligned. The overlaps of
/// noisy reads held whole took most of a run's memory: between *P. kielensis*' 360,000 ONT
/// reads laid out (6.7% from the contigs), 12.9 million overlaps.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
///
/// # Panics
///
/// Panics if a read is longer than 2^31 bases.
pub fn find_overlaps_counted(
    reads: &Reads,
    params: &OverlapParams,
) -> Result<Vec<Overlap>, InvalidOverlapParamsError> {
    Ok(overlaps_between(reads, params, &[], Alignment::Counts)?.overlaps)
}

/// Returns the overlaps between `reads`, as [`find_overlaps`] does, without aligning the pairs
/// in `skip` (sorted), and the pairs none of whose chains could be aligned to the reads' ends.
///
/// Between rounds of correction the reads change at a few bases. A pair whose chains could not
/// be aligned to the reads' ends (an internal match of a shared repeat, or reads of copies far
/// apart) almost always fails again, and those pairs took most of the aligner's work; a pair
/// aligned and then found too divergent is cheap to try again, and after correction is the
/// likeliest to overlap.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
///
/// # Panics
///
/// Panics if `skip` is not sorted, or a read is longer than 2^31 bases.
pub fn find_overlaps_skipping(
    reads: &Reads,
    params: &OverlapParams,
    skip: &[Pair],
) -> Result<Overlaps, InvalidOverlapParamsError> {
    overlaps_between(reads, params, skip, Alignment::Whole)
}

/// [`find_overlaps_skipping`], each overlap keeping `alignment` of its alignment.
fn overlaps_between(
    reads: &Reads,
    params: &OverlapParams,
    skip: &[Pair],
    alignment: Alignment,
) -> Result<Overlaps, InvalidOverlapParamsError> {
    if !(0.0..1.0).contains(&params.max_divergence) {
        return Err(InvalidOverlapParamsError::Divergence(params.max_divergence));
    }
    if let Some(d) = params.max_piece_divergence
        && !(0.0..1.0).contains(&d)
    {
        return Err(InvalidOverlapParamsError::Divergence(d));
    }
    if params.min_seeds == 0 {
        return Err(InvalidOverlapParamsError::MinSeeds);
    }
    assert!(skip.is_sorted(), "the pairs to skip are sorted");
    let sketch = SketchParams::new(params.kmer_len, params.density, 1)?;
    let index = Index::new(reads, sketch);
    log::debug!(
        "overlap: {} reads, {} minimizers indexed",
        reads.len(),
        index.entries.len()
    );
    let per_read: Vec<(Vec<Overlap>, Vec<Pair>)> = reads
        .indices()
        .collect::<Vec<_>>()
        .par_iter()
        .map_init(
            || (Scratch::new(), Sketcher::new(sketch), Vec::new()),
            |(scratch, sketcher, anchors), &query| {
                // Sketched again rather than kept from the index: it is cheap, and every read's
                // minimizers held at once took as much memory as the index.
                sketcher
                    .sketch(reads.get(query), anchors)
                    .expect("read lengths were checked when added");
                let from = skip.partition_point(|p| p.query < query);
                let to = from + skip[from..].partition_point(|p| p.query == query);
                let skipped = &skip[from..to];
                // The first chain to align, of the best first, is the pair's one overlap.
                let mut overlaps: Vec<Overlap> = Vec::new();
                // Targets with a chain tried, and with one aligned to the reads' ends.
                let (mut done, mut tried, mut aligned) =
                    (BTreeSet::new(), BTreeSet::new(), BTreeSet::new());
                let query_len = reads.read_len(query);
                for chain in chains(query_len, anchors, &index, Some(query), params) {
                    if done.contains(&chain.target)
                        || skipped
                            .binary_search_by_key(&chain.target, |p| p.target)
                            .is_ok()
                    {
                        continue;
                    }
                    tried.insert(chain.target);
                    let target = reads.get(chain.target);
                    match scratch.align(
                        target,
                        reads.get(query),
                        query,
                        &chain,
                        params,
                        params.extension,
                    ) {
                        Ok(mut o) => {
                            done.insert(chain.target);
                            aligned.insert(chain.target);
                            if alignment == Alignment::Counts {
                                o.keep_counts();
                            }
                            overlaps.push(o);
                        }
                        Err(NoOverlap::Rejected) => {
                            aligned.insert(chain.target);
                        }
                        Err(NoOverlap::Unaligned) => {}
                    }
                }
                overlaps.sort_by_key(|o| o.target);
                let unaligned = tried
                    .difference(&aligned)
                    .map(|&target| Pair { query, target })
                    .collect();
                (overlaps, unaligned)
            },
        )
        .collect();
    let mut overlaps = Vec::with_capacity(per_read.iter().map(|(o, _)| o.len()).sum());
    let mut unaligned = Vec::with_capacity(per_read.iter().map(|(_, u)| u.len()).sum());
    for (o, u) in per_read {
        overlaps.extend(o);
        unaligned.extend(u);
    }
    Ok(Overlaps {
        overlaps,
        unaligned,
    })
}

/// Returns where each of `queries` lies on `targets` (contigs, say): its best alignment, as an
/// [`Overlap`] whose `query` indexes `queries` and `target` indexes `targets`, or `None`.
///
/// Chains are tried best first, as [`find_overlaps`] tries them, and the first that aligns
/// within [`OverlapParams::max_divergence`] over at least [`OverlapParams::min_overlap`] bases
/// is the mapping. Unlike an overlap, a mapping need not reach the reads' ends: an extension
/// that cannot be aligned (a burst of errors at a read's end, an adapter) ends at its seed. A
/// query that runs off a target's end (a circular contig's start) maps to one side only. With
/// [`OverlapParams::max_chain_splits`], a chain that does not align is split at its widest gap
/// between seeds, and its longer part tried (split again if it does not align), then its
/// shorter, before the next chain: a read
/// across a stretch that differs from the target (another version of a locus) maps up to that
/// stretch, not to a copy of the rest elsewhere.
///
/// The result is the same for any number of threads; work runs on the caller's rayon pool.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
///
/// # Panics
///
/// Panics if a read is longer than 2^31 bases.
pub fn map_reads(
    targets: &Reads,
    queries: &Reads,
    params: &OverlapParams,
) -> Result<Vec<Option<Overlap>>, InvalidOverlapParamsError> {
    if !(0.0..1.0).contains(&params.max_divergence) {
        return Err(InvalidOverlapParamsError::Divergence(params.max_divergence));
    }
    if let Some(d) = params.max_piece_divergence
        && !(0.0..1.0).contains(&d)
    {
        return Err(InvalidOverlapParamsError::Divergence(d));
    }
    if params.min_seeds == 0 {
        return Err(InvalidOverlapParamsError::MinSeeds);
    }
    let sketch = SketchParams::new(params.kmer_len, params.density, 1)?;
    let index = Index::new(targets, sketch);
    Ok(queries
        .indices()
        .collect::<Vec<_>>()
        .par_iter()
        .map_init(
            || (Scratch::new(), Sketcher::new(sketch), Vec::new()),
            |(scratch, sketcher, anchors), &query| {
                let bases = queries.get(query);
                sketcher
                    .sketch(bases, anchors)
                    .expect("read lengths were checked when added");
                chains(queries.read_len(query), anchors, &index, None, params)
                    .into_iter()
                    .find_map(|chain| {
                        let target = targets.get(chain.target);
                        scratch
                            .align(target, bases, query, &chain, params, Extension::Clipped)
                            .ok()
                            .or_else(|| {
                                let splits = params.max_chain_splits;
                                scratch.map_parts(target, bases, query, chain, params, splits)
                            })
                    })
            },
        )
        .collect())
}

/// A chain is split only at a gap between seeds at least this wide, on either read: narrower
/// ones are bursts of errors, which the alignment between them spans.
const MIN_SPLIT_GAP: u32 = 300;

/// One minimizer occurrence in the index, in 16 bytes; ordered by hash, read, start and strand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Entry {
    hash: u64,
    read: ReadIndex,
    /// The start shifted left one bit, and 1 for [`Strand::Reverse`].
    place: u32,
}

impl Entry {
    fn new(hash: u64, read: ReadIndex, start: u32, strand: Strand) -> Self {
        assert!(
            start <= u32::MAX >> 1,
            "a minimizer at base {start} is past the index's reach"
        );
        Self {
            hash,
            read,
            place: start << 1 | u32::from(strand.is_reverse()),
        }
    }

    fn start(self) -> u32 {
        self.place >> 1
    }

    fn strand(self) -> Strand {
        if self.place & 1 == 1 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
}

/// Every read's minimizers, sorted by hash.
struct Index {
    entries: Vec<Entry>,
    /// Where each bucket of hashes starts in `entries`, and one past the last: bucket `b`
    /// holds the hashes `h` with `h >> shift == b`. A lookup searches its bucket only; on the
    /// fly's gathered reads the search of all `entries` took 13% of the overlaps' time.
    starts: Vec<usize>,
    /// The bits below a bucket's.
    shift: u32,
}

/// The index aims at this many entries per bucket (between 4 and 16): a short search, and a
/// table of starts an eighth of the entries' count.
const ENTRIES_PER_BUCKET: usize = 8;

/// Returns the bucket of an indexed `hash` (see [`Index::starts`]).
///
/// # Panics
///
/// Panics if the bucket does not fit a `usize`, which cannot happen: there are fewer buckets
/// than entries.
fn bucket(hash: u64, shift: u32) -> usize {
    usize::try_from(hash >> shift).expect("there are fewer buckets than entries")
}

impl Index {
    fn new(reads: &Reads, sketch: SketchParams) -> Self {
        let indices: Vec<ReadIndex> = reads.indices().collect();
        let sketch_read = |sketcher: &mut Sketcher, anchors: &mut Vec<Anchor>, read: ReadIndex| {
            sketcher
                .sketch(reads.get(read), anchors)
                .expect("read lengths were checked when added");
        };
        // Sketched twice, to count and then to fill: every read's entries held apart and then
        // copied into one vector took twice the index's memory.
        let counts: Vec<usize> = indices
            .par_iter()
            .map_init(
                || (Sketcher::new(sketch), Vec::new()),
                |(sketcher, anchors), &read| {
                    sketch_read(sketcher, anchors, read);
                    anchors.len()
                },
            )
            .collect();
        let placeholder = Entry {
            hash: 0,
            read: ReadIndex(0),
            place: 0,
        };
        let mut entries = vec![placeholder; counts.iter().sum()];
        let mut slots: Vec<&mut [Entry]> = Vec::with_capacity(counts.len());
        let mut rest = &mut entries[..];
        for &count in &counts {
            let (slot, tail) = rest.split_at_mut(count);
            slots.push(slot);
            rest = tail;
        }
        indices.par_iter().zip(slots).for_each_init(
            || (Sketcher::new(sketch), Vec::new()),
            |(sketcher, anchors), (&read, slot)| {
                sketch_read(sketcher, anchors, read);
                debug_assert_eq!(slot.len(), anchors.len(), "sketching is deterministic");
                for (entry, a) in slot.iter_mut().zip(anchors.iter()) {
                    *entry = Entry::new(a.hash, read, a.raw_start, a.strand);
                }
            },
        );
        entries.par_sort_unstable();
        // Universe minimizers' hashes are spread evenly below the density's threshold: as many
        // buckets as fit the largest hash's bits, about ENTRIES_PER_BUCKET entries each.
        let bits = (entries.len() / ENTRIES_PER_BUCKET)
            .max(1)
            .next_power_of_two()
            .trailing_zeros();
        let span = entries
            .last()
            .map_or(0, |e| u64::BITS - e.hash.leading_zeros());
        let shift = span.saturating_sub(bits);
        let buckets = entries.last().map_or(0, |e| bucket(e.hash, shift)) + 1;
        let mut starts = vec![0; buckets + 1];
        for e in &entries {
            starts[bucket(e.hash, shift) + 1] += 1;
        }
        for b in 1..starts.len() {
            starts[b] += starts[b - 1];
        }
        Self {
            entries,
            starts,
            shift,
        }
    }

    /// Returns the occurrences of `hash`.
    fn lookup(&self, hash: u64) -> &[Entry] {
        // A hash past the largest indexed lands past the last bucket.
        let b = usize::try_from(hash >> self.shift).unwrap_or(usize::MAX - 1);
        let (Some(&start), Some(&end)) = (
            self.starts.get(b),
            b.checked_add(1).and_then(|next| self.starts.get(next)),
        ) else {
            return &[];
        };
        let in_bucket = &self.entries[start..end];
        let from = in_bucket.partition_point(|e| e.hash < hash);
        let to = from + in_bucket[from..].partition_point(|e| e.hash == hash);
        &in_bucket[from..to]
    }
}

/// A seed: a minimizer at `target` in the target and `query` in the query as aligned (the
/// reverse complement's coordinates on [`Strand::Reverse`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Seed {
    target: u32,
    query: u32,
}

/// A chain of seeds between a query and one target, on one relative strand.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Chain {
    target: ReadIndex,
    strand: Strand,
    seeds: Vec<Seed>,
}

/// Seeds further apart than this on either read are not chained.
const MAX_SEED_GAP: u32 = 20_000;

/// Each seed looks this many seeds back for its predecessor.
const MAX_PREDECESSORS: usize = 64;

/// Two chained seeds may leave their diagonal by this many bases (the default
/// [`OverlapParams::chain_slack`]), plus a twentieth of the gap between them (indels and
/// homopolymer-length differences).
pub const DIAGONAL_SLACK: u32 = 50;

/// Chains are taken from one pair's seeds on one strand at most this many times (the default
/// [`OverlapParams::chains_per_pair`]), each from the seeds the ones before left: the best
/// chain may be a repeat inside both reads, whose extension fails, while the overlap is the
/// next.
pub const CHAINS_PER_PAIR: usize = 3;

/// A seed of a query on a target read and strand, packed into one number whose order is that
/// of (target, strand, seed): the target's index in bits 65 and up, the strand in bit 64
/// (forward 0), the seed's target position in bits 32 to 63 and its query position below.
/// Sorting these took half the time sorting the tuples did.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Hit(u128);

impl Hit {
    fn new(target: ReadIndex, strand: Strand, seed: Seed) -> Self {
        Self(
            u128::from(target.0) << 65
                | u128::from(strand.is_reverse()) << 64
                | u128::from(seed.target) << 32
                | u128::from(seed.query),
        )
    }

    /// The target and strand, as one number in the same order.
    fn pair(self) -> u64 {
        u64::try_from(self.0 >> 64).expect("a read index and a strand fit in 64 bits")
    }

    fn target(self) -> ReadIndex {
        ReadIndex(u32::try_from(self.pair() >> 1).expect("a read index fits in u32"))
    }

    fn strand(self) -> Strand {
        if self.pair() & 1 == 1 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }

    fn seed(self) -> Seed {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "the seed's two 32-bit halves"
        )]
        Seed {
            target: (self.0 >> 32) as u32,
            query: self.0 as u32,
        }
    }
}

/// Returns the best chains of a query of `query_len` bases, whose minimizers are `anchors`, to
/// the indexed reads (only those after `after`, when given), at most
/// [`OverlapParams::max_candidates`], most seeds first, then in target order.
fn chains(
    query_len: u32,
    anchors: &[Anchor],
    index: &Index,
    after: Option<ReadIndex>,
    params: &OverlapParams,
) -> Vec<Chain> {
    let mut hits: Vec<Hit> = Vec::new();
    for anchor in anchors {
        let occurrences = index.lookup(anchor.hash);
        if occurrences.len() > usize::try_from(params.max_occurrences).unwrap_or(usize::MAX) {
            continue;
        }
        for e in occurrences
            .iter()
            .filter(|e| after.is_none_or(|query| e.read > query))
        {
            let strand = anchor.strand.compose(e.strand());
            let query_pos = if strand.is_reverse() {
                query_len - anchor.raw_end
            } else {
                anchor.raw_start
            };
            let seed = Seed {
                target: e.start(),
                query: query_pos,
            };
            hits.push(Hit::new(e.read, strand, seed));
        }
    }
    hits.sort_unstable();
    // A pair's chains hold at most its seeds: pairs are chained most seeds first, until the
    // candidates kept have more seeds each than the next pair has in all, as no chain of it or
    // of a later pair could then be kept (an equal one could, on its target).
    let mut groups: Vec<&[Hit]> = hits.chunk_by(|a, b| a.pair() == b.pair()).collect();
    groups.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].pair().cmp(&b[0].pair())));
    let mut chains: Vec<Chain> = Vec::new();
    // The seed counts of the best `max_candidates` chains so far, least on top.
    let mut best: BinaryHeap<Reverse<usize>> = BinaryHeap::new();
    for group in groups {
        if best.len() >= params.max_candidates
            && best
                .peek()
                .is_some_and(|&Reverse(least)| least > group.len())
        {
            break;
        }
        let mut left: Vec<Seed> = group.iter().map(|h| h.seed()).collect();
        let (target, strand) = (group[0].target(), group[0].strand());
        let mut found = 0;
        while found < params.chains_per_pair && left.len() >= params.min_seeds {
            let seeds = best_chain(&left, params.chain_slack);
            if seeds.len() < params.min_seeds {
                break;
            }
            // Both are sorted, and the chain is a subsequence of what is left.
            let mut used = seeds.iter().peekable();
            left.retain(|s| {
                if used.peek() == Some(&s) {
                    used.next();
                    false
                } else {
                    true
                }
            });
            best.push(Reverse(seeds.len()));
            if best.len() > params.max_candidates {
                best.pop();
            }
            chains.push(Chain {
                target,
                strand,
                seeds,
            });
            found += 1;
        }
    }
    chains.sort_by(|a, b| {
        b.seeds
            .len()
            .cmp(&a.seeds.len())
            .then(a.target.cmp(&b.target))
            .then(a.strand.cmp(&b.strand))
    });
    chains.truncate(params.max_candidates);
    chains
}

/// Returns the longest co-linear chain of `seeds` (sorted), ties to the earliest end.
fn best_chain(seeds: &[Seed], slack: u32) -> Vec<Seed> {
    let mut score = vec![1_usize; seeds.len()];
    // The best score of any seed up to each one: once the predecessors left to look at cannot
    // tie the best found, the look back stops (on a clean chain, after two or three seeds
    // rather than 64), with the same result.
    let mut best_before = vec![1_usize; seeds.len()];
    let mut previous: Vec<Option<usize>> = vec![None; seeds.len()];
    for i in 0..seeds.len() {
        let s = seeds[i];
        let mut best: Option<(usize, u32, usize)> = None; // (score, drift, j)
        for j in (i.saturating_sub(MAX_PREDECESSORS)..i).rev() {
            if best.is_some_and(|(sc, _, _)| best_before[j] + 1 < sc) {
                break;
            }
            let p = seeds[j];
            if p.target >= s.target || p.query >= s.query {
                continue;
            }
            let (dt, dq) = (s.target - p.target, s.query - p.query);
            let gap = dt.max(dq);
            let drift = dt.abs_diff(dq);
            if gap > MAX_SEED_GAP || drift > slack.saturating_add(gap / 20) {
                continue;
            }
            let candidate = (score[j] + 1, drift, j);
            if best.is_none_or(|(sc, dr, _)| candidate.0 > sc || (candidate.0 == sc && drift < dr))
            {
                best = Some(candidate);
            }
        }
        if let Some((sc, _, j)) = best {
            score[i] = sc;
            previous[i] = Some(j);
        }
        best_before[i] = if i == 0 {
            score[i]
        } else {
            best_before[i - 1].max(score[i])
        };
    }
    let Some(end) = (0..seeds.len()).max_by(|&a, &b| score[a].cmp(&score[b]).then(b.cmp(&a)))
    else {
        return Vec::new();
    };
    let mut chain = vec![seeds[end]];
    let mut at = end;
    while let Some(p) = previous[at] {
        chain.push(seeds[p]);
        at = p;
    }
    chain.reverse();
    chain
}

/// Seeds closer than this to the last one aligned are skipped: fewer, longer pieces.
const MIN_PIECE: u32 = 200;

/// A piece or an extension may hold this many edits beyond its share of the divergence (and,
/// for a piece, its seeds' drift from one diagonal): a short tandem repeat's length can differ
/// by tens of bases between two reads of one copy. The overlap as a whole must still be within
/// [`OverlapParams::max_divergence`].
const EXTRA_EDITS: usize = 32;

/// No piece or extension is aligned with more edits than this: the aligner's memory grows with
/// the square of the cap.
const MAX_EDITS: usize = 2_000;

/// An extension is given up once its last 128 edits moved it fewer than 512 bases along each read
/// (1,024 together): it has run into sequence the other read lacks, past a shared repeat's end or
/// into another copy's. Unrelated sequence moves it about 520 bases (together) in 128 edits; an
/// overlap's extension moves it hundreds of bases an edit, and crosses an indel of up to some 120
/// bases. Without this, every failing extension ran to its cap, and on the fly's repeat reads
/// they took 97% of the aligner's work. A 64-edit window (512 bases) saves another 38% of the
/// CPU but drops reads with a burst of errors (one read with some 100 inserted bases in 200 lost
/// all its overlaps, and stayed uncorrected).
const EXTENSION_STALL: Stall = Stall {
    edits: 128,
    advance: 1_024,
};

/// Why a chain is no overlap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum NoOverlap {
    /// An extension or a piece between seeds could not be aligned (the cap, or a stall).
    Unaligned,
    /// Aligned to the reads' ends, but too divergent or too short.
    Rejected,
}

/// Per-thread buffers for aligning chains.
struct Scratch {
    aligner: DiffAligner,
    oriented: Vec<u8>,
    target_rev: Vec<u8>,
    query_rev: Vec<u8>,
    left: Vec<Edit>,
    right: Vec<Edit>,
    piece: Vec<Edit>,
}

impl Scratch {
    fn new() -> Self {
        Self {
            aligner: DiffAligner::new(0),
            oriented: Vec::new(),
            target_rev: Vec::new(),
            query_rev: Vec::new(),
            left: Vec::new(),
            right: Vec::new(),
            piece: Vec::new(),
        }
    }

    /// Returns a cap on edits over `len` columns whose ends lie `drift` bases off one diagonal,
    /// at `divergence`.
    fn cap(len: usize, drift: usize, divergence: f64) -> usize {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a non-negative share of a read length"
        )]
        let share = (len as f64 * divergence).ceil() as usize;
        (share + drift + EXTRA_EDITS).min(MAX_EDITS)
    }

    /// Splits `chain` (one that does not align) at its widest gap between seeds, if at least
    /// [`MIN_SPLIT_GAP`] wide, and aligns its longer part with clipped extensions, then its
    /// shorter; a part that does not align is split in turn, with `splits - 1` splits left.
    fn map_parts(
        &mut self,
        target_seq: &[u8],
        query_fwd: &[u8],
        query: ReadIndex,
        chain: Chain,
        params: &OverlapParams,
        splits: u32,
    ) -> Option<Overlap> {
        if splits == 0 {
            return None;
        }
        let gap = |w: &[Seed]| (w[1].target - w[0].target).max(w[1].query - w[0].query);
        let (at, widest) = chain
            .seeds
            .windows(2)
            .enumerate()
            .map(|(i, w)| (i + 1, gap(w)))
            .max_by_key(|&(i, g)| (g, Reverse(i)))?;
        if widest < MIN_SPLIT_GAP {
            return None;
        }
        let Chain {
            target,
            strand,
            mut seeds,
        } = chain;
        let right = seeds.split_off(at);
        let span = |s: &[Seed]| s[s.len() - 1].target - s[0].target;
        let (first, second) = if span(&seeds) >= span(&right) {
            (seeds, right)
        } else {
            (right, seeds)
        };
        [first, second].into_iter().find_map(|seeds| {
            let part = Chain {
                target,
                strand,
                seeds,
            };
            self.align(
                target_seq,
                query_fwd,
                query,
                &part,
                params,
                Extension::Clipped,
            )
            .ok()
            .or_else(|| self.map_parts(target_seq, query_fwd, query, part, params, splits - 1))
        })
    }

    /// Aligns `chain` of `query` (whose bases are `query_fwd`) to `target_seq` from its seeds
    /// out to the reads' ends (or, with [`Extension::Clipped`], as far as the extensions get),
    /// if it is an overlap.
    fn align(
        &mut self,
        target_seq: &[u8],
        query_fwd: &[u8],
        query: ReadIndex,
        chain: &Chain,
        params: &OverlapParams,
        extension: Extension,
    ) -> Result<Overlap, NoOverlap> {
        self.oriented.clear();
        if chain.strand.is_reverse() {
            picoasm_seq::reverse_complement_into(query_fwd, &mut self.oriented);
        } else {
            self.oriented.extend_from_slice(query_fwd);
        }
        let query_seq = &self.oriented[..];
        self.aligner.set_max_lag(params.max_lag.map(pos));
        let first = chain.seeds[0];
        let last = *chain.seeds.last().expect("chains have seeds");

        // Left: from the first seed back to whichever read starts first.
        let (target_before, query_before) = (pos(first.target), pos(first.query));
        let max_edits = Self::cap(target_before.min(query_before), 0, params.max_divergence);
        // With `max_edits` edits the extension covers at most that many target bases more than
        // query bases, so a longer target (a contig, when reads are mapped onto contigs) is
        // reversed only that far, plus one base so that its cut end is never reached: the
        // alignment is the same, without copying megabases per read.
        let reach = target_before.min(query_before + max_edits + 1);
        self.target_rev.clear();
        self.target_rev.extend(
            target_seq[target_before - reach..target_before]
                .iter()
                .rev(),
        );
        self.query_rev.clear();
        self.query_rev
            .extend(query_seq[..query_before].iter().rev());
        self.aligner.set_max_edits(max_edits);
        self.left.clear();
        self.aligner.set_stall(Some(EXTENSION_STALL));
        let left = match (
            self.aligner
                .extend(&self.target_rev, &self.query_rev, &mut self.left),
            extension,
        ) {
            (Some((_, extent)), _) => extent,
            (None, Extension::Clipped) => Extent::default(),
            (None, Extension::Full) => return Err(NoOverlap::Unaligned),
        };

        // Right: from the last seed on to whichever read ends first. Before the pieces
        // between seeds, which hold most of the bases: a chain that is no overlap (a repeat's
        // copy inside both reads) most often fails at an end.
        let t = &target_seq[pos(last.target)..];
        let q = &query_seq[pos(last.query)..];
        self.aligner
            .set_max_edits(Self::cap(t.len().min(q.len()), 0, params.max_divergence));
        self.right.clear();
        let right = match (self.aligner.extend(t, q, &mut self.right), extension) {
            (Some((_, extent)), _) => extent,
            (None, Extension::Clipped) => Extent::default(),
            (None, Extension::Full) => return Err(NoOverlap::Unaligned),
        };

        // Between seeds the reads are known to share sequence: only the cap ends a piece.
        self.aligner.set_stall(None);
        let mut edits: Vec<Edit> = Vec::new();
        for e in self.left.iter().rev() {
            push(&mut edits, e.op(), e.len());
        }
        let target_start = first.target - u32::try_from(left.target).expect("within a read");
        let query_start = first.query - u32::try_from(left.query).expect("within a read");

        // Between seeds.
        let mut at = first;
        for &seed in &chain.seeds[1..] {
            if seed.target - at.target < MIN_PIECE && seed != last {
                continue;
            }
            let t = &target_seq[pos(at.target)..pos(seed.target)];
            let q = &query_seq[pos(at.query)..pos(seed.query)];
            self.aligner.set_max_edits(Self::cap(
                t.len().max(q.len()),
                t.len().abs_diff(q.len()),
                params.max_piece_divergence.unwrap_or(params.max_divergence),
            ));
            self.piece.clear();
            self.aligner
                .global(t, q, &mut self.piece)
                .ok_or(NoOverlap::Unaligned)?;
            for e in &self.piece {
                push(&mut edits, e.op(), e.len());
            }
            at = seed;
        }

        for e in &self.right {
            push(&mut edits, e.op(), e.len());
        }
        let target_end = last.target + u32::try_from(right.target).expect("within a read");
        let oriented_end = last.query + u32::try_from(right.query).expect("within a read");

        let query_len = u32::try_from(query_fwd.len()).expect("read lengths fit in u32");
        let (query_start, query_end) = if chain.strand.is_reverse() {
            (query_len - oriented_end, query_len - query_start)
        } else {
            (query_start, oriented_end)
        };
        let overlap = Overlap {
            query,
            target: chain.target,
            strand: chain.strand,
            query_start,
            query_end,
            target_start,
            target_end,
            edits,
        };
        #[expect(
            clippy::cast_precision_loss,
            reason = "a share of a read length's columns"
        )]
        let divergence = overlap.differences() as f64 / overlap.columns().max(1) as f64;
        let long_enough = query_end - query_start >= params.min_overlap
            && target_end - target_start >= params.min_overlap;
        if long_enough && divergence <= params.max_divergence {
            let mut overlap = overlap;
            overlap.edits.shrink_to_fit();
            Ok(overlap)
        } else {
            Err(NoOverlap::Rejected)
        }
    }
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// [`best_chain`] looking back over all [`MAX_PREDECESSORS`] seeds every time.
    fn best_chain_looking_back_fully(seeds: &[Seed], slack: u32) -> Vec<Seed> {
        let mut score = vec![1_usize; seeds.len()];
        let mut previous: Vec<Option<usize>> = vec![None; seeds.len()];
        for i in 0..seeds.len() {
            let s = seeds[i];
            let mut best: Option<(usize, u32, usize)> = None;
            for j in (i.saturating_sub(MAX_PREDECESSORS)..i).rev() {
                let p = seeds[j];
                if p.target >= s.target || p.query >= s.query {
                    continue;
                }
                let (dt, dq) = (s.target - p.target, s.query - p.query);
                let gap = dt.max(dq);
                let drift = dt.abs_diff(dq);
                if gap > MAX_SEED_GAP || drift > slack.saturating_add(gap / 20) {
                    continue;
                }
                let candidate = (score[j] + 1, drift, j);
                if best
                    .is_none_or(|(sc, dr, _)| candidate.0 > sc || (candidate.0 == sc && drift < dr))
                {
                    best = Some(candidate);
                }
            }
            if let Some((sc, _, j)) = best {
                score[i] = sc;
                previous[i] = Some(j);
            }
        }
        let Some(end) = (0..seeds.len()).max_by(|&a, &b| score[a].cmp(&score[b]).then(b.cmp(&a)))
        else {
            return Vec::new();
        };
        let mut chain = vec![seeds[end]];
        let mut at = end;
        while let Some(p) = previous[at] {
            chain.push(seeds[p]);
            at = p;
        }
        chain.reverse();
        chain
    }

    #[test]
    fn hits_sort_as_their_target_strand_and_seed_do() {
        let values = [0, 1, 2, 1 << 31, u32::MAX - 1, u32::MAX];
        let mut tuples = Vec::new();
        for &read in &values {
            for strand in [Strand::Forward, Strand::Reverse] {
                for &target in &values {
                    for &query in &values {
                        tuples.push((ReadIndex(read), strand, Seed { target, query }));
                    }
                }
            }
        }
        tuples.reverse();
        let mut hits: Vec<Hit> = tuples
            .iter()
            .map(|&(r, s, seed)| Hit::new(r, s, seed))
            .collect();
        tuples.sort_unstable();
        hits.sort_unstable();
        let unpacked: Vec<(ReadIndex, Strand, Seed)> = hits
            .iter()
            .map(|h| (h.target(), h.strand(), h.seed()))
            .collect();
        assert_eq!(unpacked, tuples);
    }

    #[test]
    fn the_best_chain_is_the_one_a_full_look_back_finds() {
        let mut state = 5_u64;
        let mut next = move |n: u32| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            u32::try_from((state >> 33) % u64::from(n)).expect("below n")
        };
        for case in 0..500 {
            // Chains along a few diagonals, with drift, ties and repeats' stray seeds, sorted
            // as the pairs' seeds are.
            let mut seeds: Vec<Seed> = Vec::new();
            let diagonals: Vec<u32> = (0..1 + next(4)).map(|_| next(3_000)).collect();
            for _ in 0..next(300) {
                let target = next(30_000);
                let seed = if next(5) == 0 {
                    Seed {
                        target,
                        query: next(30_000),
                    }
                } else {
                    let diagonal =
                        diagonals[usize::try_from(next(4)).expect("small") % diagonals.len()];
                    Seed {
                        target,
                        query: (target + diagonal + next(40)).saturating_sub(20),
                    }
                };
                seeds.push(seed);
            }
            seeds.sort_unstable();
            seeds.dedup();
            let slack = [0, 10, DIAGONAL_SLACK][usize::try_from(next(3)).expect("small")];
            assert_eq!(
                best_chain(&seeds, slack),
                best_chain_looking_back_fully(&seeds, slack),
                "case {case}"
            );
        }
    }

    /// A random genome from a fixed linear congruential sequence.
    fn genome(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
            })
            .collect()
    }

    fn params() -> OverlapParams {
        OverlapParams {
            density: 0.2,
            min_overlap: 500,
            ..OverlapParams::default()
        }
    }

    #[test]
    fn a_lookup_finds_exactly_the_entries_of_its_hash() {
        // Repeated pieces, so that some hashes are held many times.
        let g = genome(40_000, 11);
        let pieces: Vec<Vec<u8>> = (0..30)
            .map(|i| g[(i * 997) % 30_000..(i * 997) % 30_000 + 8_000].to_vec())
            .collect();
        let reads: Reads = pieces.iter().map(Vec::as_slice).collect();
        let p = params();
        let index = Index::new(&reads, SketchParams::new(p.kmer_len, p.density, 1).unwrap());
        assert!(index.entries.len() > 1_000);
        let scan = |hash: u64| -> Vec<Entry> {
            index
                .entries
                .iter()
                .filter(|e| e.hash == hash)
                .copied()
                .collect()
        };
        for e in index.entries.iter().step_by(7) {
            assert_eq!(index.lookup(e.hash), scan(e.hash).as_slice());
        }
        let largest = index.entries.last().unwrap().hash;
        for absent in [0, 1, largest + 1, largest / 2 + 1, u64::MAX] {
            assert_eq!(index.lookup(absent), scan(absent).as_slice());
        }
    }

    #[test]
    fn the_candidates_are_the_best_of_all_chains() {
        // Reads of a 30 kb stretch at many offsets, so each has many targets with chains of
        // many sizes, some of one size: the candidates chained with a limit are those of all
        // chains, best first, cut to the limit.
        let g = genome(40_000, 11);
        let pieces: Vec<Vec<u8>> = (0..60)
            .map(|i| g[(i * 613) % 30_000..(i * 613) % 30_000 + 8_000].to_vec())
            .collect();
        let reads: Reads = pieces.iter().map(Vec::as_slice).collect();
        let p = params();
        let index = Index::new(&reads, SketchParams::new(p.kmer_len, p.density, 1).unwrap());
        let mut sketcher = Sketcher::new(SketchParams::new(p.kmer_len, p.density, 1).unwrap());
        let unlimited = OverlapParams {
            max_candidates: usize::MAX,
            ..p
        };
        for query in reads.indices().step_by(5) {
            let mut anchors = Vec::new();
            sketcher.sketch(reads.get(query), &mut anchors).unwrap();
            let len = reads.read_len(query);
            let all = chains(len, &anchors, &index, None, &unlimited);
            assert!(all.len() > 7, "{} chains", all.len());
            for limit in [1, 2, 3, 5, 7] {
                let limited = OverlapParams {
                    max_candidates: limit,
                    ..p
                };
                let cut = chains(len, &anchors, &index, None, &limited);
                assert_eq!(cut, all[..limit], "limit {limit}");
            }
        }
    }

    #[test]
    fn counted_overlaps_are_the_whole_ones_with_their_counts_only() {
        let g = genome(40_000, 12);
        // Reads of 8 kb with a base changed every 97 and one left out every 211, from where
        // each starts: their overlaps differ in substitutions and indels.
        let pieces: Vec<Vec<u8>> = (0..40)
            .map(|i| {
                let from = (i * 787) % 30_000;
                g[from..from + 8_000]
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| (from + j) % 211 != 0)
                    .map(|(j, &b)| match ((from + j) % 97, b) {
                        (0, b'A') => b'C',
                        (0, _) => b'A',
                        _ => b,
                    })
                    .collect()
            })
            .collect();
        let reads: Reads = pieces.iter().map(Vec::as_slice).collect();
        let whole = find_overlaps(&reads, &params()).expect("valid");
        let counted = find_overlaps_counted(&reads, &params()).expect("valid");
        assert!(whole.len() > 40, "{} overlaps", whole.len());
        assert!(whole.iter().any(|o| o.edits.len() > 2));
        assert_eq!(whole.len(), counted.len());
        for (w, c) in whole.iter().zip(&counted) {
            assert!(c.edits.len() <= 2);
            assert_eq!(
                (w.columns(), w.differences()),
                (c.columns(), c.differences())
            );
            let mut kept = w.clone();
            kept.keep_counts();
            assert_eq!(&kept, c);
        }
    }

    #[test]
    fn overlapping_reads_overlap_on_both_strands() {
        let g = genome(6_000, 7);
        let a = g[0..3_000].to_vec();
        let b = g[2_000..5_000].to_vec();
        let mut c = picoasm_seq::reverse_complement(&g[2_500..6_000]);
        c[1_000] = if c[1_000] == b'A' { b'C' } else { b'A' };
        let reads: Reads = [&a[..], &b[..], &c[..]].into_iter().collect();
        let overlaps = find_overlaps(&reads, &params()).expect("valid settings");
        let pairs: Vec<_> = overlaps
            .iter()
            .map(|o| (o.query.0, o.target.0, o.strand))
            .collect();
        assert_eq!(
            pairs,
            [
                (0, 1, Strand::Forward),
                (0, 2, Strand::Reverse),
                (1, 2, Strand::Reverse)
            ]
        );
        let ab = &overlaps[0];
        assert_eq!((ab.query_start, ab.query_end), (2_000, 3_000));
        assert_eq!((ab.target_start, ab.target_end), (0, 1_000));
        assert_eq!(ab.differences(), 0);
        // b (2,000..5,000) against c (the reverse of 2,500..6,000): b's 2,500..5,000 is c's
        // last 2,500 bases, reverse-complemented, with the one substitution.
        let bc = &overlaps[2];
        assert_eq!((bc.query_start, bc.query_end), (500, 3_000));
        assert_eq!((bc.target_start, bc.target_end), (1_000, 3_500));
        assert_eq!(bc.differences(), 1);
    }

    #[test]
    fn a_shared_repeat_inside_both_reads_is_no_overlap() {
        let repeat = genome(2_000, 11);
        let mut a = genome(2_000, 12);
        a.extend_from_slice(&repeat);
        a.extend(genome(2_000, 13));
        let mut b = genome(2_000, 14);
        b.extend_from_slice(&repeat);
        b.extend(genome(2_000, 15));
        let reads: Reads = [&a[..], &b[..]].into_iter().collect();
        assert!(find_overlaps(&reads, &params()).expect("valid").is_empty());
    }

    #[test]
    fn pairs_that_cannot_be_aligned_are_reported_and_skipped() {
        // a and b share only a repeat inside both; b and c overlap.
        let repeat = genome(2_000, 51);
        let mut a = genome(2_000, 52);
        a.extend_from_slice(&repeat);
        a.extend(genome(2_000, 53));
        let tail = genome(2_000, 54);
        let mut b = genome(2_000, 55);
        b.extend_from_slice(&repeat);
        b.extend_from_slice(&tail);
        let mut c = tail.clone();
        c.extend(genome(2_000, 56));
        let reads: Reads = [&a[..], &b[..], &c[..]].into_iter().collect();
        let pair = |query, target| Pair {
            query: ReadIndex(query),
            target: ReadIndex(target),
        };
        let found = find_overlaps_skipping(&reads, &params(), &[]).expect("valid");
        assert_eq!(found.unaligned, [pair(0, 1)]);
        assert_eq!(found.overlaps.len(), 1);
        assert_eq!(
            (found.overlaps[0].query, found.overlaps[0].target),
            (ReadIndex(1), ReadIndex(2))
        );

        // Skipped pairs are not aligned, so neither overlap nor are reported.
        let again = find_overlaps_skipping(&reads, &params(), &found.unaligned).expect("valid");
        assert!(again.unaligned.is_empty());
        assert_eq!(again.overlaps, found.overlaps);
        let none =
            find_overlaps_skipping(&reads, &params(), &[pair(0, 1), pair(1, 2)]).expect("valid");
        assert_eq!(none, Overlaps::default());
    }

    #[test]
    fn frequent_minimizers_seed_nothing() {
        let g = genome(3_000, 21);
        let copies: Vec<Vec<u8>> = (0..5).map(|_| g.clone()).collect();
        let reads: Reads = copies.iter().map(Vec::as_slice).collect();
        let few = OverlapParams {
            max_occurrences: 4,
            ..params()
        };
        assert!(find_overlaps(&reads, &few).expect("valid").is_empty());
        assert_eq!(find_overlaps(&reads, &params()).expect("valid").len(), 10);
    }

    #[test]
    fn a_tandem_repeat_length_difference_keeps_the_overlap() {
        // One read has 25 more bases of a CA repeat than the other.
        let left = genome(5_000, 31);
        let right = genome(5_000, 32);
        let mut a = left.clone();
        a.extend(b"CA".repeat(20));
        a.extend_from_slice(&right);
        let mut b = left[1_000..].to_vec();
        b.extend(b"CA".repeat(20));
        b.extend(b"CACACACACACACACACACACACAC");
        b.extend_from_slice(&right);
        let reads: Reads = [&a[..], &b[..]].into_iter().collect();
        let overlaps = find_overlaps(&reads, &params()).expect("valid");
        assert_eq!(overlaps.len(), 1);
        assert_eq!(overlaps[0].differences(), 25);
    }

    #[test]
    fn a_duplicate_inside_both_reads_does_not_hide_their_overlap() {
        // A 2.5 kb duplicate in the middle of each read, and a 1 kb overlap at their ends.
        let dup = genome(2_500, 41);
        let shared = genome(1_000, 42);
        let mut a = genome(1_000, 43);
        a.extend_from_slice(&dup);
        a.extend(genome(1_000, 44));
        a.extend_from_slice(&shared);
        let mut b = shared.clone();
        b.extend(genome(1_000, 45));
        b.extend_from_slice(&dup);
        b.extend(genome(1_000, 46));
        let reads: Reads = [&a[..], &b[..]].into_iter().collect();
        let overlaps = find_overlaps(&reads, &params()).expect("valid");
        assert_eq!(overlaps.len(), 1);
        let o = &overlaps[0];
        assert_eq!((o.query_start, o.query_end), (4_500, 5_500));
        assert_eq!((o.target_start, o.target_end), (0, 1_000));
    }

    #[test]
    fn reads_map_onto_contigs_with_their_ends_clipped() {
        let contig = genome(20_000, 61);
        let targets: Reads = [&contig[..]].into_iter().collect();
        // Inside the contig on the forward strand; on the reverse strand with 300 unrelated
        // bases at its end; and unrelated.
        let a = contig[2_000..7_000].to_vec();
        let mut b = picoasm_seq::reverse_complement(&contig[10_000..16_000]);
        b.extend(genome(300, 62));
        let c = genome(4_000, 63);
        let queries: Reads = [&a[..], &b[..], &c[..]].into_iter().collect();
        let mapped = map_reads(&targets, &queries, &params()).expect("valid settings");
        assert_eq!(mapped.len(), 3);
        let a_on = mapped[0].as_ref().expect("a maps");
        assert_eq!(a_on.strand, Strand::Forward);
        assert_eq!((a_on.target_start, a_on.target_end), (2_000, 7_000));
        assert_eq!((a_on.query_start, a_on.query_end), (0, 5_000));
        let b_on = mapped[1].as_ref().expect("b maps");
        assert_eq!(b_on.strand, Strand::Reverse);
        // The unrelated tail (the query's last 300 bases, before the contig's 10,000 on its
        // reverse strand) is clipped, give or take the last seed's distance from it.
        assert_eq!((b_on.query_start, b_on.target_end), (0, 16_000));
        assert!((5_800..=6_010).contains(&b_on.query_end), "{b_on:?}");
        assert!((9_990..=10_200).contains(&b_on.target_start), "{b_on:?}");
        assert!(mapped[2].is_none());
    }

    #[test]
    fn a_read_across_a_stretch_unlike_the_contig_maps_up_to_it() {
        // The read has 2 kb of other bases where the contig has 2 kb of its own: its chain
        // across them does not align, and its longer part (the 6 kb before) does.
        let contig = genome(20_000, 71);
        let targets: Reads = [&contig[..]].into_iter().collect();
        let read = [
            &contig[2_000..8_000],
            &genome(2_000, 72),
            &contig[10_000..13_000],
        ]
        .concat();
        let queries: Reads = [&read[..]].into_iter().collect();
        let mapped = map_reads(&targets, &queries, &params()).expect("valid settings");
        assert!(mapped[0].is_none(), "no chain aligns whole");
        let split = OverlapParams {
            max_chain_splits: 2,
            ..params()
        };
        let mapped = map_reads(&targets, &queries, &split).expect("valid settings");
        let on = mapped[0].as_ref().expect("the read maps");
        assert_eq!((on.target_start, on.query_start), (2_000, 0));
        assert!((7_800..=8_010).contains(&on.target_end), "{on:?}");
    }

    #[test]
    fn bad_settings_are_refused() {
        let reads = Reads::new();
        let bad = OverlapParams {
            max_divergence: 1.5,
            ..OverlapParams::default()
        };
        assert!(find_overlaps(&reads, &bad).is_err());
        let bad = OverlapParams {
            kmer_len: 40,
            ..OverlapParams::default()
        };
        assert!(find_overlaps(&reads, &bad).is_err());
    }
}
