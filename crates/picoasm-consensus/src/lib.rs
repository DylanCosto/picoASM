//! Per-window base-level consensus for the picoasm bacterial assembler.
//!
//! The assembler cuts every unitig into short windows between consecutive anchors and collects,
//! for each window, the segments of the reads that cross it, already oriented to the unitig's
//! strand (a [`Window`]). A [`ConsensusBackend`] turns a batch of windows into one consensus
//! sequence per window. [`PoaBackend`] is the CPU backend. [`PairAligner`] aligns one sequence
//! to another with the same aligner.
//!
//! # Algorithm
//!
//! [`PoaBackend`] builds a partial-order alignment (POA) graph per window, in the spirit of
//! spoa and abPOA:
//!
//! 1. Bases are normalised: `ACGT` in either case become upper case, every other byte becomes
//!    `N`, which mismatches everything (itself included). Empty segments are dropped.
//! 2. The backbone, the segment at the lower-median position when the segments are sorted by
//!    (length, index), seeds the graph as a chain between a virtual source and sink. The
//!    remaining segments follow in order of increasing length difference from it (ties to the
//!    lower index), so truncated or overlong segments are added last and never shape the
//!    backbone.
//! 3. Each segment is aligned globally to the whole graph (source to sink) with affine gaps
//!    ([`PoaParams`]), then threaded in: aligned bases reuse the node they align to (or a node
//!    with the same base in its alignment column), mismatches and insertions add nodes, and
//!    every edge on the segment's path gains one unit of weight.
//! 4. The consensus is the source-to-sink path with the highest total node score (a
//!    heaviest path in the spirit of pbdagcon). A node scores the number of segments through it minus
//!    the number of segments that skip its alignment column, so a base enters the consensus
//!    when more segments have it than lack the column, and of the alternative bases of a column
//!    the most common wins. Lone insertions, minority deletions and minority substitutions all
//!    lose, including at the window ends, and a truncated segment only adds support to the part
//!    it covers. Ties go to the heavier edge, then to the backbone's edge, then to the lower
//!    node id, so where the votes split evenly (a substitution or an indel alike) the consensus
//!    follows the backbone.
//!
//! The alignment is banded: each graph node is aligned only to the segment positions within
//! [`PoaParams::band_width`] of the best-scoring positions of its predecessors, so the band
//! follows indels instead of a fixed diagonal. The segment is realigned with a doubled band, up
//! to a full (unbanded) alignment, whenever the traced path shows the band may have cut off a
//! better one: the path touches its row's band edge, passes a cell one of whose candidates lay
//! outside a predecessor row's band, or ends while some sink predecessor's band stopped short
//! of the last column. The last two count only when a simple score bound (matches along the
//! longest path to that row, insertions for the rest) says the cut-off path could have scored
//! higher, so the dead ends that truncated segments add mid-graph do not force a full
//! alignment. If a wider band does not fit a fixed cap of 48 MiB of matrices per
//! thread ([`PoaBackend::MATRIX_BYTES_PER_THREAD`]), the best narrower alignment is kept, or
//! the segment is left out (and does not vote). Nothing panics on any input.
//!
//! # Cost
//!
//! Aligning a segment of length `m` to a graph of `V` nodes and `E` edges costs
//! `O((V + E) * (2w + 1))` time for a band half-width `w` that never needed widening, and the
//! same order of memory. For a window of `n` segments of length `L` at HiFi error rates, `V`
//! stays close to `L`, so a window costs about `O(n * L * w)`. HiFi reads are accurate (about
//! 99.9%), so a handful of segments already gives a high-quality consensus; the assembler caps
//! depth at 12 segments per window, which keeps the cost small and bounded.
//!
//! # Limits
//!
//! - Segments carry no base qualities: every segment votes with weight one.
//! - Votes are counted over all segments of the window, so a segment that does not cover a
//!   column counts against it. This suits picoasm's polishing windows, where every segment spans the window
//!   from anchor to anchor; a window made mostly of partial segments loses the parts that fewer
//!   than half of its segments cover.
//! - The consensus is a single path. At a true heterozygous site it follows the majority
//!   allele; phasing must split the reads into haplotypes before consensus.
//! - With two segments, or an even split at a site, the backbone's allele wins; the result is
//!   deterministic but not necessarily the truth. A two-segment window yields the backbone. At
//!   an exact even split of an indel with otherwise equal lengths, the lower median is a
//!   segment of the shorter allele, so such ties lean towards deletions.
//! - The band checks are heuristics: a path that leaves the band and rejoins it far from any
//!   edge is not detected, so on graphs with long indels a narrow band can, rarely, settle for
//!   a slightly suboptimal alignment. The default half-width of 32 makes this negligible for
//!   HiFi windows.
//! - Output is deterministic: it depends only on the window, including the order of its
//!   segments, and the parameters, never on the batch, the thread count or scheduling.
//!   Permuting a window's segments can change the backbone and so the tie-breaks; fill windows
//!   in a canonical order (see [`Window`]).
//!
//! # Examples
//!
//! ```
//! use picoasm_consensus::{ConsensusBackend, PoaBackend, PoaParams, Window};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let backend = PoaBackend::new(PoaParams::default())?;
//! let window = Window::new(vec![
//!     b"ACGTTGCATGCAAGTC".to_vec(),
//!     b"ACGTTGCATCAAGTC".to_vec(),
//!     b"ACGTTGCATGCAAGTC".to_vec(),
//!     b"ACGTAGCATGCAAGTC".to_vec(),
//! ]);
//! let consensus = backend.consensus(&[window, Window::default()]);
//! assert_eq!(consensus, vec![b"ACGTTGCATGCAAGTC".to_vec(), Vec::new()]);
//! # Ok(())
//! # }
//! ```

mod align;
mod backend;
mod graph;
mod pair;
mod params;
mod poa;
mod window;

pub use backend::ConsensusBackend;
pub use pair::{PairAligner, PairOp};
pub use params::{InvalidPoaParamsError, PoaParams};
pub use poa::PoaBackend;
pub use window::Window;
