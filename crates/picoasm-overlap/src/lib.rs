//! Base-level read overlaps, read correction and the string-graph layout, as hifiasm builds them.
//!
//! Reads are overlapped base by base, each read is corrected only from reads of its own copy
//! of a repeat, and the corrected reads are laid out in a string graph whose unitigs are the
//! draft contigs.
//!
//! # Examples
//!
//! ```
//! use picoasm_overlap::{OverlapParams, Reads, find_overlaps};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let mut state = 1_u64;
//! let genome: Vec<u8> = (0..4_000)
//!     .map(|_| {
//!         state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
//!         b"ACGT"[(state >> 62) as usize]
//!     })
//!     .collect();
//! let reads: Reads = [&genome[..2_500], &genome[1_000..]].into_iter().collect();
//! let params = OverlapParams { density: 0.2, ..OverlapParams::default() };
//! let overlaps = find_overlaps(&reads, &params)?;
//! assert_eq!(overlaps.len(), 1);
//! assert_eq!((overlaps[0].query_start, overlaps[0].query_end), (1_000, 2_500));
//! # Ok(())
//! # }
//! ```

mod align;
mod correct;
mod layout;
mod overlap;
mod reads;
mod rounds;

pub use align::{DiffAligner, Edit, EditOp, Extent, Stall};
pub use correct::{Correction, CorrectionParams, correct};
pub use layout::{
    Contig, Ends, Fate, Layout, LayoutParams, MAX_DIVERGENCE, MAX_RESCUE_DIVERGENCE, Placed, Stop,
    Usage, lay_out, usage,
};
pub use overlap::{
    CHAINS_PER_PAIR, DIAGONAL_SLACK, Extension, InvalidOverlapParamsError, Overlap, OverlapParams,
    Overlaps, Pair, find_overlaps, find_overlaps_counted, find_overlaps_skipping, map_reads,
};
pub use reads::{ReadIndex, Reads};
pub use rounds::{Corrected, RoundStats, correct_in_rounds};
