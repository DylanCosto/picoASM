//! Sequence primitives for the picoasm bacterial assembler.
//!
//! This crate turns raw reads into the sparse, strand-independent features the rest of the
//! assembler works with: homopolymer-compressed sequence, [`Anchor`]s (universe minimizers chosen
//! in homopolymer-compressed space) and [`KMinMer`]s (windows of consecutive anchors).
//!
//! # Algorithm
//!
//! 1. [`Hpc`] collapses each homopolymer run to one base and remembers where every run starts in
//!    the raw read, because HiFi errors are mostly run-length errors.
//! 2. [`Sketcher`] rolls the 2-bit encodings of every compressed k-mer and of its reverse
//!    complement, hashes the smaller one with the invertible [`mix64`], and keeps the k-mer as an
//!    [`Anchor`] when the hash falls below `density * 2^64` (a universe minimizer: the choice
//!    depends only on the k-mer, never on its neighbours). Anchors carry raw coordinates.
//! 3. [`kminmers`] slides a window of [`SketchParams::window`] anchors over a read and hashes
//!    each window's anchor hashes in both directions; the smaller hash is canonical.
//!
//! # Limits
//!
//! - Anchors are at most [`MAX_KMER_LEN`] (31) compressed bases, and read positions are `u32`,
//!   so reads longer than `u32::MAX` bases are rejected with [`ReadTooLongError`].
//! - Bases other than `A`, `C`, `G` and `T` (in either case) break k-mers: no anchor spans them.
//! - Reverse-complement palindromic k-mers and palindromic anchor-hash windows have no
//!   orientation and are skipped, so a k-min-mer's index in a read's list can differ from its
//!   [`KMinMer::first_anchor`].
//! - The k-min-mer hash ignores anchor strands, as in mdBG, so windows that differ only in
//!   anchor strands (inverted repeats at anchor scale) share a hash.
//!
//! # Examples
//!
//! ```
//! use picoasm::seq::{SketchParams, Sketcher, kminmers};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let params = SketchParams::new(5, 0.5, 3)?;
//! let mut sketcher = Sketcher::new(params);
//! let mut anchors = Vec::new();
//! sketcher.sketch(b"ACGTTGCATGCAAGTCCGATCGATTACG", &mut anchors)?;
//!
//! let mut windows = Vec::new();
//! kminmers(&anchors, params.window(), &mut windows);
//! // Each k-min-mer covers `window` anchors from its first one. Palindromic windows are
//! // skipped, so index in `windows` and `first_anchor` can differ.
//! for w in &windows {
//!     let first = usize::try_from(w.first_anchor)?;
//!     assert!(first + params.window() <= anchors.len());
//! }
//! # Ok(())
//! # }
//! ```

mod base;
mod hash;
mod hpc;
mod kminmer;
mod read_id;
mod read_len;
mod sketch;
mod strand;

pub use base::{complement, encode, reverse_complement, reverse_complement_into};
pub use hash::{combine, mix64};
pub use hpc::Hpc;
pub use kminmer::{KMinMer, kminmers};
pub use read_id::ReadId;
pub use read_len::{ReadTooLongError, checked_read_len};
pub use sketch::{Anchor, InvalidSketchParamsError, MAX_KMER_LEN, SketchParams, Sketcher};
pub use strand::Strand;
