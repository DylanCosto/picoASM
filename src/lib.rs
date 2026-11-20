//! Bacterial genome assembly from Oxford Nanopore or PacBio HiFi reads.
//!
//! A bacterial isolate is a few megabases in a handful of replicons, most of them circular: a
//! chromosome and plasmids from about a kilobase up. Its reads are overlapped base by base
//! and laid out in a string graph ([`crate::overlap`]), whose unitigs are the draft contigs; a
//! path that closes on itself is a whole circular replicon. The drafts are then polished from
//! all the reads by partial-order consensus in windows ([`crate::consensus`]). See `README.md`
//! for the benchmark.
//!
//! # Examples
//!
//! ```
//! use picoasm::{BactParams, ReadKind, assemble};
//! use picoasm::overlap::Reads;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let reads = Reads::new();
//! let assembly = assemble(&reads, &BactParams::new(ReadKind::Ont))?;
//! assert!(assembly.contigs.is_empty());
//! # Ok(())
//! # }
//! ```
//!
//! The building blocks are public modules: [`seq`] (bases, homopolymer compression, anchors and
//! k-min-mers), [`io`] (read input and assembly output), [`overlap`] (base-level overlap assembly)
//! and [`consensus`] (per-window base-level consensus).

pub mod consensus;
pub mod io;
pub mod overlap;
pub mod seq;

mod assemble;
mod bridge;
mod check;
mod circle;
mod context;
mod correct;
mod hairpin;
mod hang;
mod homopolymer;
mod join;
mod locus;
mod motif;
mod params;
mod polish;
mod repeat;
mod select;
mod short;
mod small;
mod trim;
mod vote;

pub use assemble::{
    AssembleError, Assembly, Contig, MAX_EXCISED_DEPTH_SHARE, MIN_DEPTH_SHARE,
    MIN_LINEAR_DEPTH_SHARE, MIN_WINDOW_LEN, assemble,
};
pub use check::{BIN as CHECK_BIN, Check, Flag, FlaggedBin};
pub use homopolymer::{Evidence, Run};
pub use params::{BactParams, ReadKind};
pub use short::{ShortParams, ShortStats, polish_short};
