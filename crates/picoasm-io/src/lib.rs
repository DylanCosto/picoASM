//! Read input and assembly output for the picoasm bacterial assembler.
//!
//! # Input
//!
//! A [`ReadSource`] names the input files of an assembly. Every assembler pass calls
//! [`ReadSource::open`] to get a fresh [`ReadStream`], which yields [`ReadRecord`]s with
//! [`ReadId`](picoasm_seq::ReadId)s `0, 1, 2, ...` in input order across the files, in the order
//! the files were given. Because the files are re-read the same way every time, every pass sees
//! the same id for the same read.
//!
//! The format of each file is detected from its content, never from its extension:
//!
//! 1. The compression layer is sniffed from the first bytes: a gzip member whose extra field
//!    carries the `BC` subfield is BGZF and is decompressed with pasta's
//!    `bgzf::io::MultithreadedReader` (worker count from [`InputOptions::bgzf_workers`]); any
//!    other gzip stream goes through flate2's `MultiGzDecoder` (one thread); anything else is
//!    read as is.
//! 2. The record format is sniffed from the first decompressed bytes: `>` is FASTA, `@` is
//!    FASTQ and the magic number `BAM\1` is BAM. SAM text (a first line starting `@HD`, `@SQ`,
//!    `@RG`, `@PG` or `@CO` and a tab) and CRAM are recognized and rejected with
//!    [`ReadError::UnsupportedFormat`]; anything else is [`ReadError::UnrecognizedFormat`]. An
//!    empty file yields no reads. [`ReadSource::detect_formats`] runs both steps on every input,
//!    so a pass can fail fast before reading any records.
//!
//! FASTA is parsed by pasta's streaming FASTA reader (multi-line sequences are joined). FASTQ is
//! parsed by this crate's [`fastq`] module, which accepts only 4-line records and rejects
//! malformed or truncated input with the offending line number. BAM records are decoded with
//! pasta; secondary and supplementary records are skipped, and records flagged as reverse
//! complemented get their original orientation back, so aligned and unaligned BAM give the
//! same reads. The SAM header text of a BAM file is skipped, not parsed.
//!
//! Read names are the first whitespace-separated word of the FASTA or FASTQ header line (the
//! BAM `QNAME`). Sequences are kept byte for byte (case and non-`ACGT` bytes included; BAM
//! sequences decode to upper case). Reads longer than `u32::MAX` bases are rejected, as are
//! inputs with more than 2³² reads. Empty reads are kept. Text lines are read with a cap
//! (`u32::MAX` bytes plus the line ending), so corrupt input fails with an error after buffering
//! at most that much rather than exhausting memory on one endless line.
//!
//! Parallel callers read [`ReadBatch`]es with [`ReadStream::next_batch`]: each batch holds a
//! contiguous, ordered run of ids, so a caller can process a batch with rayon and still emit
//! results in id order. A batch reuses its buffers between fills, and
//! [`ReadBatch::heap_bytes`] reports the memory it really holds, for charging against a memory
//! budget. [`ReadStream`] is also an [`Iterator`] over owned records.
//!
//! # Output
//!
//! [`FastaWriter`] writes FASTA with a configurable line width. [`GfaWriter`] writes GFA 1
//! header, segment, link and path lines from plain borrowed data ([`GfaSegment`], [`GfaLink`],
//! [`GfaPath`]), so it does not depend on any graph type. Both validate their input and never
//! write a partial line.
//!
//! # Examples
//!
//! ```no_run
//! use picoasm_io::{ReadBatch, ReadSource};
//!
//! # fn main() -> Result<(), picoasm_io::ReadError> {
//! let source = ReadSource::new(["reads.bam", "more_reads.fq.gz"]);
//! let mut stream = source.open()?;
//! let mut batch = ReadBatch::new();
//! while stream.next_batch(&mut batch)? {
//!     for read in batch.records() {
//!         println!("{} {} {}", read.id, String::from_utf8_lossy(&read.name), read.seq.len());
//!     }
//! }
//! # Ok(())
//! # }
//! ```
//!
//! # Limits
//!
//! - Multi-line FASTQ, SAM text and CRAM are not supported. Headerless SAM is not recognized as
//!   SAM and is reported as an unrecognized format.
//! - A single overlong line may still buffer up to 4 GiB before it is rejected.
//! - Plain (non-BGZF) gzip decompresses on one thread.
//! - BGZF detection looks only at the first gzip member's first extra subfield, which is where
//!   every BGZF writer puts `BC`.

mod bam;
mod batch;
mod fasta_writer;
pub mod fastq;
mod format;
mod gfa;
mod read_error;
mod record;
mod source;
mod write_error;

pub use batch::{DEFAULT_BATCH_BASES, MAX_BATCH_RECORDS, ReadBatch};
pub use fasta_writer::{DEFAULT_LINE_WIDTH, FastaWriter};
pub use format::{Compression, InputFormat, RecordFormat, detect_format};
pub use gfa::{
    GfaLink, GfaPath, GfaSegment, GfaTag, GfaWriter, OrientedSegment, Overlap, TagValue,
};
pub use read_error::ReadError;
pub use record::ReadRecord;
pub use source::{InputOptions, ReadSource, ReadStream};
pub use write_error::WriteError;

// Streams are handed to reader threads by the pipeline, and sources are shared between passes.
const _: () = {
    const fn assert_send<T: Send>() {}
    const fn assert_send_sync<T: Send + Sync>() {}
    assert_send::<ReadStream>();
    assert_send_sync::<ReadSource>();
    assert_send_sync::<ReadBatch>();
    assert_send_sync::<ReadError>();
    assert_send_sync::<WriteError>();
};
