//! Errors from reading input files.

use std::{
    io,
    path::{Path, PathBuf},
};

use crate::io::fastq;

/// An error from opening, detecting or parsing an input file.
///
/// Every variant names the file. Record numbers are 1-based and count the records of that file
/// (for BAM, skipped secondary and supplementary records included).
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReadError {
    /// The file could not be opened.
    #[error("cannot open {}", path.display())]
    Open {
        /// The input file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// Reading or decompressing the file failed outside a record (for example while sniffing).
    #[error("cannot read {}", path.display())]
    Io {
        /// The input file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// Not even one BGZF decompression thread could be started.
    #[error("cannot start BGZF decompression threads for {}", path.display())]
    Threads {
        /// The input file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// The decompressed content is not FASTA, FASTQ or BAM.
    #[error(
        "unrecognized input format in {} (first byte {first_byte:#04x}); expected FASTA, FASTQ \
         or BAM, optionally gzip- or BGZF-compressed",
        path.display()
    )]
    UnrecognizedFormat {
        /// The input file.
        path: PathBuf,
        /// The first decompressed byte.
        first_byte: u8,
    },
    /// The file holds a known format that is not supported, such as SAM text or CRAM.
    #[error(
        "{format} input is not supported ({}); convert it to BAM, FASTQ or FASTA",
        path.display()
    )]
    UnsupportedFormat {
        /// The input file.
        path: PathBuf,
        /// The name of the detected format, such as `SAM` or `CRAM`.
        format: &'static str,
    },
    /// A FASTA record could not be read.
    #[error("cannot read FASTA record {record} of {}", path.display())]
    Fasta {
        /// The input file.
        path: PathBuf,
        /// The 1-based record number in the file.
        record: u64,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A FASTQ record is malformed or could not be read; the source has the line number.
    #[error("cannot read FASTQ record {record} of {}", path.display())]
    Fastq {
        /// The input file.
        path: PathBuf,
        /// The 1-based record number in the file.
        record: u64,
        /// The cause, with its line number.
        #[source]
        source: fastq::ParseError,
    },
    /// The BAM header is invalid or could not be read.
    #[error("cannot read the BAM header of {}", path.display())]
    BamHeader {
        /// The input file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A BAM record is invalid or could not be read.
    #[error("cannot read BAM record {record} of {}", path.display())]
    Bam {
        /// The input file.
        path: PathBuf,
        /// The 1-based record number in the file.
        record: u64,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A read is longer than `u32::MAX` bases, the longest read positions can address.
    #[error(
        "read {name:?} (record {record} of {}) has at least {len} bases, more than the limit of \
         {}",
        path.display(),
        u32::MAX
    )]
    ReadTooLong {
        /// The input file.
        path: PathBuf,
        /// The 1-based record number in the file.
        record: u64,
        /// The read name, lossily converted to UTF-8.
        name: String,
        /// The read length, or a lower bound on it: reading a sequence stops one base past the
        /// limit, so it is never buffered whole.
        len: usize,
    },
    /// The input holds more than 2³² reads, more than [`ReadId`](crate::seq::ReadId) can number.
    #[error("more than 2^32 reads in the input (reached in {})", path.display())]
    TooManyReads {
        /// The file in which the limit was reached.
        path: PathBuf,
    },
}

impl ReadError {
    /// Returns the input file the error refers to.
    #[must_use]
    pub fn path(&self) -> &Path {
        match self {
            Self::Open { path, .. }
            | Self::Io { path, .. }
            | Self::Threads { path, .. }
            | Self::UnrecognizedFormat { path, .. }
            | Self::UnsupportedFormat { path, .. }
            | Self::Fasta { path, .. }
            | Self::Fastq { path, .. }
            | Self::BamHeader { path, .. }
            | Self::Bam { path, .. }
            | Self::ReadTooLong { path, .. }
            | Self::TooManyReads { path } => path,
        }
    }
}
