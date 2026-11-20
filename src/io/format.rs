//! Content-based detection of input compression and record format.

use std::{
    fmt,
    fs::File,
    io::{BufRead, BufReader},
    num::NonZero,
    path::Path,
};

use flate2::bufread::MultiGzDecoder;
use pasta::bgzf;

use crate::io::ReadError;

/// Capacity of the buffers between the file, the decompressor and the parsers.
const BUFFER_CAPACITY: usize = 1 << 20;

/// The gzip magic bytes.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// The BAM magic number (SAM specification §4.2).
pub(crate) const BAM_MAGIC: [u8; 4] = *b"BAM\x01";

/// The CRAM magic bytes (CRAM specification §6).
const CRAM_MAGIC: [u8; 4] = *b"CRAM";

/// The SAM header record types that can start a SAM file, each followed by a tab.
///
/// htslib recognizes SAM text by the same prefixes. A FASTQ header such as `@HD\t...` would need a
/// read named `HD` followed by a tab, which no sequencer writes.
const SAM_HEADER_PREFIXES: [&[u8]; 5] = [b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"];

/// How an input file is compressed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Compression {
    /// Not compressed.
    None,
    /// gzip that is not BGZF, decompressed on one thread.
    Gzip,
    /// BGZF (blocked gzip, as used by BAM and `bgzip`), decompressed on several threads.
    Bgzf,
}

/// What kind of records an input file holds once decompressed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum RecordFormat {
    /// FASTA, single- or multi-line.
    Fasta,
    /// FASTQ with 4-line records.
    Fastq,
    /// BAM, aligned or unaligned.
    Bam,
    /// No data at all; the file yields no reads.
    Empty,
}

/// The detected format of one input file.
///
/// Marked `#[non_exhaustive]` so that later details (such as a BAM sort order) can be added;
/// build one with [`InputFormat::new`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct InputFormat {
    /// The compression layer.
    pub compression: Compression,
    /// The record format inside it.
    pub records: RecordFormat,
}

impl InputFormat {
    /// Creates a format from its compression layer and record format.
    #[must_use]
    pub const fn new(compression: Compression, records: RecordFormat) -> Self {
        Self {
            compression,
            records,
        }
    }
}

impl fmt::Display for Compression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::None => "uncompressed",
            Self::Gzip => "gzip",
            Self::Bgzf => "BGZF",
        })
    }
}

impl fmt::Display for RecordFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Fasta => "FASTA",
            Self::Fastq => "FASTQ",
            Self::Bam => "BAM",
            Self::Empty => "empty",
        })
    }
}

impl fmt::Display for InputFormat {
    /// Formats as, for example, `FASTQ (gzip)` or `BAM (BGZF)`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.records, self.compression)
    }
}

/// A decompressed input stream.
pub(crate) type Decoded = Box<dyn BufRead + Send>;

/// An opened input file: its format and its decompressed content, positioned at the start.
pub(crate) struct Opened {
    pub(crate) format: InputFormat,
    pub(crate) reader: Decoded,
}

/// Detects the format of the file at `path` from its first bytes.
///
/// Opens the file, sniffs the compression and decompresses just enough to sniff the record
/// format. Useful to check every input before a long pass starts; [`ReadSource::detect_formats`]
/// does this for all inputs.
///
/// [`ReadSource::detect_formats`]: crate::io::ReadSource::detect_formats
///
/// # Errors
///
/// Returns [`ReadError::Open`] when the file cannot be opened, [`ReadError::Io`] when its first
/// bytes cannot be read or decompressed, [`ReadError::UnsupportedFormat`] for SAM text or CRAM,
/// and [`ReadError::UnrecognizedFormat`] for any other content that is not FASTA, FASTQ or BAM
/// (a BAM file must start with the full BAM magic number).
pub fn detect_format(path: impl AsRef<Path>) -> Result<InputFormat, ReadError> {
    open(path.as_ref(), NonZero::<usize>::MIN).map(|opened| opened.format)
}

/// Opens `path` and returns its format and decompressed content.
pub(crate) fn open(path: &Path, bgzf_workers: NonZero<usize>) -> Result<Opened, ReadError> {
    let file = File::open(path).map_err(|source| ReadError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    let mut raw = BufReader::with_capacity(BUFFER_CAPACITY, file);
    let head = raw.fill_buf().map_err(|source| ReadError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let compression = sniff_compression(head);
    let mut reader: Decoded = match compression {
        Compression::None => Box::new(raw),
        Compression::Gzip => Box::new(BufReader::with_capacity(
            BUFFER_CAPACITY,
            MultiGzDecoder::new(raw),
        )),
        Compression::Bgzf => Box::new(
            bgzf::io::MultithreadedReader::try_with_worker_count(bgzf_workers, raw).map_err(
                |error| ReadError::Threads {
                    path: path.to_path_buf(),
                    source: error.into_error(),
                },
            )?,
        ),
    };
    let head = reader.fill_buf().map_err(|source| ReadError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let records = sniff_records(head).map_err(|error| match error {
        SniffError::Unrecognized { first_byte } => ReadError::UnrecognizedFormat {
            path: path.to_path_buf(),
            first_byte,
        },
        SniffError::Unsupported { format } => ReadError::UnsupportedFormat {
            path: path.to_path_buf(),
            format,
        },
    })?;
    Ok(Opened {
        format: InputFormat::new(compression, records),
        reader,
    })
}

/// Classifies the first raw bytes of a file.
fn sniff_compression(head: &[u8]) -> Compression {
    if !head.starts_with(&GZIP_MAGIC) {
        Compression::None
    } else if is_bgzf_header(head) {
        Compression::Bgzf
    } else {
        Compression::Gzip
    }
}

/// Returns `true` when `head` starts with a gzip member header carrying the BGZF `BC` subfield.
///
/// Layout (SAM specification §4.1): ID1 ID2 CM FLG MTIME(4) XFL OS XLEN(2), then the first
/// extra subfield SI1 = `B`, SI2 = `C`, SLEN = 2.
fn is_bgzf_header(head: &[u8]) -> bool {
    const FEXTRA: u8 = 0x04;
    const DEFLATE: u8 = 8;
    head.len() >= 16
        && head[2] == DEFLATE
        && head[3] & FEXTRA != 0
        && head[12..16] == [b'B', b'C', 2, 0]
}

/// Why the first decompressed bytes of a file are not a supported record format.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SniffError {
    /// No known format starts this way.
    Unrecognized { first_byte: u8 },
    /// A known format that this crate does not read.
    Unsupported { format: &'static str },
}

/// Classifies the first decompressed bytes of a file.
///
/// BAM needs its full 4-byte magic number. When the decompressor hands over fewer than 4 bytes
/// (possible only for a tiny first gzip member), a prefix of the magic number is accepted and the
/// BAM header reader checks the rest.
fn sniff_records(head: &[u8]) -> Result<RecordFormat, SniffError> {
    let Some(&first_byte) = head.first() else {
        return Ok(RecordFormat::Empty);
    };
    let has_prefix = |magic: &[u8]| {
        if head.len() >= magic.len() {
            head.starts_with(magic)
        } else {
            magic.starts_with(head)
        }
    };
    if first_byte == b'@' && SAM_HEADER_PREFIXES.iter().any(|p| head.starts_with(p)) {
        return Err(SniffError::Unsupported { format: "SAM" });
    }
    match first_byte {
        b'>' => Ok(RecordFormat::Fasta),
        b'@' => Ok(RecordFormat::Fastq),
        b'B' if has_prefix(&BAM_MAGIC) => Ok(RecordFormat::Bam),
        b'C' if has_prefix(&CRAM_MAGIC) => Err(SniffError::Unsupported { format: "CRAM" }),
        _ => Err(SniffError::Unrecognized { first_byte }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sniffs_compression() {
        assert_eq!(sniff_compression(b">r\nACGT\n"), Compression::None);
        assert_eq!(
            sniff_compression(&[0x1f, 0x8b, 8, 0, 0, 0]),
            Compression::Gzip
        );
        let bgzf = [
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0,
        ];
        assert_eq!(sniff_compression(&bgzf), Compression::Bgzf);
        assert_eq!(sniff_compression(&bgzf[..15]), Compression::Gzip);
    }

    #[test]
    fn sniffs_records() {
        let unrecognized = |first_byte| Err(SniffError::Unrecognized { first_byte });
        assert_eq!(sniff_records(b""), Ok(RecordFormat::Empty));
        assert_eq!(sniff_records(b">r"), Ok(RecordFormat::Fasta));
        assert_eq!(sniff_records(b"@r"), Ok(RecordFormat::Fastq));
        assert_eq!(sniff_records(b"@HDX\tdesc"), Ok(RecordFormat::Fastq));
        assert_eq!(sniff_records(b"BAM\x01\x00"), Ok(RecordFormat::Bam));
        assert_eq!(sniff_records(b"BA"), Ok(RecordFormat::Bam));
        assert_eq!(sniff_records(b"BAM\x02"), unrecognized(b'B'));
        assert_eq!(sniff_records(b"Bx"), unrecognized(b'B'));
        assert_eq!(sniff_records(b"ACGT"), unrecognized(b'A'));
        for sam in [&b"@HD\tVN:1.6"[..], b"@SQ\tSN:chr1", b"@CO\tnote"] {
            assert_eq!(
                sniff_records(sam),
                Err(SniffError::Unsupported { format: "SAM" }),
                "{sam:?}"
            );
        }
        assert_eq!(
            sniff_records(b"CRAM\x03\x00"),
            Err(SniffError::Unsupported { format: "CRAM" })
        );
    }

    #[test]
    fn displays_formats() {
        let format = InputFormat::new(Compression::Gzip, RecordFormat::Fastq);
        assert_eq!(format.to_string(), "FASTQ (gzip)");
    }
}
