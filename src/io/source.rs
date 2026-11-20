//! Input files and the read streams that pass over them.

use std::{
    io::{self, BufRead, Read},
    iter::FusedIterator,
    num::NonZero,
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::seq::ReadId;
use pasta::fasta;

use crate::io::{
    ReadBatch, ReadError, ReadRecord,
    bam::BamReads,
    fastq,
    format::{self, Decoded, InputFormat, RecordFormat},
};

/// Options for reading inputs.
///
/// Built with [`InputOptions::new`] (or [`Default`]) and the `with_*` methods.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct InputOptions {
    bgzf_workers: NonZero<usize>,
}

impl InputOptions {
    /// Returns the default options: one BGZF decompression worker.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            bgzf_workers: NonZero::<usize>::MIN,
        }
    }

    /// Sets the number of BGZF decompression worker threads per open file.
    ///
    /// A BGZF file also uses one extra thread to read compressed blocks. Plain gzip always
    /// decompresses on the calling thread.
    #[must_use]
    pub const fn with_bgzf_workers(mut self, workers: NonZero<usize>) -> Self {
        self.bgzf_workers = workers;
        self
    }

    /// Returns the number of BGZF decompression worker threads per open file.
    #[must_use]
    pub const fn bgzf_workers(&self) -> NonZero<usize> {
        self.bgzf_workers
    }
}

impl Default for InputOptions {
    fn default() -> Self {
        Self::new()
    }
}

/// The input files of an assembly, re-readable once per pass.
///
/// Cloning is cheap: the path list is shared.
#[derive(Clone, Debug)]
pub struct ReadSource {
    paths: Arc<[PathBuf]>,
    options: InputOptions,
}

impl ReadSource {
    /// Creates a source over `paths`, read in the given order, with default options.
    ///
    /// Nothing is opened until [`Self::open`] or [`Self::detect_formats`].
    pub fn new<I, P>(paths: I) -> Self
    where
        I: IntoIterator<Item = P>,
        P: Into<PathBuf>,
    {
        Self {
            paths: paths.into_iter().map(Into::into).collect(),
            options: InputOptions::new(),
        }
    }

    /// Replaces the options.
    #[must_use]
    pub fn with_options(mut self, options: InputOptions) -> Self {
        self.options = options;
        self
    }

    /// Returns the input paths, in reading order.
    #[must_use]
    pub fn paths(&self) -> &[PathBuf] {
        &self.paths
    }

    /// Returns the options.
    #[must_use]
    pub fn options(&self) -> InputOptions {
        self.options
    }

    /// Detects the format of every input, in order, without reading any records.
    ///
    /// Call this before a long pass to fail fast on a missing or unsupported file.
    ///
    /// # Errors
    ///
    /// Returns the first error from [`detect_format`](crate::io::detect_format).
    pub fn detect_formats(&self) -> Result<Vec<InputFormat>, ReadError> {
        self.paths.iter().map(format::detect_format).collect()
    }

    /// Starts a pass over all inputs, with read ids starting at 0.
    ///
    /// The first file is opened (and its format detected) immediately; later files are opened
    /// when the stream reaches them.
    ///
    /// # Errors
    ///
    /// Returns an error when the first file cannot be opened or its format is not recognized.
    pub fn open(&self) -> Result<ReadStream, ReadError> {
        let mut stream = ReadStream {
            paths: Arc::clone(&self.paths),
            options: self.options,
            next_file: 0,
            current: None,
            next_id: 0,
            finished: false,
        };
        if !stream.paths.is_empty() {
            stream.current = Some(OpenFile::open(&stream.paths[0], stream.options)?);
            stream.next_file = 1;
        }
        Ok(stream)
    }
}

/// One pass over the inputs of a [`ReadSource`], yielding reads in id order.
///
/// Read records one at a time with [`Self::read_record`] (reusing a buffer) or through the
/// [`Iterator`] impl (allocating), or in batches with [`Self::next_batch`]. After an error the
/// stream is finished: every later call returns no reads.
pub struct ReadStream {
    paths: Arc<[PathBuf]>,
    options: InputOptions,
    /// Index into `paths` of the next file to open.
    next_file: usize,
    current: Option<OpenFile>,
    /// The id of the next read; may reach 2³², which is then reported as too many reads.
    next_id: u64,
    finished: bool,
}

impl ReadStream {
    /// Reads the next read into `record`, reusing its buffers.
    ///
    /// Returns `Ok(false)` once every input is exhausted; `record` is then unspecified.
    ///
    /// # Errors
    ///
    /// Returns a [`ReadError`] when a file cannot be opened or parsed, when a read is longer than
    /// `u32::MAX` bases, or when the input holds more than 2³² reads. The stream is finished
    /// afterwards.
    pub fn read_record(&mut self, record: &mut ReadRecord) -> Result<bool, ReadError> {
        if self.finished {
            return Ok(false);
        }
        let result = self.read_record_inner(record);
        if !matches!(result, Ok(true)) {
            self.finished = true;
            self.current = None;
        }
        result
    }

    /// Replaces the contents of `batch` with the next reads, returning `Ok(false)` when there
    /// are none left.
    ///
    /// Reads are added until the batch holds at least [`ReadBatch::max_bases`] bases (or
    /// [`MAX_BATCH_RECORDS`](crate::io::MAX_BATCH_RECORDS) reads), or the input ends; a batch filled
    /// with `Ok(true)` holds at least one read. Ids within a batch are contiguous and follow the
    /// previous batch's.
    ///
    /// Spare record buffers left over from earlier, larger batches are freed, so the batch's
    /// [`heap_bytes`](ReadBatch::heap_bytes) follows its contents.
    ///
    /// # Errors
    ///
    /// As for [`Self::read_record`]. The batch then holds the reads read before the error.
    pub fn next_batch(&mut self, batch: &mut ReadBatch) -> Result<bool, ReadError> {
        batch.clear();
        let filled = self.fill_batch(batch);
        batch.release_spare_slots();
        filled.map(|()| !batch.is_empty())
    }

    /// Adds reads to the empty `batch` until it is full or the input ends.
    fn fill_batch(&mut self, batch: &mut ReadBatch) -> Result<(), ReadError> {
        while !batch.is_full() {
            let slot = batch.next_slot();
            if !self.read_record(slot)? {
                break;
            }
            batch.commit_slot();
        }
        Ok(())
    }

    /// Returns the number of reads yielded so far, which is also the id of the next read.
    #[must_use]
    pub fn reads_read(&self) -> u64 {
        self.next_id
    }

    /// Returns the path of the file being read.
    ///
    /// [`ReadSource::open`] opens the first file at once, so this is `None` only for a source
    /// without files, once every file is exhausted, and after an error. Between the last read of
    /// one file and the first read of the next it still names the earlier file.
    #[must_use]
    pub fn current_path(&self) -> Option<&Path> {
        self.current.as_ref().map(|file| file.path.as_path())
    }

    /// Returns the detected format of the file being read, under the same rules as
    /// [`Self::current_path`].
    #[must_use]
    pub fn current_format(&self) -> Option<InputFormat> {
        self.current.as_ref().map(|file| file.format)
    }

    fn read_record_inner(&mut self, record: &mut ReadRecord) -> Result<bool, ReadError> {
        loop {
            let file = match &mut self.current {
                Some(file) => file,
                slot @ None => {
                    let Some(path) = self.paths.get(self.next_file) else {
                        return Ok(false);
                    };
                    let file = OpenFile::open(path, self.options)?;
                    self.next_file += 1;
                    slot.insert(file)
                }
            };
            if !file.read(&mut record.name, &mut record.seq)? {
                log::debug!(
                    "finished {}: {} records",
                    file.path.display(),
                    file.records_read()
                );
                self.current = None;
                continue;
            }
            if record.seq.len() > MAX_READ_LEN {
                return Err(ReadError::ReadTooLong {
                    path: file.path.clone(),
                    record: file.records_read(),
                    name: String::from_utf8_lossy(&record.name).into_owned(),
                    len: record.seq.len(),
                });
            }
            let id = u32::try_from(self.next_id).map_err(|_| ReadError::TooManyReads {
                path: file.path.clone(),
            })?;
            record.id = ReadId(id);
            self.next_id += 1;
            return Ok(true);
        }
    }
}

impl std::fmt::Debug for ReadStream {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReadStream")
            .field("current_path", &self.current_path())
            .field("current_format", &self.current_format())
            .field("reads_read", &self.next_id)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl Iterator for ReadStream {
    type Item = Result<ReadRecord, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = ReadRecord::default();
        match self.read_record(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

impl FusedIterator for ReadStream {}

/// The longest read, in bases: positions are `u32`.
const MAX_READ_LEN: usize = u32::MAX as usize;

/// The longest FASTA header line, in bytes: as long as the longest read.
const MAX_FASTA_HEADER_LEN: usize = MAX_READ_LEN;

/// An open input file with its record parser.
struct OpenFile {
    path: PathBuf,
    format: InputFormat,
    parser: Parser,
}

/// The record parser of one file.
enum Parser {
    Fasta {
        reader: fasta::io::Reader<Decoded>,
        header: Vec<u8>,
        records_read: u64,
    },
    Fastq {
        reader: fastq::Reader<Decoded>,
        quality: Vec<u8>,
        records_read: u64,
    },
    Bam(BamReads<Decoded>),
    Empty,
}

impl OpenFile {
    fn open(path: &Path, options: InputOptions) -> Result<Self, ReadError> {
        let opened = format::open(path, options.bgzf_workers())?;
        log::debug!("reading {} as {}", path.display(), opened.format);
        let parser = match opened.format.records {
            RecordFormat::Fasta => Parser::Fasta {
                reader: fasta::io::Reader::new(opened.reader),
                header: Vec::new(),
                records_read: 0,
            },
            RecordFormat::Fastq => Parser::Fastq {
                reader: fastq::Reader::new(opened.reader),
                quality: Vec::new(),
                records_read: 0,
            },
            RecordFormat::Bam => Parser::Bam(BamReads::new(opened.reader).map_err(|source| {
                ReadError::BamHeader {
                    path: path.to_path_buf(),
                    source,
                }
            })?),
            RecordFormat::Empty => {
                log::warn!("input {} is empty", path.display());
                Parser::Empty
            }
        };
        Ok(Self {
            path: path.to_path_buf(),
            format: opened.format,
            parser,
        })
    }

    /// Returns the number of records read from this file so far.
    fn records_read(&self) -> u64 {
        match &self.parser {
            Parser::Fasta { records_read, .. } | Parser::Fastq { records_read, .. } => {
                *records_read
            }
            Parser::Bam(reads) => reads.records_read(),
            Parser::Empty => 0,
        }
    }

    /// Reads the next read's name and bases, returning `Ok(false)` at the end of the file.
    fn read(&mut self, name: &mut Vec<u8>, seq: &mut Vec<u8>) -> Result<bool, ReadError> {
        match &mut self.parser {
            Parser::Fasta {
                reader,
                header,
                records_read,
            } => {
                let record = *records_read + 1;
                let fasta_error = |source| ReadError::Fasta {
                    path: self.path.clone(),
                    record,
                    source,
                };
                if !read_fasta_header(reader.get_mut(), header, name).map_err(fasta_error)? {
                    return Ok(false);
                }
                seq.clear();
                // One byte past the limit is enough to see that a read is too long.
                let limit = u64::try_from(MAX_READ_LEN).unwrap_or(u64::MAX) + 1;
                reader
                    .sequence_reader()
                    .take(limit)
                    .read_to_end(seq)
                    .map_err(fasta_error)?;
                *records_read = record;
                Ok(true)
            }
            Parser::Fastq {
                reader,
                quality,
                records_read,
            } => {
                let record = *records_read + 1;
                let more =
                    reader
                        .read_record(name, seq, quality)
                        .map_err(|source| ReadError::Fastq {
                            path: self.path.clone(),
                            record,
                            source,
                        })?;
                if more {
                    *records_read = record;
                }
                Ok(more)
            }
            Parser::Bam(reads) => reads.read(name, seq).map_err(|source| ReadError::Bam {
                path: self.path.clone(),
                record: reads.records_read() + 1,
                source,
            }),
            Parser::Empty => Ok(false),
        }
    }
}

/// Reads a FASTA header line into `line` and its name (the first word after `>`) into `name`.
///
/// Returns `Ok(false)` at the end of the input. pasta's sequence reader leaves the input at a
/// `>` or at the end, so the line always starts with `>` unless the file is malformed. At most
/// [`MAX_FASTA_HEADER_LEN`] + 2 bytes are buffered.
fn read_fasta_header<R: BufRead>(
    reader: &mut R,
    line: &mut Vec<u8>,
    name: &mut Vec<u8>,
) -> io::Result<bool> {
    let invalid = |message: &str| io::Error::new(io::ErrorKind::InvalidData, message.to_owned());
    line.clear();
    let limit = u64::try_from(MAX_FASTA_HEADER_LEN)
        .unwrap_or(u64::MAX)
        .saturating_add(2);
    if reader.take(limit).read_until(b'\n', line)? == 0 {
        return Ok(false);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
        if line.last() == Some(&b'\r') {
            line.pop();
        }
    }
    if line.len() > MAX_FASTA_HEADER_LEN {
        return Err(invalid(&format!(
            "header line is longer than {MAX_FASTA_HEADER_LEN} bytes"
        )));
    }
    let Some(rest) = line.strip_prefix(b">") else {
        return Err(invalid("expected a header line starting with '>'"));
    };
    let name_end = rest
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(rest.len());
    if name_end == 0 {
        return Err(invalid("header line has no name"));
    }
    name.clear();
    name.extend_from_slice(&rest[..name_end]);
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(data: &[u8]) -> io::Result<Option<Vec<u8>>> {
        let (mut line, mut name) = (Vec::new(), Vec::new());
        let more = read_fasta_header(&mut &data[..], &mut line, &mut name)?;
        Ok(more.then_some(name))
    }

    #[test]
    fn parses_fasta_headers() {
        assert_eq!(header(b">r1 desc\nACGT").ok(), Some(Some(b"r1".to_vec())));
        assert_eq!(header(b">r1\r\n").ok(), Some(Some(b"r1".to_vec())));
        assert_eq!(header(b">r1\tx").ok(), Some(Some(b"r1".to_vec())));
        assert_eq!(header(b"").ok(), Some(None));
        for bad in [&b">\nACGT"[..], b"> r1\n", b"r1\n"] {
            assert!(header(bad).is_err(), "{bad:?}");
        }
    }
}
