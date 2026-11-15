//! A streaming parser for 4-line FASTQ records.
//!
//! Each record is exactly four lines: `@name [description]`, the sequence, a separator line
//! starting with `+` (anything after the `+` is ignored), and a quality line with one byte per
//! base. Line endings may be `\n` or `\r\n`, and the last line may lack its newline. Empty lines
//! are skipped only where a record header is expected, so trailing blank lines are accepted.
//!
//! Malformed input is an error carrying the 1-based line number: a header without `@` or
//! without a name, a missing separator (which is also how multi-line FASTQ shows up), a
//! quality line whose length differs from the sequence, or input that ends inside a record.
//! Quality bytes are not otherwise validated.
//!
//! No line is buffered past [`Reader::max_line_len`] bytes (by default [`DEFAULT_MAX_LINE_LEN`],
//! the longest read that `u32` coordinates can address), so corrupt or misidentified input fails
//! with [`ParseError::LineTooLong`] instead of being read whole into memory.
//!
//! The module depends only on `std` and `thiserror`, so it can move to another crate as is.
//!
//! # Examples
//!
//! ```
//! use picoasm_io::fastq::Reader;
//!
//! # fn main() -> Result<(), picoasm_io::fastq::ParseError> {
//! let data = b"@r1 first read\nACGT\n+\nIIII\n@r2\nGG\n+r2\n!!\n";
//! let mut reader = Reader::new(&data[..]);
//! let (mut name, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
//!
//! assert!(reader.read_record(&mut name, &mut seq, &mut qual)?);
//! assert_eq!((name.as_slice(), seq.as_slice()), (&b"r1"[..], &b"ACGT"[..]));
//! assert!(reader.read_record(&mut name, &mut seq, &mut qual)?);
//! assert_eq!(qual, b"!!");
//! assert!(!reader.read_record(&mut name, &mut seq, &mut qual)?);
//! # Ok(())
//! # }
//! ```

use std::{
    fmt,
    io::{self, BufRead, Read},
};

/// The default longest line, in bytes without the line ending: `u32::MAX`.
pub const DEFAULT_MAX_LINE_LEN: usize = u32::MAX as usize;

/// A line of a FASTQ record after the header.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum RecordLine {
    /// The sequence line.
    Sequence,
    /// The `+` separator line.
    Separator,
    /// The quality line.
    Quality,
}

impl fmt::Display for RecordLine {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Sequence => "sequence",
            Self::Separator => "separator",
            Self::Quality => "quality",
        })
    }
}

/// An error from parsing FASTQ, with the 1-based line number where it was found.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ParseError {
    /// The underlying reader failed.
    #[error("i/o error at line {line}")]
    Io {
        /// The line being read.
        line: u64,
        /// The reader's error.
        #[source]
        source: io::Error,
    },
    /// A record header does not start with `@`.
    #[error("line {line}: expected a record header starting with '@'")]
    MissingHeaderPrefix {
        /// The header line.
        line: u64,
    },
    /// A record header has nothing (or only whitespace) after the `@`.
    #[error("line {line}: record header has no name")]
    MissingName {
        /// The header line.
        line: u64,
    },
    /// The line after the sequence does not start with `+`.
    #[error("line {line}: expected a '+' separator line (multi-line FASTQ is not supported)")]
    MissingSeparator {
        /// The line that should be the separator.
        line: u64,
    },
    /// The input ends inside a record.
    #[error("line {line}: truncated record, the input ends before its {missing} line")]
    Truncated {
        /// The line number the missing line would have had.
        line: u64,
        /// The first missing line of the record.
        missing: RecordLine,
    },
    /// A line is longer than [`Reader::max_line_len`]; reading stopped at the limit.
    #[error("line {line} is longer than the limit of {max_len} bytes")]
    LineTooLong {
        /// The overlong line.
        line: u64,
        /// The limit it exceeds.
        max_len: usize,
    },
    /// The quality line and the sequence line differ in length.
    #[error(
        "line {line}: quality line has {quality_len} values but the sequence has {sequence_len} \
         bases"
    )]
    QualityLength {
        /// The quality line.
        line: u64,
        /// The number of bases in the sequence line.
        sequence_len: usize,
        /// The number of bytes in the quality line.
        quality_len: usize,
    },
}

impl ParseError {
    /// Returns the 1-based line number the error refers to.
    #[must_use]
    pub fn line(&self) -> u64 {
        match self {
            Self::Io { line, .. }
            | Self::MissingHeaderPrefix { line }
            | Self::MissingName { line }
            | Self::MissingSeparator { line }
            | Self::Truncated { line, .. }
            | Self::LineTooLong { line, .. }
            | Self::QualityLength { line, .. } => *line,
        }
    }
}

/// A streaming FASTQ reader over any [`BufRead`].
///
/// Reading reuses caller-owned buffers, so a loop over a whole file allocates only when a
/// record is longer than any before it. Wrap unbuffered sources in a large
/// [`std::io::BufReader`]; reading is line by line.
#[derive(Debug)]
pub struct Reader<R> {
    inner: R,
    line: u64,
    max_line_len: usize,
    header: Vec<u8>,
    separator: Vec<u8>,
}

impl<R: BufRead> Reader<R> {
    /// Creates a reader positioned at the start of a record, with lines of at most
    /// [`DEFAULT_MAX_LINE_LEN`] bytes.
    #[must_use]
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            line: 0,
            max_line_len: DEFAULT_MAX_LINE_LEN,
            header: Vec::new(),
            separator: Vec::new(),
        }
    }

    /// Sets the longest accepted line, in bytes without the line ending.
    ///
    /// A longer line fails with [`ParseError::LineTooLong`] after at most `max_line_len + 2`
    /// bytes of it are buffered.
    #[must_use]
    pub fn with_max_line_len(mut self, max_line_len: usize) -> Self {
        self.max_line_len = max_line_len;
        self
    }

    /// Returns the longest accepted line, in bytes without the line ending.
    #[must_use]
    pub fn max_line_len(&self) -> usize {
        self.max_line_len
    }

    /// Returns the number of lines consumed so far.
    #[must_use]
    pub fn lines_read(&self) -> u64 {
        self.line
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Reads the next record, replacing the contents of `name`, `sequence` and `quality`.
    ///
    /// `name` receives the header up to the first whitespace (the description is dropped).
    /// Returns `Ok(false)` at the end of the input, when the buffers' contents are unspecified.
    ///
    /// # Errors
    ///
    /// Returns a [`ParseError`] for malformed or truncated records and for read errors. The
    /// reader should not be used after an error.
    pub fn read_record(
        &mut self,
        name: &mut Vec<u8>,
        sequence: &mut Vec<u8>,
        quality: &mut Vec<u8>,
    ) -> Result<bool, ParseError> {
        let mut lines = LineReader {
            inner: &mut self.inner,
            line: &mut self.line,
            max_len: self.max_line_len,
        };
        loop {
            if !lines.read(&mut self.header)? {
                return Ok(false);
            }
            if !self.header.is_empty() {
                break;
            }
        }
        let header_line = *lines.line;
        let Some(rest) = self.header.strip_prefix(b"@") else {
            return Err(ParseError::MissingHeaderPrefix { line: header_line });
        };
        let name_end = rest
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(rest.len());
        if name_end == 0 {
            return Err(ParseError::MissingName { line: header_line });
        }
        name.clear();
        name.extend_from_slice(&rest[..name_end]);

        if !lines.read(sequence)? {
            return Err(truncated(header_line + 1, RecordLine::Sequence));
        }
        if !lines.read(&mut self.separator)? {
            return Err(truncated(header_line + 2, RecordLine::Separator));
        }
        if !self.separator.starts_with(b"+") {
            return Err(ParseError::MissingSeparator { line: *lines.line });
        }
        if !lines.read(quality)? {
            return Err(truncated(header_line + 3, RecordLine::Quality));
        }
        if quality.len() != sequence.len() {
            return Err(ParseError::QualityLength {
                line: *lines.line,
                sequence_len: sequence.len(),
                quality_len: quality.len(),
            });
        }
        Ok(true)
    }
}

fn truncated(line: u64, missing: RecordLine) -> ParseError {
    ParseError::Truncated { line, missing }
}

/// Reads bounded lines and counts them.
struct LineReader<'reader, R> {
    inner: &'reader mut R,
    /// The number of lines read so far.
    line: &'reader mut u64,
    /// The longest accepted line, without its line ending.
    max_len: usize,
}

impl<R: BufRead> LineReader<'_, R> {
    /// Replaces `buf` with the next line, without its `\n` or `\r\n`, and counts it.
    ///
    /// Returns `Ok(false)` at the end of the input. Buffers at most `max_len + 2` bytes, enough
    /// for a line of exactly `max_len` bytes and its `\r\n`.
    fn read(&mut self, buf: &mut Vec<u8>) -> Result<bool, ParseError> {
        let this_line = *self.line + 1;
        buf.clear();
        let limit = u64::try_from(self.max_len)
            .unwrap_or(u64::MAX)
            .saturating_add(2);
        let read = (&mut *self.inner)
            .take(limit)
            .read_until(b'\n', buf)
            .map_err(|source| ParseError::Io {
                line: this_line,
                source,
            })?;
        if read == 0 {
            return Ok(false);
        }
        *self.line = this_line;
        if buf.last() == Some(&b'\n') {
            buf.pop();
            if buf.last() == Some(&b'\r') {
                buf.pop();
            }
        }
        // Without a newline within the limit, `buf` holds `max_len + 2` bytes and fails here too.
        if buf.len() > self.max_len {
            return Err(ParseError::LineTooLong {
                line: this_line,
                max_len: self.max_len,
            });
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A parsed `(name, sequence)` pair.
    type NamedSequence = (Vec<u8>, Vec<u8>);

    /// Reads every record of `data`, returning `(name, sequence)` pairs or the first error.
    fn parse(data: &[u8]) -> Result<Vec<NamedSequence>, ParseError> {
        let mut reader = Reader::new(data);
        let (mut name, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
        let mut records = Vec::new();
        while reader.read_record(&mut name, &mut seq, &mut qual)? {
            records.push((name.clone(), seq.clone()));
        }
        Ok(records)
    }

    #[test]
    fn parses_records_with_crlf_and_trailing_blank_lines() {
        let data = b"@a desc\r\nACGN\r\n+\r\nIIII\r\n@b\nacgt\n+b\n!!!!\n\n\n";
        let records = parse(data).unwrap();
        assert_eq!(
            records,
            vec![
                (b"a".to_vec(), b"ACGN".to_vec()),
                (b"b".to_vec(), b"acgt".to_vec())
            ]
        );
    }

    #[test]
    fn accepts_missing_final_newline_and_empty_reads() {
        let records = parse(b"@e\n\n+\n\n@a\nAC\n+\nII").unwrap();
        assert_eq!(records.len(), 2);
        assert!(records[0].1.is_empty());
        assert_eq!(records[1].1, b"AC");
    }

    #[test]
    fn reports_truncation_with_line_numbers() {
        let cases: [(&[u8], u64, RecordLine); 3] = [
            (b"@a\nAC\n+\nII\n@b\n", 6, RecordLine::Sequence),
            (b"@a\nAC\n", 3, RecordLine::Separator),
            (b"@a\nAC\n+\n", 4, RecordLine::Quality),
        ];
        for (data, want_line, want_missing) in cases {
            match parse(data) {
                Err(ParseError::Truncated { line, missing }) => {
                    assert_eq!((line, missing), (want_line, want_missing), "{data:?}");
                }
                other => panic!("expected truncation for {data:?}, got {other:?}"),
            }
        }
    }

    #[test]
    fn rejects_malformed_records() {
        let header = parse(b"@a\nAC\n+\nII\nb\nAC\n+\nII\n").unwrap_err();
        assert!(matches!(
            header,
            ParseError::MissingHeaderPrefix { line: 5 }
        ));

        let name = parse(b"@ desc\nAC\n+\nII\n").unwrap_err();
        assert!(matches!(name, ParseError::MissingName { line: 1 }));

        let multiline = parse(b"@a\nAC\nGT\n+\nIIII\n").unwrap_err();
        assert!(matches!(
            multiline,
            ParseError::MissingSeparator { line: 3 }
        ));

        let quality = parse(b"@a\nACG\n+\nII\n").unwrap_err();
        assert!(matches!(
            quality,
            ParseError::QualityLength {
                line: 4,
                sequence_len: 3,
                quality_len: 2
            }
        ));
        assert_eq!(quality.line(), 4);
    }

    /// Reads every record of `data` with a line limit of `max_line_len`.
    fn parse_with_limit(data: &[u8], max_line_len: usize) -> Result<usize, ParseError> {
        let mut reader = Reader::new(data).with_max_line_len(max_line_len);
        let (mut name, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
        let mut records = 0;
        while reader.read_record(&mut name, &mut seq, &mut qual)? {
            records += 1;
        }
        Ok(records)
    }

    #[test]
    fn accepts_lines_up_to_the_limit() {
        assert_eq!(parse_with_limit(b"@a\nACGT\n+\nIIII\n", 4).ok(), Some(1));
        assert_eq!(
            parse_with_limit(b"@a\r\nACGT\r\n+\r\nIIII\r\n", 4).ok(),
            Some(1)
        );
        assert_eq!(parse_with_limit(b"@a\nACGT\n+\nIIII", 4).ok(), Some(1));
    }

    #[test]
    fn stops_reading_an_overlong_line_at_the_limit() {
        let cases: [(&[u8], u64); 4] = [
            (b"@a\nACGTA\n+\nIIIII\n", 2),
            (b"@a\nACGTA", 2),
            (b"@a\nAC\n+\nII\n@abcdef\nA\n+\nI\n", 5),
            (b"@a\nACGT\r\r\n+\nIIII\n", 2),
        ];
        for (data, want_line) in cases {
            match parse_with_limit(data, 4) {
                Err(ParseError::LineTooLong { line, max_len: 4 }) => {
                    assert_eq!(line, want_line, "{data:?}");
                }
                other => panic!("expected an overlong line for {data:?}, got {other:?}"),
            }
        }

        // Reading stops at the limit: the rest of a huge line stays unread.
        let mut data = b"@a\n".to_vec();
        data.extend(std::iter::repeat_n(b'A', 10_000));
        let mut source = &data[..];
        let mut reader = Reader::new(&mut source).with_max_line_len(100);
        let (mut name, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
        let error = reader
            .read_record(&mut name, &mut seq, &mut qual)
            .unwrap_err();
        assert!(
            matches!(error, ParseError::LineTooLong { line: 2, .. }),
            "{error:?}"
        );
        assert!(seq.len() <= 102, "buffered {} bytes", seq.len());
        drop(reader);
        assert_eq!(source.len(), data.len() - 3 - 102);
    }
}
