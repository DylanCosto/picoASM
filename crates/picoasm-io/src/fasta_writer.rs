//! FASTA output.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use crate::{WriteError, write_error::io_error};

/// The default FASTA line width.
pub const DEFAULT_LINE_WIDTH: usize = 80;

/// Writes FASTA records with sequences wrapped at a fixed line width.
///
/// The writer does no buffering of its own: wrap unbuffered outputs in a [`BufWriter`], as
/// [`FastaWriter::create`] does. Call [`FastaWriter::finish`] to flush and see flush errors;
/// dropping the writer flushes only what the inner writer flushes on drop.
///
/// # Examples
///
/// ```
/// use picoasm_io::FastaWriter;
///
/// # fn main() -> Result<(), picoasm_io::WriteError> {
/// let mut writer = FastaWriter::with_line_width(Vec::new(), 4);
/// writer.write_record("ctg1 len=6", "ACGTAC")?;
/// assert_eq!(writer.finish()?, b">ctg1 len=6\nACGT\nAC\n");
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct FastaWriter<W> {
    inner: W,
    line_width: usize,
    /// The output file, named in I/O errors, when made by [`FastaWriter::create`].
    path: Option<PathBuf>,
}

impl FastaWriter<BufWriter<File>> {
    /// Creates (or truncates) the file at `path` and returns a buffered writer with the default
    /// line width.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::Create`] when the file cannot be created.
    pub fn create(path: impl AsRef<Path>) -> Result<Self, WriteError> {
        let path = path.as_ref();
        let file = File::create(path).map_err(|source| WriteError::Create {
            path: path.to_path_buf(),
            source,
        })?;
        let mut writer = Self::new(BufWriter::new(file));
        writer.path = Some(path.to_path_buf());
        Ok(writer)
    }
}

impl<W: Write> FastaWriter<W> {
    /// Creates a writer with lines of [`DEFAULT_LINE_WIDTH`] bases.
    #[must_use]
    pub fn new(inner: W) -> Self {
        Self::with_line_width(inner, DEFAULT_LINE_WIDTH)
    }

    /// Creates a writer with lines of `line_width` bases; `0` writes each sequence on one line.
    #[must_use]
    pub fn with_line_width(inner: W, line_width: usize) -> Self {
        Self {
            inner,
            line_width,
            path: None,
        }
    }

    /// Returns the line width; `0` means unwrapped.
    #[must_use]
    pub fn line_width(&self) -> usize {
        self.line_width
    }

    /// Writes one record: `>` and `header`, then `sequence` wrapped at the line width.
    ///
    /// `header` is the whole definition line after `>`: the name, optionally followed by a space
    /// and a description. An empty sequence writes only the header line.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::InvalidName`] for an empty header or one that starts with
    /// whitespace or contains a line break, [`WriteError::InvalidSequence`] for a sequence that
    /// contains a line break or `>`, and [`WriteError::Io`] when writing fails. Nothing is
    /// written for an invalid record.
    pub fn write_record(
        &mut self,
        header: impl AsRef<[u8]>,
        sequence: impl AsRef<[u8]>,
    ) -> Result<(), WriteError> {
        let header = header.as_ref();
        let sequence = sequence.as_ref();
        validate_header(header)?;
        if memchr::memchr3(b'\n', b'\r', b'>', sequence).is_some() {
            return Err(WriteError::InvalidSequence {
                name: String::from_utf8_lossy(header).into_owned(),
                reason: "a FASTA sequence must not contain line breaks or '>'",
            });
        }

        write_wrapped(&mut self.inner, self.line_width, header, sequence)
            .map_err(io_error(self.path.as_deref()))
    }

    /// Returns a reference to the inner writer.
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the inner writer.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Flushes the inner writer and returns it.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::Io`] when flushing fails.
    pub fn finish(mut self) -> Result<W, WriteError> {
        self.inner.flush().map_err(io_error(self.path.as_deref()))?;
        Ok(self.inner)
    }
}

/// Writes a validated record, wrapping `sequence` at `line_width` bases (`0`: unwrapped).
fn write_wrapped<W: Write>(
    out: &mut W,
    line_width: usize,
    header: &[u8],
    sequence: &[u8],
) -> io::Result<()> {
    out.write_all(b">")?;
    out.write_all(header)?;
    out.write_all(b"\n")?;
    if sequence.is_empty() {
        return Ok(());
    }
    let width = if line_width == 0 {
        sequence.len()
    } else {
        line_width
    };
    for line in sequence.chunks(width) {
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn validate_header(header: &[u8]) -> Result<(), WriteError> {
    let reason = if header.is_empty() {
        "a FASTA header must not be empty"
    } else if header[0].is_ascii_whitespace() {
        "a FASTA header must start with the name, not whitespace"
    } else if memchr::memchr2(b'\n', b'\r', header).is_some() {
        "a FASTA header must not contain line breaks"
    } else {
        return Ok(());
    };
    Err(WriteError::InvalidName {
        name: String::from_utf8_lossy(header).into_owned(),
        reason,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(width: usize, records: &[(&str, &str)]) -> Vec<u8> {
        let mut writer = FastaWriter::with_line_width(Vec::new(), width);
        for (name, seq) in records {
            writer.write_record(name, seq).unwrap();
        }
        writer.finish().unwrap()
    }

    #[test]
    fn wraps_at_line_width() {
        assert_eq!(
            write(3, &[("a", "ACGTACG"), ("b", "ACG"), ("e", "")]),
            b">a\nACG\nTAC\nG\n>b\nACG\n>e\n"
        );
        assert_eq!(write(0, &[("a", "ACGTACG")]), b">a\nACGTACG\n");
        let long = "A".repeat(200);
        let text = write(DEFAULT_LINE_WIDTH, &[("x", &long)]);
        let lengths: Vec<usize> = text.split(|&b| b == b'\n').map(<[u8]>::len).collect();
        assert_eq!(lengths, [2, 80, 80, 40, 0]);
    }

    #[test]
    fn rejects_invalid_records_without_writing() {
        let mut writer = FastaWriter::new(Vec::new());
        for (name, seq) in [
            ("", "A"),
            (" a", "A"),
            ("a\nb", "A"),
            ("a", "AC\nGT"),
            ("a", ">"),
        ] {
            assert!(writer.write_record(name, seq).is_err(), "{name:?} {seq:?}");
        }
        assert!(writer.get_ref().is_empty());
    }
}
