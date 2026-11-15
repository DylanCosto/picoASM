//! GFA 1 output from plain borrowed data.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use picoasm_seq::Strand;

use crate::{WriteError, write_error::io_error};

/// The header tag this writer always writes: the GFA version.
const VERSION_TAG: [u8; 2] = *b"VN";

/// The value of a GFA optional field.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum TagValue<'data> {
    /// A signed integer, written with type `i`.
    Int(i64),
    /// An unsigned integer, written with type `i`.
    UInt(u64),
    /// A finite float, written with type `f` in Rust's shortest round-trip form.
    Float(f64),
    /// Printable ASCII text (spaces allowed), written with type `Z`.
    String(&'data str),
}

/// A GFA optional field such as `LN:i:1500`.
///
/// The name and value are checked when the tag is written, so building a tag never fails.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GfaTag<'data> {
    name: [u8; 2],
    value: TagValue<'data>,
}

impl<'data> GfaTag<'data> {
    /// Creates a tag. A valid name is a letter followed by a letter or digit.
    #[must_use]
    pub const fn new(name: [u8; 2], value: TagValue<'data>) -> Self {
        Self { name, value }
    }

    /// Returns an `LN:i` tag: the segment length in bases.
    #[must_use]
    pub const fn length(bases: u64) -> Self {
        Self::new(*b"LN", TagValue::UInt(bases))
    }

    /// Returns a `KC:i` tag: the k-mer count (total coverage) of a segment.
    #[must_use]
    pub const fn kmer_count(count: u64) -> Self {
        Self::new(*b"KC", TagValue::UInt(count))
    }

    /// Returns an `RC:i` tag: the read count of a segment or link.
    #[must_use]
    pub const fn read_count(count: u64) -> Self {
        Self::new(*b"RC", TagValue::UInt(count))
    }

    /// Returns a `dp:f` tag: the mean depth of a segment.
    #[must_use]
    pub const fn depth(depth: f64) -> Self {
        Self::new(*b"dp", TagValue::Float(depth))
    }

    /// Returns the two-byte tag name.
    #[must_use]
    pub const fn name(&self) -> [u8; 2] {
        self.name
    }

    /// Returns the tag value.
    #[must_use]
    pub const fn value(&self) -> TagValue<'data> {
        self.value
    }
}

/// A segment name with an orientation, as used in links and path steps.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct OrientedSegment<'data> {
    /// The segment name.
    pub name: &'data str,
    /// The orientation: [`Strand::Forward`] is `+`, [`Strand::Reverse`] is `-`.
    pub strand: Strand,
}

impl<'data> OrientedSegment<'data> {
    /// Creates an oriented segment.
    #[must_use]
    pub const fn new(name: &'data str, strand: Strand) -> Self {
        Self { name, strand }
    }
}

/// The overlap of a link or of consecutive path steps.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Overlap<'data> {
    /// Unknown, written as `*`.
    #[default]
    Unknown,
    /// An exact overlap of this many bases, written as `<n>M`.
    Matches(u64),
    /// A CIGAR string such as `10M1I5M`, written as is after validation.
    Cigar(&'data str),
}

/// A GFA `S` line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GfaSegment<'data> {
    /// The segment name.
    pub name: &'data str,
    /// The bases, or `None` to write `*` (an empty sequence is also written as `*`).
    pub sequence: Option<&'data [u8]>,
    /// Optional fields, written in order.
    pub tags: &'data [GfaTag<'data>],
}

impl<'data> GfaSegment<'data> {
    /// Creates a segment with no sequence and no tags.
    #[must_use]
    pub const fn new(name: &'data str) -> Self {
        Self {
            name,
            sequence: None,
            tags: &[],
        }
    }

    /// Sets the sequence.
    #[must_use]
    pub const fn with_sequence(mut self, sequence: &'data [u8]) -> Self {
        self.sequence = Some(sequence);
        self
    }

    /// Sets the optional fields.
    #[must_use]
    pub const fn with_tags(mut self, tags: &'data [GfaTag<'data>]) -> Self {
        self.tags = tags;
        self
    }
}

/// A GFA `L` line: `from` is followed by `to`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GfaLink<'data> {
    /// The source segment and its orientation.
    pub from: OrientedSegment<'data>,
    /// The target segment and its orientation.
    pub to: OrientedSegment<'data>,
    /// The overlap between the end of `from` and the start of `to`.
    pub overlap: Overlap<'data>,
    /// Optional fields, written in order.
    pub tags: &'data [GfaTag<'data>],
}

impl<'data> GfaLink<'data> {
    /// Creates a link with no tags.
    #[must_use]
    pub const fn new(
        from: OrientedSegment<'data>,
        to: OrientedSegment<'data>,
        overlap: Overlap<'data>,
    ) -> Self {
        Self {
            from,
            to,
            overlap,
            tags: &[],
        }
    }

    /// Sets the optional fields.
    #[must_use]
    pub const fn with_tags(mut self, tags: &'data [GfaTag<'data>]) -> Self {
        self.tags = tags;
        self
    }
}

/// A GFA `P` line: a named walk through oriented segments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GfaPath<'data> {
    /// The path name.
    pub name: &'data str,
    /// The steps, at least one.
    pub steps: &'data [OrientedSegment<'data>],
    /// Either empty (written as `*`) or one overlap per pair of consecutive steps.
    pub overlaps: &'data [Overlap<'data>],
    /// Optional fields, written in order.
    pub tags: &'data [GfaTag<'data>],
}

impl<'data> GfaPath<'data> {
    /// Creates a path with unknown overlaps and no tags.
    #[must_use]
    pub const fn new(name: &'data str, steps: &'data [OrientedSegment<'data>]) -> Self {
        Self {
            name,
            steps,
            overlaps: &[],
            tags: &[],
        }
    }

    /// Sets the overlaps between consecutive steps.
    #[must_use]
    pub const fn with_overlaps(mut self, overlaps: &'data [Overlap<'data>]) -> Self {
        self.overlaps = overlaps;
        self
    }

    /// Sets the optional fields.
    #[must_use]
    pub const fn with_tags(mut self, tags: &'data [GfaTag<'data>]) -> Self {
        self.tags = tags;
        self
    }
}

/// Writes GFA 1: a header, then `S`, `L` and `P` lines in the caller's order.
///
/// Every line is validated before any of it is written. Names must be printable ASCII without
/// spaces or commas and must not start with `*` or `=`; sequences must be letters, `=` or `.`;
/// overlaps must be CIGAR strings. The writer does no buffering of its own beyond one line:
/// wrap unbuffered outputs in a [`BufWriter`], as [`GfaWriter::create`] does, and call
/// [`GfaWriter::finish`] to flush.
///
/// # Examples
///
/// ```
/// use picoasm_io::{GfaLink, GfaSegment, GfaTag, GfaWriter, OrientedSegment, Overlap};
/// use picoasm_seq::Strand;
///
/// # fn main() -> Result<(), picoasm_io::WriteError> {
/// let mut gfa = GfaWriter::new(Vec::new());
/// gfa.write_header(&[])?;
/// gfa.write_segment(&GfaSegment::new("u1").with_tags(&[GfaTag::length(1200)]))?;
/// gfa.write_segment(&GfaSegment::new("u2").with_sequence(b"ACGT"))?;
/// let from = OrientedSegment::new("u1", Strand::Forward);
/// let to = OrientedSegment::new("u2", Strand::Reverse);
/// gfa.write_link(&GfaLink::new(from, to, Overlap::Matches(3)))?;
/// let text = gfa.finish()?;
/// assert_eq!(
///     text,
///     b"H\tVN:Z:1.0\nS\tu1\t*\tLN:i:1200\nS\tu2\tACGT\nL\tu1\t+\tu2\t-\t3M\n"
/// );
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct GfaWriter<W> {
    inner: W,
    line: Vec<u8>,
    /// The output file, named in I/O errors, when made by [`GfaWriter::create`].
    path: Option<PathBuf>,
}

impl GfaWriter<BufWriter<File>> {
    /// Creates (or truncates) the file at `path` and returns a buffered writer.
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

impl<W: Write> GfaWriter<W> {
    /// Creates a writer.
    #[must_use]
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            line: Vec::new(),
            path: None,
        }
    }

    /// Writes the header line `H\tVN:Z:1.0`, followed by `tags`.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::InvalidTag`] for an invalid tag or a `VN` tag (the writer writes
    /// the version itself) and [`WriteError::Io`] when writing fails.
    pub fn write_header(&mut self, tags: &[GfaTag<'_>]) -> Result<(), WriteError> {
        if let Some(tag) = tags.iter().find(|tag| tag.name == VERSION_TAG) {
            return Err(WriteError::InvalidTag {
                tag: describe_tag(tag),
                reason: "the header's VN tag is written by the writer (GFA 1.0)",
            });
        }
        self.line.clear();
        self.line.extend_from_slice(b"H\tVN:Z:1.0");
        push_tags(&mut self.line, tags)?;
        self.line.push(b'\n');
        self.emit_line()
    }

    /// Writes an `S` line.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::InvalidName`], [`WriteError::InvalidSequence`] or
    /// [`WriteError::InvalidTag`] for invalid input and [`WriteError::Io`] when writing fails.
    pub fn write_segment(&mut self, segment: &GfaSegment<'_>) -> Result<(), WriteError> {
        validate_name(segment.name)?;
        let sequence = match segment.sequence {
            None | Some([]) => &b"*"[..],
            Some(sequence) => {
                if !sequence
                    .iter()
                    .all(|&b| b.is_ascii_alphabetic() || b == b'=' || b == b'.')
                {
                    return Err(WriteError::InvalidSequence {
                        name: segment.name.to_owned(),
                        reason: "a GFA sequence must hold only letters, '=' or '.'",
                    });
                }
                sequence
            }
        };
        // The sequence can be megabases long, so it is written directly rather than copied
        // into the line buffer: `line` holds the text before it, then the text after it.
        self.line.clear();
        self.line.extend_from_slice(b"S\t");
        self.line.extend_from_slice(segment.name.as_bytes());
        self.line.push(b'\t');
        let prefix_len = self.line.len();
        push_tags(&mut self.line, segment.tags)?;
        self.line.push(b'\n');
        let (prefix, suffix) = self.line.split_at(prefix_len);
        write_parts(&mut self.inner, &[prefix, sequence, suffix])
            .map_err(io_error(self.path.as_deref()))
    }

    /// Writes an `L` line.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::InvalidName`], [`WriteError::InvalidOverlap`] or
    /// [`WriteError::InvalidTag`] for invalid input and [`WriteError::Io`] when writing fails.
    pub fn write_link(&mut self, link: &GfaLink<'_>) -> Result<(), WriteError> {
        validate_name(link.from.name)?;
        validate_name(link.to.name)?;
        self.line.clear();
        self.line.extend_from_slice(b"L\t");
        for (i, end) in [link.from, link.to].into_iter().enumerate() {
            if i > 0 {
                self.line.push(b'\t');
            }
            self.line.extend_from_slice(end.name.as_bytes());
            self.line.push(b'\t');
            self.line.push(strand_char(end.strand));
        }
        self.line.push(b'\t');
        push_overlap(&mut self.line, link.overlap)?;
        push_tags(&mut self.line, link.tags)?;
        self.line.push(b'\n');
        self.emit_line()
    }

    /// Writes a `P` line.
    ///
    /// # Errors
    ///
    /// Returns [`WriteError::EmptyPath`] for a path without steps,
    /// [`WriteError::PathOverlapCount`] when overlaps are given but not one per pair of steps,
    /// [`WriteError::InvalidOverlap`] for an invalid or [`Overlap::Unknown`] overlap in a
    /// non-empty overlap list, [`WriteError::InvalidName`] or [`WriteError::InvalidTag`] for
    /// other invalid input, and [`WriteError::Io`] when writing fails.
    pub fn write_path(&mut self, path: &GfaPath<'_>) -> Result<(), WriteError> {
        validate_name(path.name)?;
        if path.steps.is_empty() {
            return Err(WriteError::EmptyPath {
                name: path.name.to_owned(),
            });
        }
        if !path.overlaps.is_empty() && path.overlaps.len() + 1 != path.steps.len() {
            return Err(WriteError::PathOverlapCount {
                name: path.name.to_owned(),
                steps: path.steps.len(),
                overlaps: path.overlaps.len(),
            });
        }
        self.line.clear();
        self.line.extend_from_slice(b"P\t");
        self.line.extend_from_slice(path.name.as_bytes());
        self.line.push(b'\t');
        for (i, step) in path.steps.iter().enumerate() {
            validate_name(step.name)?;
            if i > 0 {
                self.line.push(b',');
            }
            self.line.extend_from_slice(step.name.as_bytes());
            self.line.push(strand_char(step.strand));
        }
        self.line.push(b'\t');
        if path.overlaps.is_empty() {
            self.line.push(b'*');
        }
        for (i, &overlap) in path.overlaps.iter().enumerate() {
            if overlap == Overlap::Unknown {
                return Err(WriteError::InvalidOverlap {
                    overlap: "*".to_owned(),
                    reason: "'*' cannot appear in a list of path overlaps",
                });
            }
            if i > 0 {
                self.line.push(b',');
            }
            push_overlap(&mut self.line, overlap)?;
        }
        push_tags(&mut self.line, path.tags)?;
        self.line.push(b'\n');
        self.emit_line()
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

    /// Writes the line buffer, which holds one complete line.
    fn emit_line(&mut self) -> Result<(), WriteError> {
        write_parts(&mut self.inner, &[&self.line]).map_err(io_error(self.path.as_deref()))
    }
}

/// Writes `parts` in order.
fn write_parts<W: Write>(out: &mut W, parts: &[&[u8]]) -> io::Result<()> {
    parts.iter().try_for_each(|part| out.write_all(part))
}

/// Returns a tag roughly as it would be written, for error messages.
fn describe_tag(tag: &GfaTag<'_>) -> String {
    format!("{}:{:?}", String::from_utf8_lossy(&tag.name), tag.value)
}

const fn strand_char(strand: Strand) -> u8 {
    match strand {
        Strand::Forward => b'+',
        Strand::Reverse => b'-',
    }
}

/// Checks a segment or path name: printable ASCII without spaces or commas, not starting with
/// `*` or `=` (GFA 1 name grammar, minus commas so path step lists stay unambiguous).
fn validate_name(name: &str) -> Result<(), WriteError> {
    let reason = match name.as_bytes() {
        [] => "a GFA name must not be empty",
        [b'*' | b'=', ..] => "a GFA name must not start with '*' or '='",
        bytes if !bytes.iter().all(u8::is_ascii_graphic) => {
            "a GFA name must be printable ASCII without spaces"
        }
        bytes if bytes.contains(&b',') => "a GFA name must not contain ','",
        _ => return Ok(()),
    };
    Err(WriteError::InvalidName {
        name: name.to_owned(),
        reason,
    })
}

/// Appends the overlap text after validating it.
fn push_overlap(line: &mut Vec<u8>, overlap: Overlap<'_>) -> Result<(), WriteError> {
    match overlap {
        Overlap::Unknown => line.push(b'*'),
        Overlap::Matches(bases) => write!(line, "{bases}M")?,
        Overlap::Cigar(cigar) => {
            if !is_cigar(cigar.as_bytes()) {
                return Err(WriteError::InvalidOverlap {
                    overlap: cigar.to_owned(),
                    reason: "expected a CIGAR string such as 10M2I5M",
                });
            }
            line.extend_from_slice(cigar.as_bytes());
        }
    }
    Ok(())
}

/// Returns `true` for `([0-9]+[MIDNSHPX=])+`.
fn is_cigar(cigar: &[u8]) -> bool {
    let mut digits = 0_usize;
    for &b in cigar {
        if b.is_ascii_digit() {
            digits += 1;
        } else if digits > 0 && b"MIDNSHPX=".contains(&b) {
            digits = 0;
        } else {
            return false;
        }
    }
    !cigar.is_empty() && digits == 0
}

/// Appends `\t`-separated tags after validating each.
fn push_tags(line: &mut Vec<u8>, tags: &[GfaTag<'_>]) -> Result<(), WriteError> {
    for tag in tags {
        let [first, second] = tag.name;
        let invalid = |reason| WriteError::InvalidTag {
            tag: describe_tag(tag),
            reason,
        };
        if !first.is_ascii_alphabetic() || !second.is_ascii_alphanumeric() {
            return Err(invalid(
                "a tag name is a letter followed by a letter or digit",
            ));
        }
        line.push(b'\t');
        line.extend_from_slice(&tag.name);
        match tag.value {
            TagValue::Int(value) => write!(line, ":i:{value}")?,
            TagValue::UInt(value) => write!(line, ":i:{value}")?,
            TagValue::Float(value) => {
                if !value.is_finite() {
                    return Err(invalid("a float tag must be finite"));
                }
                write!(line, ":f:{value}")?;
            }
            TagValue::String(text) => {
                if text.is_empty() || !text.bytes().all(|b| b == b' ' || b.is_ascii_graphic()) {
                    return Err(invalid("a string tag must be non-empty printable ASCII"));
                }
                line.extend_from_slice(b":Z:");
                line.extend_from_slice(text.as_bytes());
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_names() {
        for name in ["utg000001l", "a+b", "x.1|2"] {
            assert!(validate_name(name).is_ok(), "{name}");
        }
        for name in ["", "*a", "=a", "a b", "a\tb", "a,b", "é"] {
            assert!(validate_name(name).is_err(), "{name}");
        }
    }

    #[test]
    fn validates_cigars() {
        for cigar in ["0M", "10M2I5M", "3=1X"] {
            assert!(is_cigar(cigar.as_bytes()), "{cigar}");
        }
        for cigar in ["", "M", "10", "10Q", "*", "5M3"] {
            assert!(!is_cigar(cigar.as_bytes()), "{cigar}");
        }
    }

    #[test]
    fn writes_tags_and_rejects_bad_ones() {
        let mut line = Vec::new();
        let tags = [
            GfaTag::length(10),
            GfaTag::kmer_count(7),
            GfaTag::depth(2.5),
            GfaTag::new(*b"xi", TagValue::Int(-3)),
            GfaTag::new(*b"zz", TagValue::String("a b")),
        ];
        push_tags(&mut line, &tags).unwrap();
        assert_eq!(line, b"\tLN:i:10\tKC:i:7\tdp:f:2.5\txi:i:-3\tzz:Z:a b");

        let mut gfa = GfaWriter::new(Vec::new());
        let version = GfaTag::new(*b"VN", TagValue::String("1.0"));
        assert!(matches!(
            gfa.write_header(&[version]),
            Err(WriteError::InvalidTag { .. })
        ));
        assert!(gfa.get_ref().is_empty());

        for bad in [
            GfaTag::new(*b"1a", TagValue::Int(0)),
            GfaTag::new(*b"a_", TagValue::Int(0)),
            GfaTag::depth(f64::NAN),
            GfaTag::new(*b"zz", TagValue::String("")),
            GfaTag::new(*b"zz", TagValue::String("a\tb")),
        ] {
            assert!(push_tags(&mut Vec::new(), &[bad]).is_err(), "{bad:?}");
        }
    }
}
