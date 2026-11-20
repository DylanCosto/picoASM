//! Errors from the FASTA and GFA writers.

use std::{
    fmt, io,
    path::{Path, PathBuf},
};

/// An error from [`FastaWriter`](crate::io::FastaWriter) or [`GfaWriter`](crate::io::GfaWriter).
///
/// Validation errors are reported before anything of the offending record or line is written,
/// so the output stays well-formed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum WriteError {
    /// The output file could not be created.
    #[error("cannot create {}", path.display())]
    Create {
        /// The output file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// Writing or flushing failed.
    #[error("cannot write {}", OutputName(path.as_deref()))]
    Io {
        /// The output file, when the writer was made by `create`; `None` for a writer over a
        /// caller-supplied stream.
        path: Option<PathBuf>,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A record, segment or path name is not allowed by the format.
    #[error("invalid name {name:?}: {reason}")]
    InvalidName {
        /// The name, lossily converted to UTF-8.
        name: String,
        /// Why it was rejected.
        reason: &'static str,
    },
    /// A sequence holds bytes the format cannot represent.
    #[error("invalid sequence for {name:?}: {reason}")]
    InvalidSequence {
        /// The name of the record or segment.
        name: String,
        /// Why it was rejected.
        reason: &'static str,
    },
    /// A GFA tag has an invalid name or value.
    #[error("invalid GFA tag {tag:?}: {reason}")]
    InvalidTag {
        /// The tag, roughly as it would have been written.
        tag: String,
        /// Why it was rejected.
        reason: &'static str,
    },
    /// A GFA overlap is not a CIGAR string.
    #[error("invalid GFA overlap {overlap:?}: {reason}")]
    InvalidOverlap {
        /// The overlap text.
        overlap: String,
        /// Why it was rejected.
        reason: &'static str,
    },
    /// A GFA path has overlaps, but not one per pair of consecutive steps.
    #[error(
        "GFA path {name:?} has {steps} steps but {overlaps} overlaps (expected 0 or steps - 1)"
    )]
    PathOverlapCount {
        /// The path name.
        name: String,
        /// The number of steps.
        steps: usize,
        /// The number of overlaps given.
        overlaps: usize,
    },
    /// A GFA path has no steps.
    #[error("GFA path {name:?} has no steps")]
    EmptyPath {
        /// The path name.
        name: String,
    },
}

/// Converts an I/O error from an output with no known path.
impl From<io::Error> for WriteError {
    fn from(source: io::Error) -> Self {
        Self::Io { path: None, source }
    }
}

/// Returns a function that wraps an I/O error from the output at `path`.
pub(crate) fn io_error(path: Option<&Path>) -> impl FnOnce(io::Error) -> WriteError + '_ {
    move |source| WriteError::Io {
        path: path.map(Path::to_path_buf),
        source,
    }
}

/// Displays an output path, or a placeholder when there is none.
struct OutputName<'path>(Option<&'path Path>);

impl fmt::Display for OutputName<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Some(path) => path.display().fmt(f),
            None => f.write_str("the output"),
        }
    }
}
