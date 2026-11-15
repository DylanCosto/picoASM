//! Read identifiers.

use std::fmt;

/// The 0-based index of a read in input order, across all input files.
///
/// Every pass over the input assigns the same id to the same read.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct ReadId(pub u32);

impl ReadId {
    /// Returns the id as an index into a per-read table.
    #[must_use]
    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<u32> for ReadId {
    fn from(id: u32) -> Self {
        Self(id)
    }
}

impl From<ReadId> for u32 {
    fn from(id: ReadId) -> Self {
        id.0
    }
}

impl fmt::Display for ReadId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
