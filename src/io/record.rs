//! One input read.

use crate::seq::ReadId;

/// One input read: its id, name and bases.
///
/// Records produced by a [`ReadStream`](crate::io::ReadStream) have `seq.len() <= u32::MAX` and ids
/// assigned in input order. The fields are public so callers can build records for tests or
/// move the buffers out with [`std::mem::take`].
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct ReadRecord {
    /// The 0-based index of the read in input order across all input files.
    pub id: ReadId,
    /// The read name: the first word of the FASTA/FASTQ header, or the BAM `QNAME`.
    ///
    /// Empty when a BAM record has no name (`*`).
    pub name: Vec<u8>,
    /// The bases, exactly as stored in the input (any case, non-`ACGT` bytes kept).
    pub seq: Vec<u8>,
}

impl ReadRecord {
    /// Creates a record from its parts.
    #[must_use]
    pub fn new(id: ReadId, name: Vec<u8>, seq: Vec<u8>) -> Self {
        Self { id, name, seq }
    }

    /// Clears the name and sequence, keeping their allocations, and resets the id to 0.
    pub fn clear(&mut self) {
        self.id = ReadId::default();
        self.name.clear();
        self.seq.clear();
    }
}
