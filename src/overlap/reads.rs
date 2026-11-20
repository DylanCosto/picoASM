//! The reads an overlap assembly works on, held together in one buffer.

use std::ops::Range;

/// A read's index in a [`Reads`] set.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct ReadIndex(pub u32);

impl ReadIndex {
    /// Returns the index as a `usize`.
    #[must_use]
    pub fn index(self) -> usize {
        usize::try_from(self.0).expect("read indices fit in usize")
    }
}

/// A set of reads, their bases in one buffer, read by [`ReadIndex`].
#[derive(Clone, Debug)]
pub struct Reads {
    bases: Vec<u8>,
    /// Read `i`'s bases are `bases[starts[i]..starts[i + 1]]`.
    starts: Vec<usize>,
}

impl Reads {
    /// Returns an empty set.
    #[must_use]
    pub fn new() -> Self {
        Self {
            bases: Vec::new(),
            starts: vec![0],
        }
    }

    /// Returns an empty set with room for `reads` reads of `bases` bases in all.
    #[must_use]
    pub fn with_capacity(reads: usize, bases: usize) -> Self {
        let mut starts = Vec::with_capacity(reads + 1);
        starts.push(0);
        Self {
            bases: Vec::with_capacity(bases),
            starts,
        }
    }

    /// Frees the room held beyond the reads added (a set grown read by read holds up to twice
    /// its bases).
    pub fn shrink_to_fit(&mut self) {
        self.bases.shrink_to_fit();
        self.starts.shrink_to_fit();
    }

    /// Adds a read (its bases upper-cased, any other than `ACGT` as `N`) and returns its index.
    ///
    /// # Panics
    ///
    /// Panics if the set already holds `u32::MAX` reads, or the read is longer than `u32::MAX`
    /// bases.
    pub fn push(&mut self, read: &[u8]) -> ReadIndex {
        assert!(
            u32::try_from(read.len()).is_ok(),
            "read of {} bases is too long",
            read.len()
        );
        let index = ReadIndex(u32::try_from(self.len()).expect("fewer than u32::MAX reads"));
        self.bases.extend(read.iter().map(|&b| match b {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            _ => b'N',
        }));
        self.starts.push(self.bases.len());
        index
    }

    /// Returns the number of reads.
    #[must_use]
    pub fn len(&self) -> usize {
        self.starts.len().saturating_sub(1)
    }

    /// Returns whether the set holds no reads.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns a read's bases.
    ///
    /// # Panics
    ///
    /// Panics if `read` is not in the set.
    #[must_use]
    pub fn get(&self, read: ReadIndex) -> &[u8] {
        &self.bases[self.range(read)]
    }

    /// Returns a read's length.
    ///
    /// # Panics
    ///
    /// Panics if `read` is not in the set.
    #[must_use]
    pub fn read_len(&self, read: ReadIndex) -> u32 {
        u32::try_from(self.range(read).len()).expect("checked when the read was added")
    }

    /// Returns the total number of bases.
    #[must_use]
    pub fn total_bases(&self) -> usize {
        self.bases.len()
    }

    /// Returns the indices of all reads, in order.
    pub fn indices(&self) -> impl Iterator<Item = ReadIndex> + use<> {
        (0..u32::try_from(self.len()).expect("fewer than u32::MAX reads")).map(ReadIndex)
    }

    fn range(&self, read: ReadIndex) -> Range<usize> {
        let i = read.index();
        self.starts[i]..self.starts[i + 1]
    }
}

impl Default for Reads {
    fn default() -> Self {
        Self::new()
    }
}

impl<'r> FromIterator<&'r [u8]> for Reads {
    fn from_iter<I: IntoIterator<Item = &'r [u8]>>(iter: I) -> Self {
        let mut reads = Self::new();
        for read in iter {
            reads.push(read);
        }
        reads
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_come_back_upper_cased_in_order_with_other_bases_as_n() {
        let reads: Reads = [&b"acgt"[..], b"", b"GRa"].into_iter().collect();
        assert_eq!(reads.len(), 3);
        assert_eq!(reads.get(ReadIndex(0)), b"ACGT");
        assert_eq!(reads.get(ReadIndex(1)), b"");
        assert_eq!(reads.get(ReadIndex(2)), b"GNA");
        assert_eq!(reads.read_len(ReadIndex(2)), 3);
        assert_eq!(reads.total_bases(), 7);
        assert_eq!(reads.indices().collect::<Vec<_>>().len(), 3);
    }

    #[test]
    fn room_given_or_freed_keeps_the_reads() {
        let mut reads = Reads::with_capacity(2, 100);
        reads.push(b"ACGT");
        reads.push(b"GG");
        reads.shrink_to_fit();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads.get(ReadIndex(1)), b"GG");
        assert_eq!(reads.total_bases(), 6);
    }
}
