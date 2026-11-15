//! Reusable batches of consecutive reads.

use std::{fmt, slice};

use picoasm_seq::ReadId;

use crate::ReadRecord;

/// The default base budget of a [`ReadBatch`]: 64 Mbp.
///
/// Large enough to keep a rayon pool busy on HiFi reads (thousands of reads per batch), small
/// enough that a batch and its per-read results stay well under a gigabyte.
pub const DEFAULT_BATCH_BASES: usize = 64_000_000;

/// The most records a [`ReadBatch`] holds, whatever its base budget.
///
/// Bounds the batch when the input holds very many short or empty reads.
pub const MAX_BATCH_RECORDS: usize = 1 << 20;

/// A record buffer keeps at most this many times the capacity its current content needs.
///
/// Bounds how far a reused slot's capacity can drift above its read: without it, every slot
/// would end up as large as the longest read it ever held. Lower values save memory but
/// reallocate more often.
const MAX_SLACK: usize = 4;

/// Buffers up to this capacity are never shrunk, so short names and reads do not reallocate.
const MIN_RETAINED_CAPACITY: usize = 256;

/// Record tables up to this many slots are never shrunk.
const MIN_RETAINED_SLOTS: usize = 1024;

/// A reusable batch of reads with contiguous, increasing ids.
///
/// Filled by [`ReadStream::next_batch`](crate::ReadStream::next_batch). Record `i` of a batch
/// has id `first_id + i`, so a caller can process [`records`](Self::records) in parallel (for
/// example with rayon's `par_iter` on the slice) and still emit results in id order.
///
/// Refilling a batch reuses the name and sequence buffers of its records, so a long-lived batch
/// performs few allocations once warmed up. Memory still follows the contents: a buffer whose
/// capacity exceeds four times what its read needs (and 256 bytes) is shrunk when the read is
/// added, and [`ReadStream::next_batch`](crate::ReadStream::next_batch) frees the buffers of
/// slots the new batch does not use. [`heap_bytes`](Self::heap_bytes) reports the result,
/// which is what a memory budget should charge for the batch.
///
/// Cloning copies only the current reads, never spare buffers.
pub struct ReadBatch {
    /// Record slots; only the first `len` hold the current batch, the rest are spare buffers.
    pub(crate) slots: Vec<ReadRecord>,
    pub(crate) len: usize,
    pub(crate) bases: usize,
    max_bases: usize,
}

impl ReadBatch {
    /// Creates an empty batch with the default base budget, [`DEFAULT_BATCH_BASES`].
    #[must_use]
    pub fn new() -> Self {
        Self::with_max_bases(DEFAULT_BATCH_BASES)
    }

    /// Creates an empty batch that is filled until it holds at least `max_bases` bases.
    ///
    /// A filled batch always holds at least one read, so `0` gives one read per batch.
    #[must_use]
    pub fn with_max_bases(max_bases: usize) -> Self {
        Self {
            slots: Vec::new(),
            len: 0,
            bases: 0,
            max_bases,
        }
    }

    /// Returns the base budget: filling stops once the batch holds at least this many bases.
    #[must_use]
    pub fn max_bases(&self) -> usize {
        self.max_bases
    }

    /// Returns the reads of the batch, in increasing id order.
    #[must_use]
    pub fn records(&self) -> &[ReadRecord] {
        &self.slots[..self.len]
    }

    /// Returns an iterator over the reads of the batch, in increasing id order.
    pub fn iter(&self) -> slice::Iter<'_, ReadRecord> {
        self.records().iter()
    }

    /// Returns the number of reads in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` when the batch holds no reads.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns the total number of bases in the batch.
    #[must_use]
    pub fn bases(&self) -> usize {
        self.bases
    }

    /// Returns the id of the first read, or `None` for an empty batch.
    #[must_use]
    pub fn first_id(&self) -> Option<ReadId> {
        self.records().first().map(|record| record.id)
    }

    /// Returns the heap memory the batch holds, in bytes: every record buffer's capacity
    /// (spare ones included) plus the record table itself.
    ///
    /// After [`ReadStream::next_batch`](crate::ReadStream::next_batch) each buffer holds at
    /// most four times its content or 1 KiB, whichever is larger, so this is at most four times
    /// the bytes of the names and [`bases`](Self::bases), plus about 2.3 KiB per read, plus a
    /// record table of at most about 230 KiB for small batches. For HiFi batches it is usually
    /// within a factor of two of `bases`. Takes time linear in the number of reads.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        let table = self.slots.capacity() * size_of::<ReadRecord>();
        let buffers: usize = self
            .slots
            .iter()
            .map(|slot| slot.name.capacity() + slot.seq.capacity())
            .sum();
        table + buffers
    }

    /// Empties the batch, keeping its buffers for reuse.
    pub fn clear(&mut self) {
        self.len = 0;
        self.bases = 0;
    }

    /// Frees the slots past the current reads, and the record table's excess capacity.
    pub(crate) fn release_spare_slots(&mut self) {
        self.slots.truncate(self.len);
        let wanted = self.len.max(MIN_RETAINED_SLOTS);
        if self.slots.capacity() > MAX_SLACK * wanted {
            self.slots.shrink_to(wanted);
        }
    }

    /// Returns `true` when the batch has reached its base or record budget.
    pub(crate) fn is_full(&self) -> bool {
        self.len >= MAX_BATCH_RECORDS || (self.len > 0 && self.bases >= self.max_bases)
    }

    /// Returns the next free slot, clearing it for reuse.
    pub(crate) fn next_slot(&mut self) -> &mut ReadRecord {
        if self.len == self.slots.len() {
            self.slots.push(ReadRecord::default());
        }
        let slot = &mut self.slots[self.len];
        slot.clear();
        slot
    }

    /// Marks the slot returned by [`Self::next_slot`] as part of the batch, shrinking its
    /// buffers if they are far larger than the read.
    pub(crate) fn commit_slot(&mut self) {
        let slot = &mut self.slots[self.len];
        cap_slack(&mut slot.name);
        cap_slack(&mut slot.seq);
        self.bases += slot.seq.len();
        self.len += 1;
    }
}

/// Shrinks `buf` when its capacity exceeds [`MAX_SLACK`] times what its content needs.
fn cap_slack(buf: &mut Vec<u8>) {
    let wanted = buf.len().max(MIN_RETAINED_CAPACITY);
    if buf.capacity() > MAX_SLACK * wanted {
        buf.shrink_to(wanted);
    }
}

impl Clone for ReadBatch {
    /// Copies the current reads; spare buffers are not copied.
    fn clone(&self) -> Self {
        Self {
            slots: self.records().to_vec(),
            len: self.len,
            bases: self.bases,
            max_bases: self.max_bases,
        }
    }

    /// Copies the current reads of `source`, reusing this batch's record buffers.
    fn clone_from(&mut self, source: &Self) {
        let records = source.records();
        self.slots.truncate(records.len());
        for (slot, record) in self.slots.iter_mut().zip(records) {
            slot.id = record.id;
            slot.name.clone_from(&record.name);
            slot.seq.clone_from(&record.seq);
        }
        let reused = self.slots.len();
        self.slots.extend_from_slice(&records[reused..]);
        self.len = source.len;
        self.bases = source.bases;
        self.max_bases = source.max_bases;
    }
}

impl Default for ReadBatch {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Debug for ReadBatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Spare slots are an implementation detail; show only the batch itself.
        f.debug_struct("ReadBatch")
            .field("first_id", &self.first_id())
            .field("len", &self.len)
            .field("bases", &self.bases)
            .field("max_bases", &self.max_bases)
            .finish_non_exhaustive()
    }
}

impl<'a> IntoIterator for &'a ReadBatch {
    type Item = &'a ReadRecord;
    type IntoIter = slice::Iter<'a, ReadRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
