//! Read lengths in `u32` coordinates.

use thiserror::Error;

/// The error returned when a read is too long for `u32` coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
#[error("read of {len} bases is longer than the supported maximum of {max}", max = u32::MAX)]
pub struct ReadTooLongError {
    len: usize,
}

impl ReadTooLongError {
    /// Returns the length of the rejected read.
    #[must_use]
    pub const fn read_len(&self) -> usize {
        self.len
    }
}

/// Returns a read length of `len` bases as a `u32`, the coordinate type of every picoasm crate.
///
/// Readers and sketchers call this once per read, so later `u32` arithmetic on positions within
/// the read cannot overflow.
///
/// # Errors
///
/// Returns [`ReadTooLongError`] if `len` is greater than `u32::MAX`.
pub fn checked_read_len(len: usize) -> Result<u32, ReadTooLongError> {
    u32::try_from(len).map_err(|_| ReadTooLongError { len })
}
