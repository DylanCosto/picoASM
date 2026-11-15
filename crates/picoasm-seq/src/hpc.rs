//! Homopolymer compression.

use crate::{ReadTooLongError, checked_read_len};

/// A homopolymer-compressed read: each run of one byte becomes a single byte.
///
/// HiFi errors are mostly homopolymer-length errors, so features chosen in compressed space are
/// far more stable between reads. The compressed sequence keeps, for every compressed position,
/// where its run starts in the raw read, so features map back to raw coordinates.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Hpc {
    bases: Vec<u8>,
    /// `run_starts[i]` is the raw start of run `i`; the final entry is the raw length.
    run_starts: Vec<u32>,
}

impl Hpc {
    /// Returns an empty compressed sequence.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Compresses `raw` into `self`, reusing its buffers.
    ///
    /// Runs are compared case-insensitively, and bases are stored uppercase. Bytes that are not
    /// `A`, `C`, `G` or `T` are kept as `N`, so they still break k-mers downstream.
    ///
    /// # Errors
    ///
    /// Returns an error if `raw` is longer than `u32::MAX` bases; `self` is then left empty.
    pub fn compress(&mut self, raw: &[u8]) -> Result<(), ReadTooLongError> {
        self.bases.clear();
        self.run_starts.clear();
        let len = checked_read_len(raw.len())?;

        let mut previous = None;
        for (pos, &byte) in (0..len).zip(raw) {
            let base = normalize(byte);
            if previous != Some(base) {
                self.bases.push(base);
                self.run_starts.push(pos);
                previous = Some(base);
            }
        }
        self.run_starts.push(len);
        Ok(())
    }

    /// Returns the compressed bases.
    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// Returns the number of compressed bases.
    #[must_use]
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// Returns `true` if there are no compressed bases.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// Returns the raw position where compressed position `i` starts.
    ///
    /// `i == self.len()` returns the raw length, so `raw_start(i + n)` is the raw end of the
    /// compressed range `i..i + n`.
    ///
    /// # Panics
    ///
    /// Panics if `i > self.len()`.
    #[must_use]
    pub fn raw_start(&self, i: usize) -> u32 {
        self.run_starts[i]
    }
}

const fn normalize(byte: u8) -> u8 {
    match byte {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compresses_runs_and_maps_back() {
        let mut hpc = Hpc::new();
        hpc.compress(b"AAcGGGtTNNa").unwrap();
        assert_eq!(hpc.bases(), b"ACGTNA");
        let starts: Vec<u32> = (0..=hpc.len()).map(|i| hpc.raw_start(i)).collect();
        assert_eq!(starts, [0, 2, 3, 6, 8, 10, 11]);
    }

    #[test]
    fn empty_read() {
        let mut hpc = Hpc::new();
        hpc.compress(b"").unwrap();
        assert!(hpc.is_empty());
        assert_eq!(hpc.raw_start(0), 0);
    }

    #[test]
    fn reuse_replaces_previous_read() {
        let mut hpc = Hpc::new();
        hpc.compress(b"GGGGTTTTAC").unwrap();
        hpc.compress(b"CCA").unwrap();
        assert_eq!(hpc.bases(), b"CA");
        assert_eq!(hpc.raw_start(hpc.len()), 3);
    }
}
