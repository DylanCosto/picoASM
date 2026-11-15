//! Anchors: universe minimizers chosen in homopolymer-compressed space.

use thiserror::Error;

use crate::{ReadTooLongError, Strand, base, checked_read_len, mix64};

/// The longest anchor k-mer: 31 bases fit in 62 bits of a `u64`.
pub const MAX_KMER_LEN: u8 = 31;

/// Parameters that decide which k-mers become anchors and how anchors form k-min-mers.
///
/// Defaults follow rust-mdbg's settings for human HiFi reads and are not yet tuned for picoasm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SketchParams {
    kmer_len: u8,
    density: f64,
    window: u16,
}

impl SketchParams {
    /// Default anchor k-mer length, in homopolymer-compressed bases.
    pub const DEFAULT_KMER_LEN: u8 = 14;
    /// Default fraction of compressed k-mers kept as anchors.
    pub const DEFAULT_DENSITY: f64 = 0.003;
    /// Default number of consecutive anchors in a k-min-mer.
    pub const DEFAULT_WINDOW: u16 = 21;

    /// Validates and returns sketch parameters.
    ///
    /// `kmer_len` is the anchor length in compressed bases, `density` the expected fraction of
    /// compressed k-mers kept as anchors, and `window` the number of consecutive anchors in a
    /// k-min-mer. A window of 1 makes each anchor its own k-min-mer; graph construction needs at
    /// least 2 so that consecutive k-min-mers overlap.
    ///
    /// # Errors
    ///
    /// Returns an error if `kmer_len` is not in `1..=31`, `density` is not in `(0, 1]`, or
    /// `window` is zero.
    pub fn new(kmer_len: u8, density: f64, window: u16) -> Result<Self, InvalidSketchParamsError> {
        if !(1..=MAX_KMER_LEN).contains(&kmer_len) {
            return Err(InvalidSketchParamsError::KmerLen(kmer_len));
        }
        if !(density > 0.0 && density <= 1.0) {
            return Err(InvalidSketchParamsError::Density(density));
        }
        if window == 0 {
            return Err(InvalidSketchParamsError::Window);
        }
        Ok(Self {
            kmer_len,
            density,
            window,
        })
    }

    /// Returns the anchor k-mer length, in compressed bases.
    #[must_use]
    pub const fn kmer_len(&self) -> u8 {
        self.kmer_len
    }

    /// Returns the expected fraction of compressed k-mers kept as anchors.
    #[must_use]
    pub const fn density(&self) -> f64 {
        self.density
    }

    /// Returns the number of consecutive anchors in a k-min-mer.
    #[must_use]
    pub const fn window(&self) -> usize {
        self.window as usize
    }

    /// Returns the exclusive hash threshold `density * 2^64`, saturated to `u64::MAX`.
    ///
    /// A k-mer is an anchor when its canonical hash is below the threshold, except at density
    /// 1.0, where every hash is kept (including `u64::MAX`, which the saturated threshold cannot
    /// express). [`SketchParams::keeps_hash`] applies both rules.
    #[must_use]
    pub fn hash_threshold(&self) -> u64 {
        if self.density >= 1.0 {
            u64::MAX
        } else {
            // 2^64 as an exact f64.
            const TWO_POW_64: f64 = 18_446_744_073_709_551_616.0;
            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "density is in (0, 1), so the product is in [0, 2^64) and the cast saturates"
            )]
            let threshold = (self.density * TWO_POW_64) as u64;
            threshold
        }
    }

    /// Returns `true` if a k-mer whose canonical hash is `hash` is kept as an anchor.
    #[must_use]
    pub fn keeps_hash(&self, hash: u64) -> bool {
        HashFilter::new(self).keeps(hash)
    }
}

// The density is validated to lie in (0, 1], so it is never NaN and equality is reflexive.
impl Eq for SketchParams {}

impl Default for SketchParams {
    fn default() -> Self {
        Self {
            kmer_len: Self::DEFAULT_KMER_LEN,
            density: Self::DEFAULT_DENSITY,
            window: Self::DEFAULT_WINDOW,
        }
    }
}

/// The error returned by [`SketchParams::new`] for out-of-range parameters.
#[derive(Clone, Copy, Debug, PartialEq, Error)]
#[non_exhaustive]
pub enum InvalidSketchParamsError {
    /// The anchor k-mer length is not in `1..=31`.
    #[error("anchor k-mer length must be between 1 and {MAX_KMER_LEN}, got {0}")]
    KmerLen(u8),
    /// The density is not in `(0, 1]`.
    #[error("anchor density must be in (0, 1], got {0}")]
    Density(f64),
    /// The k-min-mer window is zero.
    #[error("k-min-mer window must be at least 1")]
    Window,
}

/// A k-mer selected as a sparse, strand-independent feature of a read.
///
/// The derived order compares `hash`, then `strand`, then the raw span.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Anchor {
    /// Canonical hash: the hash of the smaller 2-bit encoding of the k-mer and its reverse
    /// complement.
    pub hash: u64,
    /// [`Strand::Forward`] if the k-mer as read is the canonical (smaller) encoding.
    pub strand: Strand,
    /// Raw read position of the k-mer's first base.
    pub raw_start: u32,
    /// Raw read position one past the end of the k-mer's last homopolymer run.
    pub raw_end: u32,
}

/// Picks anchors from reads, reusing its buffers between reads.
#[derive(Clone, Debug)]
pub struct Sketcher {
    params: SketchParams,
    filter: HashFilter,
}

/// Which canonical k-mer hashes become anchors, precomputed from the density.
#[derive(Clone, Copy, Debug)]
enum HashFilter {
    /// Density 1.0: every hash.
    All,
    /// Hashes strictly below the threshold.
    Below(u64),
}

impl HashFilter {
    fn new(params: &SketchParams) -> Self {
        if params.density() >= 1.0 {
            Self::All
        } else {
            Self::Below(params.hash_threshold())
        }
    }

    const fn keeps(self, hash: u64) -> bool {
        match self {
            Self::All => true,
            Self::Below(threshold) => hash < threshold,
        }
    }
}

impl Sketcher {
    /// Returns a sketcher for `params`.
    #[must_use]
    pub fn new(params: SketchParams) -> Self {
        Self {
            params,
            filter: HashFilter::new(&params),
        }
    }

    /// Returns the parameters this sketcher uses.
    #[must_use]
    pub const fn params(&self) -> SketchParams {
        self.params
    }

    /// Replaces the contents of `anchors` with the anchors of `read`, in increasing
    /// [`Anchor::raw_start`] order.
    ///
    /// The k-mers are those of the read's homopolymer-compressed form (as [`crate::Hpc`] gives it),
    /// compressed on the fly: each anchor's raw end, where the run after its last base starts,
    /// is filled in when that run starts.
    ///
    /// # Errors
    ///
    /// Returns an error if `read` is longer than `u32::MAX` bases.
    pub fn sketch(
        &mut self,
        read: &[u8],
        anchors: &mut Vec<Anchor>,
    ) -> Result<(), ReadTooLongError> {
        /// Raw starts of the last compressed bases kept, by compressed index modulo this: a
        /// k-mer spans at most [`MAX_KMER_LEN`] of them.
        const RING: usize = 32;
        const _: () = assert!(
            RING > MAX_KMER_LEN as usize,
            "the ring holds a k-mer's runs"
        );
        anchors.clear();
        let len = checked_read_len(read.len())?;

        let k = usize::from(self.params.kmer_len);
        let mask = (1_u64 << (2 * k)) - 1;
        let rc_shift = 2 * (k - 1);
        let mut forward = 0_u64;
        let mut reverse = 0_u64;
        let mut valid = 0_usize;
        let mut starts = [0_u32; RING];
        // Compressed bases so far, and the code of the last.
        let (mut compressed, mut previous) = (0_usize, NO_CODE);
        // The anchor kept last, whose raw end is where the next run starts.
        let mut pending: Option<(u64, Strand, u32)> = None;

        for (pos, &byte) in (0..len).zip(read) {
            let code = CODES[usize::from(byte)];
            if code == previous {
                continue;
            }
            previous = code;
            if let Some((hash, strand, raw_start)) = pending.take() {
                anchors.push(Anchor {
                    hash,
                    strand,
                    raw_start,
                    raw_end: pos,
                });
            }
            let end = compressed;
            starts[end % RING] = pos;
            compressed += 1;
            if code == N_CODE {
                valid = 0;
                continue;
            }
            let code = u64::from(code);
            forward = ((forward << 2) | code) & mask;
            reverse = (reverse >> 2) | ((3 - code) << rc_shift);
            valid += 1;
            if valid < k {
                continue;
            }

            // Hash the canonical encoding rather than taking the smaller of two hashes: the
            // minimum of two uniform hashes would double the anchor density. A reverse-
            // complement palindrome has no orientation to record. (The smaller encoding is
            // taken without a branch: which one it is, is a coin toss each base.)
            if forward == reverse {
                continue;
            }
            let hash = mix64(forward.min(reverse));
            if !self.filter.keeps(hash) {
                continue;
            }
            let strand = if forward < reverse {
                Strand::Forward
            } else {
                Strand::Reverse
            };
            pending = Some((hash, strand, starts[(end + 1 - k) % RING]));
        }
        if let Some((hash, strand, raw_start)) = pending {
            anchors.push(Anchor {
                hash,
                strand,
                raw_start,
                raw_end: len,
            });
        }
        Ok(())
    }
}

/// The code [`CODES`] gives anything but A, C, G and T: as [`crate::Hpc`] does, all of it is
/// `N`.
const N_CODE: u8 = 4;

/// A code no byte has, before a read's first.
const NO_CODE: u8 = u8::MAX;

/// Each byte's 2-bit code (A 0, C 1, G 2, T 3, either case), or [`N_CODE`].
const CODES: [u8; 256] = {
    let mut codes = [N_CODE; 256];
    let mut byte = 0;
    while byte < 256 {
        #[expect(clippy::cast_possible_truncation, reason = "below 256")]
        if let Some(code) = base::encode(byte as u8) {
            codes[byte] = code;
        }
        byte += 1;
    }
    codes
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Hpc, reverse_complement};

    /// The anchors of `read` found as they were before compressing on the fly: the whole read
    /// compressed first ([`Hpc`]), then its k-mers.
    fn anchors_of_compressed(read: &[u8], params: SketchParams) -> Vec<Anchor> {
        let mut hpc = Hpc::new();
        hpc.compress(read).unwrap();
        let filter = HashFilter::new(&params);
        let k = usize::from(params.kmer_len);
        let mask = (1_u64 << (2 * k)) - 1;
        let rc_shift = 2 * (k - 1);
        let (mut forward, mut reverse, mut valid) = (0_u64, 0_u64, 0_usize);
        let mut anchors = Vec::new();
        for (end, &b) in hpc.bases().iter().enumerate() {
            let Some(code) = base::encode(b) else {
                valid = 0;
                continue;
            };
            let code = u64::from(code);
            forward = ((forward << 2) | code) & mask;
            reverse = (reverse >> 2) | ((3 - code) << rc_shift);
            valid += 1;
            if valid < k {
                continue;
            }
            let (canonical, strand) = match forward.cmp(&reverse) {
                std::cmp::Ordering::Less => (forward, Strand::Forward),
                std::cmp::Ordering::Greater => (reverse, Strand::Reverse),
                std::cmp::Ordering::Equal => continue,
            };
            let hash = mix64(canonical);
            if filter.keeps(hash) {
                anchors.push(Anchor {
                    hash,
                    strand,
                    raw_start: hpc.raw_start(end + 1 - k),
                    raw_end: hpc.raw_start(end + 1),
                });
            }
        }
        anchors
    }

    #[test]
    fn anchors_are_those_of_the_compressed_read() {
        // Runs, lowercase, N and other bytes, at every k and several densities.
        let mut state = 9_u64;
        let mut next = move |n: u64| {
            state = mix64(state.wrapping_add(0x9e37_79b9_7f4a_7c15));
            state % n
        };
        for case in 0..300 {
            let mut read = Vec::new();
            for _ in 0..next(3_000) {
                let byte = match next(100) {
                    0 => b'N',
                    1 => b'R',
                    2..=9 => b"acgt"[usize::try_from(next(4)).unwrap()],
                    _ => b"ACGT"[usize::try_from(next(4)).unwrap()],
                };
                for _ in 0..1 + next(3) {
                    read.push(byte);
                }
            }
            let k = u8::try_from(1 + next(u64::from(MAX_KMER_LEN))).unwrap();
            let density = [1.0, 0.5, 0.05][usize::try_from(next(3)).unwrap()];
            let params = SketchParams::new(k, density, 1).unwrap();
            let mut anchors = Vec::new();
            Sketcher::new(params).sketch(&read, &mut anchors).unwrap();
            assert_eq!(
                anchors,
                anchors_of_compressed(&read, params),
                "case {case}, k {k}"
            );
        }
    }

    /// A deterministic pseudo-random sequence for tests.
    fn random_seq(len: usize, seed: u64) -> Vec<u8> {
        (0..len as u64)
            .map(|i| b"ACGT"[(mix64(seed ^ i) % 4) as usize])
            .collect()
    }

    #[test]
    fn rejects_bad_params() {
        assert!(SketchParams::new(0, 0.1, 3).is_err());
        assert!(SketchParams::new(32, 0.1, 3).is_err());
        assert!(SketchParams::new(15, 0.0, 3).is_err());
        assert!(SketchParams::new(15, 1.5, 3).is_err());
        assert!(SketchParams::new(15, f64::NAN, 3).is_err());
        assert!(SketchParams::new(15, 0.1, 0).is_err());
        assert!(SketchParams::new(15, 1.0, 1).is_ok());
    }

    #[test]
    fn density_is_roughly_respected() {
        let read = random_seq(20_000, 1);
        let mut sketcher = Sketcher::new(SketchParams::new(15, 0.05, 3).unwrap());
        let mut anchors = Vec::new();
        sketcher.sketch(&read, &mut anchors).unwrap();
        let mut hpc = Hpc::new();
        hpc.compress(&read).unwrap();
        let expected = hpc.len() as f64 * 0.05;
        let got = anchors.len() as f64;
        assert!(
            (got - expected).abs() < expected * 0.3,
            "got {got}, expected about {expected}"
        );
        assert!(anchors.windows(2).all(|w| w[0].raw_start < w[1].raw_start));
    }

    #[test]
    fn anchors_are_strand_symmetric() {
        let read = random_seq(5_000, 7);
        let rc = reverse_complement(&read);
        let len = u32::try_from(read.len()).unwrap();
        let mut sketcher = Sketcher::new(SketchParams::new(13, 0.1, 3).unwrap());

        let mut forward = Vec::new();
        sketcher.sketch(&read, &mut forward).unwrap();
        let mut reverse = Vec::new();
        sketcher.sketch(&rc, &mut reverse).unwrap();

        let mirrored: Vec<Anchor> = reverse
            .iter()
            .rev()
            .map(|a| Anchor {
                hash: a.hash,
                strand: !a.strand,
                raw_start: len - a.raw_end,
                raw_end: len - a.raw_start,
            })
            .collect();
        assert!(!forward.is_empty());
        assert_eq!(forward, mirrored);
    }

    #[test]
    fn homopolymer_errors_do_not_change_anchors() {
        let read = random_seq(3_000, 3);
        // Lengthen every homopolymer run by one base, as a HiFi indel error would.
        let mut stretched = Vec::new();
        for (i, &b) in read.iter().enumerate() {
            stretched.push(b);
            if read.get(i + 1) != Some(&b) {
                stretched.push(b);
            }
        }
        let mut sketcher = Sketcher::new(SketchParams::new(13, 0.1, 3).unwrap());
        let mut a = Vec::new();
        sketcher.sketch(&read, &mut a).unwrap();
        let mut b = Vec::new();
        sketcher.sketch(&stretched, &mut b).unwrap();
        let hashes = |v: &[Anchor]| v.iter().map(|x| x.hash).collect::<Vec<_>>();
        assert_eq!(hashes(&a), hashes(&b));
    }

    /// Sketches `read` with `kmer_len` and density 1.0, returning `(raw_start, raw_end, strand)`.
    fn spans(read: &[u8], kmer_len: u8) -> Vec<(u32, u32, Strand)> {
        let mut sketcher = Sketcher::new(SketchParams::new(kmer_len, 1.0, 2).unwrap());
        let mut anchors = Vec::new();
        sketcher.sketch(read, &mut anchors).unwrap();
        anchors
            .iter()
            .map(|a| (a.raw_start, a.raw_end, a.strand))
            .collect()
    }

    /// Picks every non-palindromic k-mer of an N-free, run-free read by direct encoding.
    fn reference_anchors(read: &[u8], k: usize) -> Vec<Anchor> {
        let encode_kmer = |kmer: &[u8]| {
            kmer.iter().fold(0_u64, |acc, &b| {
                (acc << 2) | u64::from(base::encode(b).unwrap())
            })
        };
        read.windows(k)
            .zip(0_u32..)
            .filter_map(|(kmer, start)| {
                let forward = encode_kmer(kmer);
                let reverse = encode_kmer(&reverse_complement(kmer));
                let (canonical, strand) = match forward.cmp(&reverse) {
                    std::cmp::Ordering::Less => (forward, Strand::Forward),
                    std::cmp::Ordering::Greater => (reverse, Strand::Reverse),
                    std::cmp::Ordering::Equal => return None,
                };
                Some(Anchor {
                    hash: mix64(canonical),
                    strand,
                    raw_start: start,
                    raw_end: start + u32::try_from(k).unwrap(),
                })
            })
            .collect()
    }

    #[test]
    fn non_acgt_breaks_kmers() {
        assert_eq!(
            spans(b"ACGANCGTA", 4),
            [(0, 4, Strand::Forward), (5, 9, Strand::Forward)]
        );
    }

    #[test]
    fn raw_end_covers_the_last_run() {
        // HPC "ACTG" with runs AA, CCC, TT, G.
        assert_eq!(
            spans(b"AACCCTTG", 2),
            [
                (0, 5, Strand::Forward),
                (2, 7, Strand::Reverse),
                (5, 8, Strand::Reverse)
            ]
        );
    }

    #[test]
    fn matches_reference_at_k_boundaries() {
        // Removing repeated bases keeps HPC positions equal to raw positions.
        let mut read = random_seq(400, 11);
        read.dedup();
        for k in [1_u8, 2, 30, 31] {
            let mut sketcher = Sketcher::new(SketchParams::new(k, 1.0, 2).unwrap());
            let mut anchors = Vec::new();
            sketcher.sketch(&read, &mut anchors).unwrap();
            assert_eq!(anchors, reference_anchors(&read, usize::from(k)), "k = {k}");
        }
    }

    #[test]
    fn skips_even_k_palindromes() {
        assert!(spans(b"ACGT", 4).is_empty());
        assert!(spans(b"AT", 2).is_empty());
        assert_eq!(spans(b"ACGTA", 4), [(1, 5, Strand::Forward)]);
    }

    #[test]
    fn lowercase_matches_uppercase() {
        let read = random_seq(2_000, 5);
        let lower = read.to_ascii_lowercase();
        let mut sketcher = Sketcher::new(SketchParams::new(9, 0.2, 3).unwrap());
        let mut a = Vec::new();
        sketcher.sketch(&read, &mut a).unwrap();
        let mut b = Vec::new();
        sketcher.sketch(&lower, &mut b).unwrap();
        assert!(!a.is_empty());
        assert_eq!(a, b);
    }

    #[test]
    fn short_or_n_only_reads_have_no_anchors() {
        assert!(spans(b"ACG", 4).is_empty());
        assert!(spans(b"", 1).is_empty());
        assert!(spans(b"NNNNNNNNNN", 1).is_empty());
        assert!(spans(b"ACNGTNCA", 3).is_empty());
    }

    #[test]
    fn full_density_keeps_every_hash() {
        let params = SketchParams::new(15, 1.0, 2).unwrap();
        assert!(params.keeps_hash(u64::MAX));
        assert!(params.keeps_hash(0));
        let half = SketchParams::new(15, 0.5, 2).unwrap();
        assert!(half.keeps_hash((1 << 63) - 1));
        assert!(!half.keeps_hash(1 << 63));
    }

    #[test]
    fn checked_read_len_bounds() {
        assert_eq!(checked_read_len(0), Ok(0));
        assert_eq!(checked_read_len(u32::MAX as usize), Ok(u32::MAX));
        if let Some(len) = (u32::MAX as usize).checked_add(1) {
            assert_eq!(checked_read_len(len).map_err(|e| e.read_len()), Err(len));
        }
    }
}
