//! Deterministic 64-bit hashing.

/// Mixes a 64-bit value with an invertible function (the `splitmix64` finalizer).
///
/// Invertibility means distinct inputs never collide, so a k-mer of up to 31 bases (62 bits) maps
/// to a unique hash.
#[must_use]
pub const fn mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Hashes a sequence of hashes, in order.
///
/// The result depends on the order of `hashes`, so a sequence and its reverse usually differ.
#[must_use]
pub fn combine(hashes: impl IntoIterator<Item = u64>) -> u64 {
    // An arbitrary odd seed keeps the empty sequence and a single zero apart.
    const SEED: u64 = 0x9e37_79b9_7f4a_7c15;
    hashes
        .into_iter()
        .fold(SEED, |acc, h| mix64(acc.rotate_left(23) ^ h))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combine_is_order_sensitive() {
        assert_ne!(combine([1, 2, 3]), combine([3, 2, 1]));
        assert_eq!(combine([1, 2, 3]), combine(vec![1, 2, 3]));
        assert_ne!(combine([]), combine([0]));
    }

    #[test]
    fn mix64_is_injective_on_small_inputs() {
        let mut seen: Vec<u64> = (0..10_000).map(mix64).collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 10_000);
    }
}
