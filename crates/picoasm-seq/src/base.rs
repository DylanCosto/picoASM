//! Nucleotide encoding and reverse complement.

/// Returns the 2-bit code of a base (`A` 0, `C` 1, `G` 2, `T` 3), in either case.
///
/// Any other byte, such as `N` or an IUPAC code, returns `None`.
#[must_use]
pub const fn encode(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Returns the complement of a base, keeping its case; other bytes become `N`.
#[must_use]
pub const fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

/// Returns the reverse complement of `seq`.
#[must_use]
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    reverse_complement_into(seq, &mut out);
    out
}

/// Replaces the contents of `out` with the reverse complement of `seq`.
///
/// Reusing `out` across calls avoids an allocation per read.
pub fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(seq.iter().rev().map(|&b| complement(b)));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_round_trips() {
        let seq = b"ACGTNacgtRA";
        let rc = reverse_complement(seq);
        assert_eq!(rc, b"TNacgtNACGT");
        // Ambiguity codes become N, so only ACGTN round-trips exactly.
        assert_eq!(
            reverse_complement(&reverse_complement(b"ACGTNacgt")),
            b"ACGTNacgt"
        );
    }

    #[test]
    fn encode_rejects_non_acgt() {
        assert_eq!(encode(b'g'), Some(2));
        assert_eq!(encode(b'N'), None);
        assert_eq!(encode(b'-'), None);
    }
}
