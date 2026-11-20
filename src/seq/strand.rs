//! Strand of a feature relative to its canonical form.

use std::ops::Not;

/// Whether a feature was seen as its canonical form or as its reverse complement.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    /// The feature appears in its canonical orientation.
    #[default]
    Forward,
    /// The feature appears as the reverse complement of its canonical form.
    Reverse,
}

impl Strand {
    /// Returns `Forward` for `Forward ∘ Forward` and `Reverse ∘ Reverse`, else `Reverse`.
    ///
    /// Composes two orientations, for example a node's strand within a read and the read's
    /// strand within a unitig.
    #[must_use]
    pub const fn compose(self, other: Self) -> Self {
        match (self, other) {
            (Self::Forward, Self::Forward) | (Self::Reverse, Self::Reverse) => Self::Forward,
            _ => Self::Reverse,
        }
    }

    /// Returns `true` for [`Strand::Reverse`].
    #[must_use]
    pub const fn is_reverse(self) -> bool {
        matches!(self, Self::Reverse)
    }
}

impl Not for Strand {
    type Output = Self;

    fn not(self) -> Self {
        match self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compose_matches_xor() {
        for a in [Strand::Forward, Strand::Reverse] {
            for b in [Strand::Forward, Strand::Reverse] {
                assert_eq!(a.compose(b).is_reverse(), a.is_reverse() ^ b.is_reverse());
            }
            assert_eq!(!!a, a);
        }
    }
}
