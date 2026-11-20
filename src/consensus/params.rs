//! Scoring and banding parameters for partial-order alignment.

use thiserror::Error;

/// Scoring and banding parameters for [`PoaBackend`](crate::consensus::PoaBackend).
///
/// Scores are small unsigned integers so that alignment scores can never overflow. A gap of
/// length `g` costs `gap_open + g * gap_extend`; set `gap_open` to zero for linear gaps.
///
/// The defaults (match 2, mismatch 4, gap open 4, gap extend 2) are abPOA's first affine gap
/// function, which suits HiFi reads. Validate custom values with [`PoaParams::validate`]; the
/// backend constructor does so.
///
/// The struct is `#[non_exhaustive]` so that parameters can be added later; outside this crate,
/// start from [`PoaParams::new`] or [`Default::default`] and adjust with the `with_*` methods.
///
/// # Examples
///
/// ```
/// use picoasm::consensus::PoaParams;
///
/// let params = PoaParams::new().with_band_width(64).with_gap_open(0);
/// assert_eq!(params.band_width, 64);
/// assert_eq!(params.validate(), Ok(()));
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct PoaParams {
    /// Score added for each matching base pair; must be positive.
    pub match_score: u8,
    /// Penalty subtracted for each mismatching base pair; must be positive.
    pub mismatch_penalty: u8,
    /// Penalty subtracted once for each gap; may be zero.
    pub gap_open: u8,
    /// Penalty subtracted for each base of a gap; must be positive.
    pub gap_extend: u8,
    /// Initial half-width of the adaptive band, in segment bases; must be positive.
    ///
    /// Each graph node is aligned against segment positions within this distance of the best
    /// cells of its predecessors. When the optimal path touches the band edge, the band is
    /// doubled and the segment realigned, so this is a speed knob rather than a hard limit.
    pub band_width: usize,
}

impl PoaParams {
    /// Default score for a matching base pair.
    pub const DEFAULT_MATCH_SCORE: u8 = 2;
    /// Default penalty for a mismatching base pair.
    pub const DEFAULT_MISMATCH_PENALTY: u8 = 4;
    /// Default gap opening penalty.
    pub const DEFAULT_GAP_OPEN: u8 = 4;
    /// Default gap extension penalty.
    pub const DEFAULT_GAP_EXTEND: u8 = 2;
    /// Default initial band half-width, in bases.
    ///
    /// HiFi windows rarely drift more than a few bases from the diagonal, and the band follows
    /// the best cell of each row, so 32 leaves a wide margin.
    pub const DEFAULT_BAND_WIDTH: usize = 32;

    /// Returns the default parameters, the same as [`Default::default`].
    #[must_use]
    pub const fn new() -> Self {
        Self {
            match_score: Self::DEFAULT_MATCH_SCORE,
            mismatch_penalty: Self::DEFAULT_MISMATCH_PENALTY,
            gap_open: Self::DEFAULT_GAP_OPEN,
            gap_extend: Self::DEFAULT_GAP_EXTEND,
            band_width: Self::DEFAULT_BAND_WIDTH,
        }
    }

    /// Returns these parameters with the match score replaced.
    #[must_use]
    pub const fn with_match_score(mut self, match_score: u8) -> Self {
        self.match_score = match_score;
        self
    }

    /// Returns these parameters with the mismatch penalty replaced.
    #[must_use]
    pub const fn with_mismatch_penalty(mut self, mismatch_penalty: u8) -> Self {
        self.mismatch_penalty = mismatch_penalty;
        self
    }

    /// Returns these parameters with the gap opening penalty replaced.
    #[must_use]
    pub const fn with_gap_open(mut self, gap_open: u8) -> Self {
        self.gap_open = gap_open;
        self
    }

    /// Returns these parameters with the gap extension penalty replaced.
    #[must_use]
    pub const fn with_gap_extend(mut self, gap_extend: u8) -> Self {
        self.gap_extend = gap_extend;
        self
    }

    /// Returns these parameters with the initial band half-width replaced.
    #[must_use]
    pub const fn with_band_width(mut self, band_width: usize) -> Self {
        self.band_width = band_width;
        self
    }

    /// Checks that the parameters are usable.
    ///
    /// # Errors
    ///
    /// Returns an error if the match score, mismatch penalty, gap extension penalty or band
    /// width is zero.
    pub const fn validate(&self) -> Result<(), InvalidPoaParamsError> {
        if self.match_score == 0 {
            return Err(InvalidPoaParamsError::MatchScore);
        }
        if self.mismatch_penalty == 0 {
            return Err(InvalidPoaParamsError::MismatchPenalty);
        }
        if self.gap_extend == 0 {
            return Err(InvalidPoaParamsError::GapExtend);
        }
        if self.band_width == 0 {
            return Err(InvalidPoaParamsError::BandWidth);
        }
        Ok(())
    }
}

impl Default for PoaParams {
    fn default() -> Self {
        Self::new()
    }
}

/// The error returned by [`PoaParams::validate`] for unusable parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Error)]
#[non_exhaustive]
pub enum InvalidPoaParamsError {
    /// The match score is zero.
    #[error("match score must be positive")]
    MatchScore,
    /// The mismatch penalty is zero.
    #[error("mismatch penalty must be positive")]
    MismatchPenalty,
    /// The gap extension penalty is zero.
    #[error("gap extension penalty must be positive")]
    GapExtend,
    /// The band width is zero.
    #[error("band width must be positive")]
    BandWidth,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_valid() {
        assert_eq!(PoaParams::new(), PoaParams::default());
        assert_eq!(PoaParams::default().validate(), Ok(()));
        let linear = PoaParams {
            gap_open: 0,
            ..PoaParams::default()
        };
        assert_eq!(linear.validate(), Ok(()));
    }

    #[test]
    fn setters_replace_one_field() {
        let params = PoaParams::new()
            .with_match_score(1)
            .with_mismatch_penalty(3)
            .with_gap_open(5)
            .with_gap_extend(6)
            .with_band_width(7);
        let want = PoaParams {
            match_score: 1,
            mismatch_penalty: 3,
            gap_open: 5,
            gap_extend: 6,
            band_width: 7,
        };
        assert_eq!(params, want);
    }

    #[test]
    fn rejects_zero_values() {
        let base = PoaParams::default();
        let cases = [
            (
                PoaParams {
                    match_score: 0,
                    ..base
                },
                InvalidPoaParamsError::MatchScore,
            ),
            (
                PoaParams {
                    mismatch_penalty: 0,
                    ..base
                },
                InvalidPoaParamsError::MismatchPenalty,
            ),
            (
                PoaParams {
                    gap_extend: 0,
                    ..base
                },
                InvalidPoaParamsError::GapExtend,
            ),
            (
                PoaParams {
                    band_width: 0,
                    ..base
                },
                InvalidPoaParamsError::BandWidth,
            ),
        ];
        for (params, want) in cases {
            assert_eq!(params.validate(), Err(want), "params {params:?}");
            assert!(!want.to_string().is_empty());
        }
    }
}
