//! Rounds of correction, each followed by finding the overlaps again, as hifiasm runs them.

use std::time::Instant;

use crate::overlap::correct::{CorrectionParams, correct};
use crate::overlap::overlap::{
    InvalidOverlapParamsError, Overlap, OverlapParams, Pair, find_overlaps_skipping,
};
use crate::overlap::reads::Reads;

/// What one round of [`correct_in_rounds`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RoundStats {
    /// Sites found, over all reads.
    pub sites: u64,
    /// Bases changed, inserted or deleted, over all reads.
    pub changes: u64,
    /// Reads whose bases changed.
    pub changed_reads: u64,
    /// The overlaps that join reads that differ at a site.
    pub disagreeing: u64,
    /// The overlaps found after the round.
    pub overlaps: u64,
}

/// The reads after the rounds of correction, and the overlaps found on them.
#[derive(Clone, Debug)]
pub struct Corrected {
    /// The corrected reads, in the input order.
    pub reads: Reads,
    /// The overlaps between the corrected reads.
    pub overlaps: Vec<Overlap>,
    /// For each overlap, whether neither read has the other's base at one of its sites (the
    /// two come from one copy).
    pub consistent: Vec<bool>,
    /// The overlaps found before the first round.
    pub first_overlaps: u64,
    /// One entry per round.
    pub rounds: Vec<RoundStats>,
}

/// Finds the overlaps between `reads`, then `rounds` times corrects the reads from them and
/// finds them again.
///
/// After a round, the pairs none of whose chains could be aligned to the reads' ends in an
/// earlier round are not aligned again (see [`find_overlaps_skipping`]); the last round's
/// overlaps are freed before the next are found. The result is the same for any number of
/// threads; work runs on the caller's rayon pool.
///
/// # Errors
///
/// Returns an error if `overlap_params` cannot be used.
///
/// # Panics
///
/// Panics if a read is longer than 2^31 bases.
pub fn correct_in_rounds(
    reads: Reads,
    overlap_params: &OverlapParams,
    correction_params: &CorrectionParams,
    rounds: u32,
) -> Result<Corrected, InvalidOverlapParamsError> {
    let mut reads = reads;
    let started = Instant::now();
    let mut found = find_overlaps_skipping(&reads, overlap_params, &[])?;
    log::debug!(
        "overlap: first search, {:.1} s",
        started.elapsed().as_secs_f64()
    );
    let first_overlaps = found.overlaps.len() as u64;
    // Pairs unaligned in some round so far.
    let mut unaligned: Vec<Pair> = std::mem::take(&mut found.unaligned);
    let mut stats = Vec::new();
    for _ in 0..rounds {
        let started = Instant::now();
        let c = correct(&reads, &found.overlaps, correction_params);
        log::debug!(
            "overlap: correction, {:.1} s",
            started.elapsed().as_secs_f64()
        );
        reads = c.reads;
        let disagreeing = c.consistent.iter().filter(|&&x| !x).count() as u64;
        let changed_reads = c.changed.iter().filter(|&&x| x).count() as u64;
        drop(std::mem::take(&mut found.overlaps));
        let started = Instant::now();
        found = find_overlaps_skipping(&reads, overlap_params, &unaligned)?;
        log::debug!("overlap: search, {:.1} s", started.elapsed().as_secs_f64());
        unaligned.append(&mut found.unaligned);
        unaligned.sort_unstable();
        unaligned.dedup();
        stats.push(RoundStats {
            sites: c.sites,
            changes: c.changes,
            changed_reads,
            disagreeing,
            overlaps: found.overlaps.len() as u64,
        });
    }
    // The sites of the corrected reads, for the overlaps between them; the bases this would
    // change are left as they are.
    let consistent = correct(&reads, &found.overlaps, correction_params).consistent;
    Ok(Corrected {
        reads,
        overlaps: found.overlaps,
        consistent,
        first_overlaps,
        rounds: stats,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A random genome from a fixed linear congruential sequence.
    fn genome(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
            })
            .collect()
    }

    #[test]
    fn rounds_correct_an_error_and_keep_the_overlaps() {
        let g = genome(8_000, 3);
        let mut tiled: Vec<Vec<u8>> = (0..6)
            .map(|i| g[i * 500..i * 500 + 5_000].to_vec())
            .collect();
        // One read with a substitution the others outvote.
        tiled[2][2_000] = if tiled[2][2_000] == b'A' { b'C' } else { b'A' };
        let reads: Reads = tiled.iter().map(Vec::as_slice).collect();
        let params = OverlapParams {
            density: 0.2,
            min_overlap: 500,
            ..OverlapParams::default()
        };
        let corrected =
            correct_in_rounds(reads, &params, &CorrectionParams::default(), 2).expect("valid");
        assert_eq!(corrected.rounds.len(), 2);
        assert_eq!(corrected.rounds[0].changes, 1);
        assert_eq!(corrected.rounds[1].changes, 0);
        for (i, read) in corrected.reads.indices().enumerate() {
            assert_eq!(corrected.reads.get(read), &g[i * 500..i * 500 + 5_000]);
        }
        // Every pair of the six reads overlaps by 2,500 bases or more.
        assert_eq!(corrected.first_overlaps, 15);
        assert_eq!(corrected.overlaps.len(), 15);
        assert!(corrected.overlaps.iter().all(|o| o.differences() == 0));
        assert!(corrected.consistent.iter().all(|&c| c));
    }
}
