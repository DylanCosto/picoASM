//! A check of each contig against the reads: the places it may be wrong.
//!
//! The reads show most of what can be wrong with an assembled contig without a reference: a
//! stretch few reads cover (a gap bridged on thin evidence, or a sample's other replicon run
//! in), and a place where many reads stop aligning with bases of theirs left over (a join that
//! is not in the genome, a locus the population holds in two versions, a repeat collapsed).
//! Each contig is cut into bins of [`BIN`] bases; a bin under [`LOW_DEPTH_DIVISOR`] of the
//! contig's median depth is shallow, and one where at least [`STOP_PERCENT`] of the reads
//! there stop is a stop. A contig with neither, and circular, is as the reads have it.

use picoasm_overlap::{Overlap, Reads};

/// Contigs are checked in bins of this many bases.
pub const BIN: u32 = 1_000;

/// A bin is shallow under the contig's median depth over this.
const LOW_DEPTH_DIVISOR: f64 = 4.0;

/// A read stops at a place when this many of its bases are left over past it.
const MIN_CLIP: u32 = 300;

/// A bin is a stop when this share (in percent) of the reads there stop in it...
const STOP_PERCENT: u64 = 20;

/// ...and at least this many do.
const MIN_STOPS: u64 = 5;

/// What the reads show of one contig.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Check {
    /// The depth of its shallowest bin.
    pub min_depth: f64,
    /// The median of its bins' depths: unlike the mean, not raised by a high-copy element
    /// inside the contig.
    pub median_depth: f64,
    /// Its bins (of [`BIN`] bases) under a quarter of its median depth.
    pub shallow_bins: usize,
    /// Its bins where many reads stop aligning (see the module docs).
    pub stop_bins: usize,
}

impl Check {
    /// Whether the reads show nothing wrong (a circular contig with no shallow bin and no stop).
    #[must_use]
    pub fn is_clean(&self, circular: bool) -> bool {
        circular && self.shallow_bins == 0 && self.stop_bins == 0
    }
}

/// Checks each of `seqs` against the reads' `mappings` onto them (see the module docs).
pub(crate) fn check(seqs: &[&[u8]], reads: &Reads, mappings: &[Option<Overlap>]) -> Vec<Check> {
    let bins = |len: usize| len.div_ceil(BIN as usize).max(1);
    let mut covered: Vec<Vec<u64>> = seqs.iter().map(|s| vec![0; bins(s.len())]).collect();
    let mut stops: Vec<Vec<u64>> = covered.clone();
    for (read, m) in reads.indices().zip(mappings) {
        let Some(m) = m else { continue };
        let c = m.target.index();
        let len = u32::try_from(seqs[c].len()).expect("contigs are shorter than 4 Gb");
        let (mut at, end) = (m.target_start, m.target_end);
        while at < end {
            let next = ((at / BIN + 1) * BIN).min(end);
            covered[c][(at / BIN) as usize] += u64::from(next - at);
            at = next;
        }
        let read_len = reads.read_len(read);
        let (left, right) = if m.strand.is_reverse() {
            (read_len - m.query_end, m.query_start)
        } else {
            (m.query_start, read_len - m.query_end)
        };
        // Not at the contig's own ends, where a read runs off a circle's cut or a line's end.
        if left >= MIN_CLIP && m.target_start >= MIN_CLIP {
            stops[c][(m.target_start / BIN) as usize] += 1;
        }
        if right >= MIN_CLIP && len - m.target_end >= MIN_CLIP {
            stops[c][((m.target_end - 1) / BIN) as usize] += 1;
        }
    }
    seqs.iter()
        .zip(covered)
        .zip(stops)
        .map(|((seq, covered), stops)| {
            let len = seq.len() as u64;
            // Each bin's mean depth; the last, shorter one, by its own length.
            #[expect(clippy::cast_precision_loss, reason = "mean depths, not exact counts")]
            let depth: Vec<f64> = covered
                .iter()
                .enumerate()
                .map(|(b, &n)| {
                    let bin_len = (len - b as u64 * u64::from(BIN)).min(u64::from(BIN));
                    n as f64 / bin_len.max(1) as f64
                })
                .collect();
            let mut sorted = depth.clone();
            sorted.sort_by(f64::total_cmp);
            let median = sorted[sorted.len() / 2];
            Check {
                min_depth: sorted[0],
                median_depth: median,
                shallow_bins: depth
                    .iter()
                    .filter(|&&d| d * LOW_DEPTH_DIVISOR < median)
                    .count(),
                stop_bins: depth
                    .iter()
                    .zip(&stops)
                    .filter(|&(&d, &s)| {
                        #[expect(
                            clippy::cast_precision_loss,
                            reason = "a count compared with a depth"
                        )]
                        let deep = s as f64 * 100.0 >= d * STOP_PERCENT as f64;
                        s >= MIN_STOPS && deep
                    })
                    .count(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use picoasm_overlap::{Edit, EditOp, ReadIndex};
    use picoasm_seq::Strand;

    use super::*;

    fn mapping(query: u32, query_len: u32, target_start: u32, target_end: u32) -> Overlap {
        Overlap {
            query: ReadIndex(query),
            target: ReadIndex(0),
            strand: Strand::Forward,
            query_start: 0,
            query_end: query_len.min(target_end - target_start),
            target_start,
            target_end,
            edits: vec![Edit::new(EditOp::Match, target_end - target_start)],
        }
    }

    #[test]
    fn a_shallow_stretch_and_a_place_reads_stop_are_found() {
        let seq = vec![b'A'; 10_000];
        let mut seqs = Vec::new();
        let mut mappings = Vec::new();
        // 20 reads across the contig, but none over 6-7 kb; and 8 reads that stop at 3.5 kb
        // with 1 kb left over.
        for i in 0..20 {
            seqs.push(vec![b'A'; 6_000]);
            mappings.push(Some(mapping(i, 6_000, 0, 6_000)));
            seqs.push(vec![b'A'; 3_000]);
            mappings.push(Some(mapping(i, 3_000, 7_000, 10_000)));
        }
        for i in 0..8 {
            seqs.push(vec![b'A'; 2_500]);
            mappings.push(Some(Overlap {
                query_end: 1_500,
                ..mapping(40 + i, 2_500, 2_000, 3_500)
            }));
        }
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let checks = check(&[seq.as_slice()], &reads, &mappings);
        assert_eq!(checks[0].shallow_bins, 1);
        assert_eq!(checks[0].stop_bins, 1);
        assert!(!checks[0].is_clean(true));
    }
}
