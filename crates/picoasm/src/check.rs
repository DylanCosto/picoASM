//! A check of each contig against the reads: the places it may be wrong.
//!
//! The reads show most of what can be wrong with an assembled contig without a reference: a
//! stretch few reads cover (a gap bridged on thin evidence, or a sample's other replicon run
//! in), and a place where many reads stop aligning with bases of theirs left over (a join that
//! is not in the genome, a locus the population holds in two versions, a repeat collapsed).
//! Each contig is cut into bins of [`BIN`] bases, the last taking the remainder (so 1-2 kb: a
//! bin of a few bases at a contig's end has too few reads to measure); a bin under
//! [`LOW_DEPTH_DIVISOR`] of the
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
    /// Its bins where the reads show something wrong, in order along the contig; a bin both
    /// shallow and a stop is listed once for each.
    pub flagged: Vec<FlaggedBin>,
}

impl Check {
    /// Whether the reads show nothing wrong (a circular contig with no shallow bin and no stop).
    #[must_use]
    pub fn is_clean(&self, circular: bool) -> bool {
        circular && self.flagged.is_empty()
    }

    /// Its bins (see the module docs) under a quarter of its median depth.
    #[must_use]
    pub fn shallow_bins(&self) -> usize {
        self.count(Flag::Shallow)
    }

    /// Its bins where many reads stop aligning (see the module docs).
    #[must_use]
    pub fn stop_bins(&self) -> usize {
        self.count(Flag::Stop)
    }

    /// Moves its flagged bins with the contig's bases, after changes to them: `moved` takes a
    /// place on the contig before the changes (a base, or the end) to where it is now.
    pub(crate) fn lift(&mut self, moved: impl Fn(u32) -> u32) {
        for bin in &mut self.flagged {
            bin.start = moved(bin.start);
            bin.end = moved(bin.end);
        }
    }

    fn count(&self, flag: Flag) -> usize {
        self.flagged.iter().filter(|b| b.flag == flag).count()
    }
}

/// What the reads show wrong in one bin.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Flag {
    /// Under a quarter of the contig's median depth: a gap bridged on thin evidence, or a
    /// sample's other replicon run in.
    Shallow,
    /// Many reads stop aligning here with bases of theirs left over: a join not in the
    /// genome, a locus the population holds in two versions, or a repeat collapsed.
    Stop,
}

/// One bin of a contig the reads flag (see [`Check::flagged`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlaggedBin {
    /// Its first base (0-based) on the contig; it moves with the bases when later base calls
    /// or short-read polishing change them.
    pub start: u32,
    /// One past its last base.
    pub end: u32,
    /// What is wrong.
    pub flag: Flag,
    /// The mean depth of the reads over it.
    pub depth: f64,
    /// The reads that stop aligning in it (see the module docs).
    pub stops: u64,
}

/// Checks each of `seqs` against the reads' `mappings` onto them (see the module docs).
pub(crate) fn check(seqs: &[&[u8]], reads: &Reads, mappings: &[Option<Overlap>]) -> Vec<Check> {
    // The last bin takes the remainder of the contig past the last whole bin.
    let bins = |len: usize| (len / BIN as usize).max(1);
    let mut covered: Vec<Vec<u64>> = seqs.iter().map(|s| vec![0; bins(s.len())]).collect();
    let mut stops: Vec<Vec<u64>> = covered.clone();
    for (read, m) in reads.indices().zip(mappings) {
        let Some(m) = m else { continue };
        let c = m.target.index();
        let len = u32::try_from(seqs[c].len()).expect("contigs are shorter than 4 Gb");
        let last = covered[c].len() - 1;
        let bin = |at: u32| ((at / BIN) as usize).min(last);
        let (mut at, end) = (m.target_start, m.target_end);
        while at < end {
            let next = ((at / BIN + 1) * BIN).min(end);
            covered[c][bin(at)] += u64::from(next - at);
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
            stops[c][bin(m.target_start)] += 1;
        }
        if right >= MIN_CLIP && len - m.target_end >= MIN_CLIP {
            stops[c][bin(m.target_end - 1)] += 1;
        }
    }
    seqs.iter()
        .zip(covered)
        .zip(stops)
        .map(|((seq, covered), stops)| {
            let len = seq.len() as u64;
            let last = covered.len() - 1;
            // Each bin's place on the contig; the last runs to its end.
            let span = |b: usize| {
                let start = b as u64 * u64::from(BIN);
                let end = if b == last {
                    len
                } else {
                    start + u64::from(BIN)
                };
                let place = |x: u64| u32::try_from(x).expect("contigs are shorter than 4 Gb");
                (place(start), place(end))
            };
            // Each bin's mean depth, by its own length.
            #[expect(clippy::cast_precision_loss, reason = "mean depths, not exact counts")]
            let depth: Vec<f64> = covered
                .iter()
                .enumerate()
                .map(|(b, &n)| {
                    let (start, end) = span(b);
                    n as f64 / f64::from((end - start).max(1))
                })
                .collect();
            let mut sorted = depth.clone();
            sorted.sort_by(f64::total_cmp);
            let median = sorted[sorted.len() / 2];
            let mut flagged = Vec::new();
            for (b, (&d, &s)) in depth.iter().zip(&stops).enumerate() {
                let (start, end) = span(b);
                let bin = |flag| FlaggedBin {
                    start,
                    end,
                    flag,
                    depth: d,
                    stops: s,
                };
                if d * LOW_DEPTH_DIVISOR < median {
                    flagged.push(bin(Flag::Shallow));
                }
                #[expect(clippy::cast_precision_loss, reason = "a count compared with a depth")]
                let deep = s as f64 * 100.0 >= d * STOP_PERCENT as f64;
                if s >= MIN_STOPS && deep {
                    flagged.push(bin(Flag::Stop));
                }
            }
            Check {
                min_depth: sorted[0],
                median_depth: median,
                flagged,
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
        assert_eq!(checks[0].shallow_bins(), 1);
        assert_eq!(checks[0].stop_bins(), 1);
        let at: Vec<(u32, u32, Flag)> = checks[0]
            .flagged
            .iter()
            .map(|b| (b.start, b.end, b.flag))
            .collect();
        assert_eq!(
            at,
            [(3_000, 4_000, Flag::Stop), (6_000, 7_000, Flag::Shallow)]
        );
        assert_eq!(checks[0].flagged[0].stops, 8);
        assert!(!checks[0].is_clean(true));
    }

    /// Reads of `len` bases laid end to end every `step` bases along a contig of `contig_len`,
    /// each mapped whole.
    fn tiled(contig_len: u32, len: u32, step: u32) -> (Vec<Vec<u8>>, Vec<Option<Overlap>>) {
        let mut seqs = Vec::new();
        let mut mappings = Vec::new();
        let mut start = 0;
        while start < contig_len {
            let end = (start + len).min(contig_len);
            let i = u32::try_from(seqs.len()).expect("few reads");
            seqs.push(vec![b'A'; (end - start) as usize]);
            mappings.push(Some(mapping(i, end - start, start, end)));
            start += step;
        }
        (seqs, mappings)
    }

    fn checked(contig_len: u32, seqs: &[Vec<u8>], mappings: &[Option<Overlap>]) -> Check {
        let seq = vec![b'A'; contig_len as usize];
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        check(&[seq.as_slice()], &reads, mappings).remove(0)
    }

    #[test]
    fn a_few_bases_past_the_last_whole_bin_are_part_of_it() {
        // 10,003 bases at 20x but the last 3 at 2x: alone they would be a shallow bin.
        let (mut seqs, mut mappings) = tiled(10_000, 2_000, 100);
        for i in 0..2 {
            seqs.push(vec![b'A'; 503]);
            mappings.push(Some(mapping(1_000 + i, 503, 9_500, 10_003)));
        }
        let c = checked(10_003, &seqs, &mappings);
        assert!(c.flagged.is_empty(), "{:?}", c.flagged);
    }

    #[test]
    fn the_last_bin_runs_to_the_contig_end_when_flagged() {
        // Reads stop short of the last 1.5 kb of 5.5 kb.
        let (seqs, mappings) = tiled(4_000, 1_000, 50);
        let c = checked(5_500, &seqs, &mappings);
        let at: Vec<(u32, u32, Flag)> =
            c.flagged.iter().map(|b| (b.start, b.end, b.flag)).collect();
        assert_eq!(at, [(4_000, 5_500, Flag::Shallow)]);
    }

    #[test]
    fn contigs_shorter_than_a_bin_or_empty_are_one_bin() {
        let (seqs, mappings) = tiled(700, 700, 70);
        let c = checked(700, &seqs, &mappings);
        assert!(c.flagged.is_empty());
        assert!(c.median_depth > 0.0);
        let empty = checked(0, &[], &[]);
        assert!(empty.flagged.is_empty());
        assert!(empty.median_depth.abs() < f64::EPSILON);
    }

    #[test]
    fn moved_bins_keep_their_order_and_flags() {
        let (seqs, mappings) = tiled(4_000, 1_000, 50);
        let mut c = checked(5_500, &seqs, &mappings);
        c.lift(|at| at + at / 1_000);
        assert_eq!((c.flagged[0].start, c.flagged[0].end), (4_004, 5_505));
        assert_eq!(c.shallow_bins(), 1);
    }
}
