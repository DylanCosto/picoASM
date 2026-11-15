//! Read trimming: each raw read cut to the stretch its overlaps vouch for.
//!
//! Oxford Nanopore reads can carry kilobases of sequence that matches nothing: a stretch of
//! very low quality (errors at 20-40%), the remains of an adapter, or a second molecule read
//! on. An overlap's alignment stops where such a stretch starts, so the read seems to overlap
//! its neighbours only part way, and once its overhang is too long to snap
//! ([`crate::hang`]) the overlap is lost, and with it the join between two contigs. As miniasm
//! does, each read is first cut to its longest stretch that at least [`MIN_COVER`] overlaps
//! cover (fewer below 24x, [`min_cover`]), and every overlap is cut to the two reads'
//! stretches; then the overhangs left are the short ones a snap can carry.

use picoasm_overlap::{Edit, EditOp, Overlap, Reads};

/// A read is kept where at least this many overlaps cover it...
pub(crate) const MIN_COVER: u32 = 3;

/// ...or one per this much depth, if fewer: at 12x a stretch only one or two reads cross has
/// fewer than three overlaps over it, and cutting the reads' ends there left gaps between
/// contigs (*L. innocua* at 12x: 50 kb left out with three, 6 kb with one).
const DEPTH_PER_COVER: f64 = 8.0;

/// The overlaps that must cover a stretch of a read for it to be kept, at `depth`.
pub(crate) fn min_cover(depth: f64) -> u32 {
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a small positive count"
    )]
    let by_depth = (depth / DEPTH_PER_COVER).floor().max(1.0) as u32;
    by_depth.min(MIN_COVER)
}

/// The reads cut to their stretches, and where each stretch lies in the read (empty for a read
/// with none at least `min_len` long, which is left out of the layout).
#[derive(Clone, Debug)]
pub(crate) struct Trimmed {
    pub(crate) reads: Reads,
    pub(crate) spans: Vec<(u32, u32)>,
    /// Each read's length before it was cut.
    lens: Vec<u32>,
}

/// Cuts each of `reads` to its longest stretch that `min_cover` of `overlaps` cover, if it is
/// at least `min_len` bases long.
pub(crate) fn trim(reads: &Reads, overlaps: &[Overlap], min_len: u32, min_cover: u32) -> Trimmed {
    let mut events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); reads.len()];
    for o in overlaps {
        events[o.query.index()].extend([(o.query_start, 1), (o.query_end, -1)]);
        events[o.target.index()].extend([(o.target_start, 1), (o.target_end, -1)]);
    }
    let spans: Vec<(u32, u32)> = events
        .iter_mut()
        .map(|e| {
            // Starts before ends at one place, so an overlap that starts where another ends
            // keeps the cover unbroken.
            e.sort_unstable_by_key(|&(at, step)| (at, -step));
            let mut best = (0, 0);
            let (mut depth, mut from) = (0_i32, None);
            for &(at, step) in e.iter() {
                depth += step;
                let covered = u32::try_from(depth).is_ok_and(|d| d >= min_cover);
                match (covered, from) {
                    (true, None) => from = Some(at),
                    (false, Some(start)) => {
                        if at - start > best.1 - best.0 {
                            best = (start, at);
                        }
                        from = None;
                    }
                    _ => {}
                }
            }
            if best.1 - best.0 >= min_len {
                best
            } else {
                (0, 0)
            }
        })
        .collect();
    let mut trimmed = Reads::with_capacity(reads.len(), reads.total_bases());
    for (read, &(start, end)) in reads.indices().zip(&spans) {
        trimmed.push(&reads.get(read)[pos(start)..pos(end)]);
    }
    Trimmed {
        reads: trimmed,
        spans,
        lens: reads.indices().map(|r| reads.read_len(r)).collect(),
    }
}

impl Trimmed {
    /// Returns `o` cut to its reads' stretches and placed on the cut reads, or `None` if less
    /// than `min_overlap` bases of either read are left. Bases are taken off both reads alike,
    /// so the cut ends are off by the indels between them; the alignment is kept only as its
    /// length (one run of matches, and the reads' difference in length as an indel).
    pub(crate) fn cut(&self, o: &Overlap, min_overlap: u32) -> Option<Overlap> {
        let (q_from, q_to) = self.spans[o.query.index()];
        let (t_from, t_to) = self.spans[o.target.index()];
        if q_to == q_from || t_to == t_from {
            return None;
        }
        let query_len = self.lens[o.query.index()];
        // The query in the orientation it is aligned in, and its stretch there.
        let (mut qs, mut qe) =
            oriented(o.strand.is_reverse(), query_len, o.query_start, o.query_end);
        let (span_s, span_e) = oriented(o.strand.is_reverse(), query_len, q_from, q_to);
        let (mut ts, mut te) = (o.target_start, o.target_end);
        let left = span_s.saturating_sub(qs).max(t_from.saturating_sub(ts));
        let right = qe.saturating_sub(span_e).max(te.saturating_sub(t_to));
        qs += left;
        ts += left;
        qe = qe.saturating_sub(right);
        te = te.saturating_sub(right);
        if qe < qs + min_overlap || te < ts + min_overlap {
            return None;
        }
        // Onto the cut reads.
        let (qs, qe) = (qs - span_s, qe - span_s);
        let (ts, te) = (ts - t_from, te - t_from);
        let cut_len = q_to - q_from;
        let (query_start, query_end) = oriented(o.strand.is_reverse(), cut_len, qs, qe);
        let (q_cols, t_cols) = (qe - qs, te - ts);
        let mut edits = vec![Edit::new(EditOp::Match, q_cols.min(t_cols))];
        match q_cols.cmp(&t_cols) {
            std::cmp::Ordering::Greater => {
                edits.push(Edit::new(EditOp::Insertion, q_cols - t_cols));
            }
            std::cmp::Ordering::Less => edits.push(Edit::new(EditOp::Deletion, t_cols - q_cols)),
            std::cmp::Ordering::Equal => {}
        }
        Some(Overlap {
            query: o.query,
            target: o.target,
            strand: o.strand,
            query_start,
            query_end,
            target_start: ts,
            target_end: te,
            edits,
        })
    }
}

/// Converts `start..end` on a read of `len` bases between its forward strand and its reverse
/// complement (the same either way).
fn oriented(reverse: bool, len: u32, start: u32, end: u32) -> (u32, u32) {
    if reverse {
        (len - end, len - start)
    } else {
        (start, end)
    }
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;
    use picoasm_overlap::ReadIndex;
    use picoasm_seq::Strand;

    fn overlap(query: u32, target: u32, strand: Strand, q: (u32, u32), t: (u32, u32)) -> Overlap {
        Overlap {
            query: ReadIndex(query),
            target: ReadIndex(target),
            strand,
            query_start: q.0,
            query_end: q.1,
            target_start: t.0,
            target_end: t.1,
            edits: vec![Edit::new(EditOp::Match, t.1 - t.0)],
        }
    }

    #[test]
    fn reads_are_cut_to_their_covered_stretch() {
        let seqs = [vec![b'A'; 10_000], vec![b'C'; 10_000]];
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        // Read 0 is covered three times over 1,000..6,000, twice past it.
        let o = [
            overlap(0, 1, Strand::Forward, (1_000, 6_000), (0, 5_000)),
            overlap(0, 1, Strand::Forward, (500, 7_000), (0, 6_500)),
            overlap(0, 1, Strand::Forward, (1_000, 8_000), (0, 7_000)),
        ];
        let t = trim(&reads, &o, 1_000, MIN_COVER);
        assert_eq!(t.spans[0], (1_000, 6_000));
        assert_eq!(t.spans[1], (0, 5_000));
        assert_eq!(t.reads.read_len(ReadIndex(0)), 5_000);
        assert_eq!(t.reads.get(ReadIndex(1))[0], b'C');
    }

    #[test]
    fn shallow_reads_need_fewer_overlaps_to_keep_a_stretch() {
        assert_eq!(min_cover(50.0), MIN_COVER);
        assert_eq!(min_cover(24.0), 3);
        assert_eq!(min_cover(20.0), 2);
        assert_eq!(min_cover(12.0), 1);
        assert_eq!(min_cover(3.0), 1);
    }

    #[test]
    fn overlaps_are_cut_on_both_reads_alike() {
        let seqs = [vec![b'A'; 10_000], vec![b'C'; 8_000]];
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let t = Trimmed {
            reads,
            spans: vec![(0, 9_000), (500, 8_000)],
            lens: vec![10_000, 8_000],
        };
        // The query's 4,000..10,000 against the target's 0..6,000: the query loses its last
        // 1,000 bases and the target its first 500.
        let o = overlap(0, 1, Strand::Forward, (4_000, 10_000), (0, 6_000));
        let c = t.cut(&o, 500).expect("long enough");
        assert_eq!((c.query_start, c.query_end), (4_500, 9_000));
        assert_eq!((c.target_start, c.target_end), (0, 4_500));
        // Reverse: the query's forward 0..6,000 is its oriented 4,000..10,000.
        let o = overlap(0, 1, Strand::Reverse, (0, 6_000), (0, 6_000));
        let c = t.cut(&o, 500).expect("long enough");
        // Oriented, the query's stretch is 1,000..10,000 and the overlap 4,000..10,000: only
        // the target's first 500 bases go, and 500 of the query's with them.
        assert_eq!((c.target_start, c.target_end), (0, 5_500));
        assert_eq!((c.query_start, c.query_end), (0, 5_500));
    }
}
