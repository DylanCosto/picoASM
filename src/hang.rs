//! Overhangs: overlaps between raw reads, snapped to the reads' ends.
//!
//! An Oxford Nanopore read's ends often hold a burst of errors or an adapter's remains, so an
//! overlap's alignment stops short of them ([`crate::overlap::Extension::Clipped`]). Left as it
//! is, such a read would be neither held inside another nor joined to the reads beyond its
//! noisy end, and every contig would stop at one. As miniasm does, an overlap whose alignment
//! stops within a short overhang of a read's end on both sides is taken to run to the ends:
//! on each side the alignment is carried on over the shorter of the two overhangs (its columns
//! counted as mismatches). One that leaves longer overhangs is two copies of a repeat inside
//! both reads, not an overlap, and is dropped.

use crate::overlap::{Edit, EditOp, Overlap, Reads};

/// Overhangs are snapped while together they are at most this many bases...
pub(crate) const MAX_HANG: u32 = 1_000;

/// ...and at most this share of the bases aligned.
pub(crate) const MAX_HANG_SHARE: f64 = 0.5;

/// Returns `o` carried on over its overhangs to the reads' ends, or `None` if they are too
/// long (see the module docs).
pub(crate) fn snap(reads: &Reads, o: &Overlap) -> Option<Overlap> {
    let query_len = reads.read_len(o.query);
    let target_len = reads.read_len(o.target);
    let (q_start, q_end) = if o.strand.is_reverse() {
        (query_len - o.query_end, query_len - o.query_start)
    } else {
        (o.query_start, o.query_end)
    };
    let left = q_start.min(o.target_start);
    let right = (query_len - q_end).min(target_len - o.target_end);
    let hang = left + right;
    let aligned = (o.target_end - o.target_start).min(q_end - q_start);
    if hang > MAX_HANG || f64::from(hang) > MAX_HANG_SHARE * f64::from(aligned) {
        return None;
    }
    let (q_start, q_end) = (q_start - left, q_end + right);
    let (query_start, query_end) = if o.strand.is_reverse() {
        (query_len - q_end, query_len - q_start)
    } else {
        (q_start, q_end)
    };
    let mut edits = Vec::with_capacity(o.edits.len() + 2);
    if left > 0 {
        edits.push(Edit::new(EditOp::Mismatch, left));
    }
    edits.extend_from_slice(&o.edits);
    if right > 0 {
        edits.push(Edit::new(EditOp::Mismatch, right));
    }
    Some(Overlap {
        query: o.query,
        target: o.target,
        strand: o.strand,
        query_start,
        query_end,
        target_start: o.target_start - left,
        target_end: o.target_end + right,
        edits,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::ReadIndex;
    use crate::seq::Strand;

    fn reads(lens: &[usize]) -> Reads {
        let seqs: Vec<Vec<u8>> = lens.iter().map(|&l| vec![b'A'; l]).collect();
        seqs.iter().map(Vec::as_slice).collect()
    }

    fn overlap(strand: Strand, q: (u32, u32), t: (u32, u32)) -> Overlap {
        Overlap {
            query: ReadIndex(0),
            target: ReadIndex(1),
            strand,
            query_start: q.0,
            query_end: q.1,
            target_start: t.0,
            target_end: t.1,
            edits: vec![Edit::new(EditOp::Match, t.1 - t.0)],
        }
    }

    #[test]
    fn short_overhangs_are_carried_to_the_ends() {
        let r = reads(&[5_000, 6_000]);
        // The query's 0..4,800 against the target's 1,200..6,000, less 100 at the query's
        // start and 200 at its end: a dovetail once snapped.
        let o = overlap(Strand::Forward, (100, 4_600), (1_300, 5_800));
        let s = snap(&r, &o).expect("short overhangs");
        assert_eq!((s.query_start, s.query_end), (0, 4_800));
        assert_eq!((s.target_start, s.target_end), (1_200, 6_000));
        assert_eq!(s.columns(), 4_800);
        assert_eq!(s.differences(), 300);
    }

    #[test]
    fn overhangs_follow_the_query_strand() {
        let r = reads(&[5_000, 6_000]);
        // Reverse: the oriented query is 5,000 - end .. 5,000 - start = 500..4,900.
        let o = overlap(Strand::Reverse, (100, 4_500), (0, 4_400));
        let s = snap(&r, &o).expect("short overhangs");
        // Left: min(500, 0) = 0. Right: min(5,000 - 4,900, 6,000 - 4,400) = 100.
        assert_eq!((s.query_start, s.query_end), (0, 4_500));
        assert_eq!((s.target_start, s.target_end), (0, 4_500));
    }

    #[test]
    fn a_repeat_inside_both_reads_is_dropped() {
        let r = reads(&[8_000, 8_000]);
        let o = overlap(Strand::Forward, (3_000, 5_000), (3_000, 5_000));
        assert_eq!(snap(&r, &o), None);
    }
}
