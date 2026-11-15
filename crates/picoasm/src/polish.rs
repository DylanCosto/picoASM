//! Polishing: the reads mapped onto the contigs, and a partial-order consensus per window.
//!
//! Each round maps every read onto the contigs ([`picoasm_overlap::map_reads`]), cuts the contigs
//! into windows of about [`BactParams::window_len`] bases, cuts each mapped read at the
//! windows' ends by its alignment, and replaces each window by the consensus of the read
//! segments that span it ([`PoaBackend`]), as Racon polishes. A window's ends are moved to a
//! place where the four bases around the cut all differ, so no homopolymer or two-base repeat
//! runs across it: reads disagree most on the lengths of those, and a cut inside one gives the
//! two windows' consensuses different shares of it.

use picoasm_consensus::{ConsensusBackend, PoaBackend, Window};
use picoasm_overlap::{EditOp, InvalidOverlapParamsError, Overlap, Reads, map_reads};
use picoasm_seq::{Strand, reverse_complement_into};

use crate::params::BactParams;

/// A contig being polished.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Draft {
    /// Its bases.
    pub(crate) seq: Vec<u8>,
    /// Whether its last base runs on into its first.
    pub(crate) circular: bool,
    /// How many bases at its start and at its end the last round could not polish (fewer
    /// than [`BactParams::min_depth`] reads), if it is linear: the one or two reads that run
    /// past a dead end in the layout. Their raw bases can keep the end from overlapping another
    /// draft; they can as well be what two ends overlap by, where few reads cross a gap.
    pub(crate) shallow_ends: [usize; 2],
    /// The mean depth of the reads mapped onto it in the last round (0 before any).
    pub(crate) depth: f64,
    /// Whether its start and its end are hairpins (see `hairpin.rs`): a linear replicon's ends.
    pub(crate) hairpins: [bool; 2],
}

impl Draft {
    /// A draft of `seq`, circular or not, with no shallow ends known.
    pub(crate) fn new(seq: Vec<u8>, circular: bool) -> Self {
        Self {
            seq,
            circular,
            shallow_ends: [0; 2],
            depth: 0.0,
            hairpins: [false; 2],
        }
    }

    /// Its bases without its shallow ends (all of them if those cover it).
    pub(crate) fn deep_part(&self) -> &[u8] {
        let [start, end] = self.shallow_ends;
        if start + end >= self.seq.len() {
            &self.seq
        } else {
            &self.seq[start..self.seq.len() - end]
        }
    }

    /// Moves the first `by` bases of a circular contig to its end.
    pub(crate) fn rotate(&mut self, by: usize) {
        if self.circular && !self.seq.is_empty() {
            let by = by % self.seq.len();
            self.seq.rotate_left(by);
        }
    }
}

/// What one round of polishing did.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct RoundStats {
    /// Reads mapped onto a contig.
    pub(crate) mapped_reads: usize,
    /// Contig bases the mapped reads cover, per contig.
    pub(crate) covered: Vec<u64>,
    /// Windows, and those with enough read segments to be polished.
    pub(crate) windows: usize,
    pub(crate) polished: usize,
}

/// What becomes of a window's bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum WindowFate {
    /// Too few reads: the draft's bases are kept.
    Kept,
    /// Most reads have no bases there: the window is removed.
    Emptied,
    /// The reads' consensus replaces them.
    Polished,
}

/// Polishes `drafts` once from `reads`, and notes the shallow ends of linear ones.
///
/// # Errors
///
/// Returns an error if the mapping settings cannot be used.
pub(crate) fn polish_round(
    drafts: &mut [Draft],
    reads: &Reads,
    params: &BactParams,
    backend: &PoaBackend,
) -> Result<RoundStats, InvalidOverlapParamsError> {
    let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
    let mappings = map_reads(&targets, reads, &params.map)?;
    let boundaries: Vec<Vec<u32>> = drafts
        .iter()
        .map(|d| window_boundaries(&d.seq, params.window_len))
        .collect();
    let mut segments: Vec<Vec<Vec<Vec<u8>>>> = boundaries
        .iter()
        .map(|b| vec![Vec::new(); b.len().saturating_sub(1)])
        .collect();
    let mut stats = RoundStats {
        covered: vec![0; drafts.len()],
        ..RoundStats::default()
    };
    let mut oriented = Vec::new();
    for (read, mapping) in reads.indices().zip(&mappings) {
        let Some(m) = mapping else { continue };
        stats.mapped_reads += 1;
        let contig = m.target.index();
        stats.covered[contig] += u64::from(m.target_end - m.target_start);
        let bases = reads.get(read);
        let oriented_start = match m.strand {
            Strand::Forward => {
                oriented.clear();
                oriented.extend_from_slice(bases);
                m.query_start
            }
            Strand::Reverse => {
                reverse_complement_into(bases, &mut oriented);
                reads.read_len(read) - m.query_end
            }
        };
        let cuts = cut_points(m, oriented_start, &boundaries[contig]);
        for pair in cuts.windows(2) {
            let ((window, from), (next, to)) = (pair[0], pair[1]);
            debug_assert_eq!(next, window + 1, "cut points are consecutive");
            let held = &mut segments[contig][window];
            if held.len() < params.max_depth {
                held.push(oriented[pos(from)..pos(to)].to_vec());
            }
        }
    }

    // One batch for all contigs; windows with too few segments keep the draft's bases, and
    // windows most of whose reads have no bases there (a stretch a join or a chimera put in
    // twice) lose theirs: the consensus leaves empty segments out, so it would keep them.
    let mut batch: Vec<Window> = Vec::new();
    let mut fates: Vec<Vec<WindowFate>> = Vec::with_capacity(segments.len());
    let mut deep: Vec<(usize, usize)> = Vec::with_capacity(segments.len());
    for windows in segments {
        stats.windows += windows.len();
        deep.push(deep_enough(&windows, params.min_depth));
        let mut these = Vec::with_capacity(windows.len());
        for held in windows {
            let empty = held.iter().filter(|s| s.is_empty()).count();
            if held.len() < params.min_depth {
                these.push(WindowFate::Kept);
            } else if 2 * empty > held.len() {
                these.push(WindowFate::Emptied);
            } else {
                these.push(WindowFate::Polished);
                batch.push(Window::new(held));
            }
        }
        fates.push(these);
    }
    stats.polished = batch.len();
    let mut consensus = backend.consensus(&batch).into_iter();
    drop(batch);
    for (((draft, fates), bounds), &(first, last)) in
        drafts.iter_mut().zip(&fates).zip(&boundaries).zip(&deep)
    {
        let mut seq = Vec::with_capacity(draft.seq.len() + draft.seq.len() / 100);
        // Where the deep windows start and end in the new bases.
        let (mut deep_start, mut deep_end) = (0, 0);
        for (w, fate) in fates.iter().enumerate() {
            if w == first {
                deep_start = seq.len();
            }
            let piece = &draft.seq[pos(bounds[w])..pos(bounds[w + 1])];
            match fate {
                WindowFate::Kept => seq.extend_from_slice(piece),
                WindowFate::Emptied => {}
                WindowFate::Polished => match consensus.next() {
                    Some(c) if !c.is_empty() => seq.extend_from_slice(&c),
                    _ => seq.extend_from_slice(piece),
                },
            }
            if w + 1 == last {
                deep_end = seq.len();
            }
        }
        // None, too, if every deep window was emptied: there are no polished bases to keep.
        draft.shallow_ends = if draft.circular || last <= first || deep_end <= deep_start {
            [0; 2]
        } else {
            [deep_start, seq.len() - deep_end]
        };
        draft.seq = seq;
    }
    for (draft, &covered) in drafts.iter_mut().zip(&stats.covered) {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a mean depth, not an exact count"
        )]
        let depth = covered as f64 / draft.seq.len().max(1) as f64;
        draft.depth = depth;
    }
    Ok(stats)
}

/// Returns the first window of a draft with enough read segments to be polished and the one
/// after the last (see [`Draft::shallow_ends`]), given each window's segments; all of them if
/// none has.
fn deep_enough(windows: &[Vec<Vec<u8>>], min_depth: usize) -> (usize, usize) {
    match (
        windows.iter().position(|w| w.len() >= min_depth),
        windows.iter().rposition(|w| w.len() >= min_depth),
    ) {
        (Some(first), Some(last)) => (first, last + 1),
        // None deep enough: the draft is left whole, for the depth filters to judge.
        _ => (0, windows.len()),
    }
}

/// Returns the window ends of `seq`: 0, then about every `window_len` bases a place where the
/// two bases either side of the cut are four different bases (or the next such place within
/// `window_len / 2`), then the length.
pub(crate) fn window_boundaries(seq: &[u8], window_len: usize) -> Vec<u32> {
    let mut out = vec![0];
    let len = seq.len();
    let mut at = window_len;
    while at + window_len / 2 < len {
        let clean = (at..(at + window_len / 2).min(len - 2)).find(|&p| {
            let w = &seq[p - 2..p + 2];
            w[0] != w[1]
                && w[0] != w[2]
                && w[0] != w[3]
                && w[1] != w[2]
                && w[1] != w[3]
                && w[2] != w[3]
        });
        let cut = clean.unwrap_or(at);
        out.push(u32::try_from(cut).expect("contigs are shorter than 4 Gb"));
        at = cut + window_len;
    }
    if len > 0 {
        out.push(u32::try_from(len).expect("contigs are shorter than 4 Gb"));
    }
    out
}

/// Returns, for each of `boundaries` (sorted target positions) that the mapping `m` spans, its
/// index and the position in the oriented read aligned to it. Bases the read has and the
/// contig lacks just before a boundary fall in the window before it.
pub(crate) fn cut_points(
    m: &Overlap,
    oriented_start: u32,
    boundaries: &[u32],
) -> Vec<(usize, u32)> {
    let mut i = boundaries.partition_point(|&b| b < m.target_start);
    let (mut t, mut q) = (m.target_start, oriented_start);
    let mut out = Vec::new();
    for e in &m.edits {
        let n = e.len();
        match e.op() {
            EditOp::Match | EditOp::Mismatch => {
                while i < boundaries.len() && boundaries[i] < t + n {
                    out.push((i, q + (boundaries[i] - t)));
                    i += 1;
                }
                t += n;
                q += n;
            }
            EditOp::Insertion => q += n,
            EditOp::Deletion => {
                while i < boundaries.len() && boundaries[i] < t + n {
                    out.push((i, q));
                    i += 1;
                }
                t += n;
            }
        }
    }
    if i < boundaries.len() && boundaries[i] == t {
        out.push((i, q));
    }
    out
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;
    use picoasm_overlap::{Edit, ReadIndex};

    #[test]
    fn shallow_windows_at_a_linear_drafts_ends_are_trimmed() {
        let windows: Vec<Vec<Vec<u8>>> = [1, 2, 30, 40, 2, 40, 9, 3, 1]
            .iter()
            .map(|&n| vec![b"ACGT".to_vec(); n])
            .collect();
        // Windows under 3 segments at either end go; the 2 inside stays.
        assert_eq!(deep_enough(&windows, 3), (2, 8));
        let none: Vec<Vec<Vec<u8>>> = vec![vec![b"A".to_vec()]; 3];
        assert_eq!(deep_enough(&none, 3), (0, 3));
    }

    #[test]
    fn windows_are_cut_where_four_bases_differ() {
        let seq = b"AAAACCCCACGTTTTTGGGGACGTAAAA".repeat(40);
        let bounds = window_boundaries(&seq, 100);
        assert_eq!(bounds[0], 0);
        assert_eq!(*bounds.last().unwrap() as usize, seq.len());
        for &b in &bounds[1..bounds.len() - 1] {
            let w = &seq[b as usize - 2..b as usize + 2];
            let mut sorted = w.to_vec();
            sorted.sort_unstable();
            sorted.dedup();
            assert_eq!(sorted.len(), 4, "cut at {b} inside {w:?}");
        }
        assert!(bounds.windows(2).all(|p| p[1] > p[0]));
    }

    #[test]
    fn cut_points_follow_the_alignment() {
        // Target 100..120 aligned to the read: 5 matches, 2 inserted, 5 matches, 3 deleted,
        // 7 matches.
        let m = Overlap {
            query: ReadIndex(0),
            target: ReadIndex(0),
            strand: Strand::Forward,
            query_start: 10,
            query_end: 29,
            target_start: 100,
            target_end: 120,
            edits: vec![
                Edit::new(EditOp::Match, 5),
                Edit::new(EditOp::Insertion, 2),
                Edit::new(EditOp::Match, 5),
                Edit::new(EditOp::Deletion, 3),
                Edit::new(EditOp::Match, 7),
            ],
        };
        let cuts = cut_points(&m, 10, &[90, 100, 105, 111, 120, 130]);
        // 100 -> 10; 105 -> after the insertion, 17; 111 -> 22; 111 is the second deleted
        // base, so the read position is where the deletion starts, 22; 120 -> the end, 29.
        assert_eq!(cuts, [(1, 10), (2, 17), (3, 22), (4, 29)]);
    }
}
