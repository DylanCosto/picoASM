//! The CPU partial-order alignment backend.

use rayon::prelude::*;

use crate::consensus::align::{Aligner, MAX_DP_BYTES};
use crate::consensus::graph::{PoaGraph, normalize_base};
use crate::consensus::{ConsensusBackend, InvalidPoaParamsError, PoaParams, Window};

/// A CPU consensus backend: banded partial-order alignment with a heaviest-path consensus.
///
/// See the [crate documentation](crate) for the algorithm, its cost and its limits. Batches are
/// processed in parallel on the caller's rayon pool (the global pool unless the call runs
/// inside [`rayon::ThreadPool::install`]); the output does not depend on the pool size.
///
/// The consensus of a window depends on the order of its segments (it picks the backbone and
/// breaks ties by it), so callers that need reproducible output must hand segments over in a
/// canonical order, as [`Window`] explains.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct PoaBackend {
    params: PoaParams,
}

impl PoaBackend {
    /// Creates a backend with validated parameters.
    ///
    /// # Errors
    ///
    /// Returns the error from [`PoaParams::validate`] if the parameters are unusable.
    pub const fn new(params: PoaParams) -> Result<Self, InvalidPoaParamsError> {
        match params.validate() {
            Ok(()) => Ok(Self { params }),
            Err(err) => Err(err),
        }
    }

    /// The most memory the alignment matrices of one thread ever hold, in bytes (48 MiB).
    ///
    /// A segment whose alignment would need more is aligned with a narrower band, or left out.
    pub const MATRIX_BYTES_PER_THREAD: usize = MAX_DP_BYTES;

    /// A conservative estimate of the other working memory per base of a window, in bytes.
    ///
    /// It covers the alignment graph (at most one node per segment base), its per-node
    /// buffers and the normalised segment copies. Unlike the matrices, this memory is not
    /// capped: it grows with the largest window a thread processes.
    pub const BYTES_PER_WINDOW_BASE: usize = 512;

    /// Returns the parameters.
    #[must_use]
    pub const fn params(&self) -> &PoaParams {
        &self.params
    }

    /// Estimates the peak memory one worker thread uses, in bytes, for windows of at most
    /// `window_bases` bases in total (summed over their segments).
    ///
    /// This is [`PoaBackend::MATRIX_BYTES_PER_THREAD`] plus
    /// [`PoaBackend::BYTES_PER_WINDOW_BASE`] per base. A batch holds one set of buffers per
    /// rayon worker, so a budget for a whole batch is this figure times the pool size, plus the
    /// windows and outputs themselves. The figure is an estimate for budgeting, not a hard
    /// limit.
    #[must_use]
    pub const fn peak_memory_per_thread(window_bases: usize) -> usize {
        Self::MATRIX_BYTES_PER_THREAD
            .saturating_add(window_bases.saturating_mul(Self::BYTES_PER_WINDOW_BASE))
    }

    /// Returns the consensus of one window; see [`PoaBackend::consensus_segments`].
    #[must_use]
    pub fn consensus_window(&self, window: &Window) -> Vec<u8> {
        self.consensus_segments(window.segments())
    }

    /// Returns the consensus of one window given as a slice of segments, on the calling thread.
    ///
    /// The output is upper-case `ACGT`, with `N` only where the consensus path runs through a
    /// non-`ACGT` base. Empty segments are ignored; with no other segment the output is empty,
    /// and with exactly one it is that segment, normalised.
    #[must_use]
    pub fn consensus_segments<S: AsRef<[u8]>>(&self, segments: &[S]) -> Vec<u8> {
        Workspace::default().consensus(segments, &self.params)
    }
}

impl ConsensusBackend for PoaBackend {
    fn consensus(&self, windows: &[Window]) -> Vec<Vec<u8>> {
        windows
            .par_iter()
            .map_init(Workspace::default, |workspace, window| {
                workspace.consensus(window.segments(), &self.params)
            })
            .collect()
    }
}

/// Buffers reused across the windows of one rayon job (one split of the batch).
///
/// [`ConsensusBackend::consensus`] builds one per job split, not one per worker thread, but a
/// thread runs one job at a time, so it holds at most one live workspace.
#[derive(Debug, Default)]
struct Workspace {
    graph: PoaGraph,
    aligner: Aligner,
    /// Normalised copies of the current window's non-empty segments (and spare buffers).
    segments: Vec<Vec<u8>>,
}

impl Workspace {
    fn consensus<S: AsRef<[u8]>>(&mut self, segments: &[S], params: &PoaParams) -> Vec<u8> {
        let Self {
            graph,
            aligner,
            segments: normalized,
        } = self;
        let mut count = 0;
        for segment in segments {
            let segment = segment.as_ref();
            if segment.is_empty() {
                continue;
            }
            if count == normalized.len() {
                normalized.push(Vec::new());
            }
            let buf = &mut normalized[count];
            buf.clear();
            buf.extend(segment.iter().map(|&base| normalize_base(base)));
            count += 1;
        }
        let segments = &normalized[..count];
        if segments.len() <= 1 {
            return segments.first().cloned().unwrap_or_default();
        }

        let order = alignment_order(segments);
        let backbone = &segments[order[0]];
        graph.reset(backbone);
        for &i in &order[1..] {
            if !graph.topo_sort() {
                log::warn!("partial-order graph has a cycle; using the median segment");
                return backbone.clone();
            }
            match aligner.align(graph, &segments[i], params) {
                Some(path) => graph.add_alignment(&segments[i], path),
                None => log::debug!(
                    "segment of {} bases does not fit the alignment matrices; skipped",
                    segments[i].len()
                ),
            }
        }
        if !graph.topo_sort() {
            log::warn!("partial-order graph has a cycle; using the median segment");
            return backbone.clone();
        }
        graph.consensus()
    }
}

/// Returns the order in which segments are aligned into the graph.
///
/// The first, the backbone, is the segment at the lower-median position when the segments are
/// sorted by (length, index), so a truncated or overlong segment never becomes the backbone.
/// Among equal-length segments this is whichever falls on that position, not necessarily the
/// lowest index. The others follow by increasing length difference from the backbone, ties to
/// the lower index, so outliers come last.
fn alignment_order(segments: &[Vec<u8>]) -> Vec<usize> {
    let mut by_len: Vec<usize> = (0..segments.len()).collect();
    by_len.sort_by_key(|&i| (segments[i].len(), i));
    let backbone = by_len[(by_len.len() - 1) / 2];
    let backbone_len = segments[backbone].len();
    let mut rest: Vec<usize> = by_len.into_iter().filter(|&i| i != backbone).collect();
    rest.sort_by_key(|&i| (segments[i].len().abs_diff(backbone_len), i));
    let mut order = Vec::with_capacity(segments.len());
    order.push(backbone);
    order.extend(rest);
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orders_from_the_median_outwards() {
        let segments: Vec<Vec<u8>> = [5, 1, 4, 4, 9].iter().map(|&n| vec![b'A'; n]).collect();
        // Sorted by (len, index): 1, 2, 3, 0, 4; the median is index 3 (length 4).
        assert_eq!(alignment_order(&segments), vec![3, 2, 0, 1, 4]);
        let two = vec![b"AC".to_vec(), b"A".to_vec()];
        assert_eq!(alignment_order(&two), vec![1, 0]);
    }

    #[test]
    fn two_segments_give_the_backbone() {
        // Equal lengths, so the first segment is the backbone. Each has one base the other
        // lacks; the even split keeps the backbone's allele at both sites.
        let first = b"ACGTTGCAGTACCATGCAAG".to_vec();
        let second = b"ACGTTGCATACCATGCTAAG".to_vec();
        let backend = PoaBackend::default();
        assert_eq!(
            backend.consensus_segments(&[first.clone(), second.clone()]),
            first
        );
        assert_eq!(backend.consensus_segments(&[second.clone(), first]), second);
    }

    #[test]
    fn an_even_indel_split_takes_the_backbone_allele() {
        let short = b"ACGTTGCATGCAAGTCCGAT".to_vec();
        let long = b"ACGTTGCATGGCAAGTCCGAT".to_vec();
        let backend = PoaBackend::default();
        // Sorted by length, the lower median is a short segment, so the backbone is short.
        let segments = [long.clone(), short.clone(), long.clone(), short.clone()];
        assert_eq!(alignment_order(&segments)[0], 3);
        assert_eq!(backend.consensus_segments(&segments), short);
        // One more long segment breaks the tie.
        let segments = [
            long.clone(),
            short.clone(),
            long.clone(),
            short,
            long.clone(),
        ];
        assert_eq!(backend.consensus_segments(&segments), long);
    }

    #[test]
    fn a_segment_that_does_not_fit_is_left_out_of_the_vote() {
        let backbone = b"ACGTTGCATGCAAGTCCGAT".to_vec();
        let mut deletion = backbone.clone();
        deletion.remove(10);
        let mut overlong = deletion.clone();
        overlong.extend(std::iter::repeat_n(b'A', 180));
        // The short segments fit in 600 cells per matrix, the 199-base one does not.
        let mut workspace = Workspace {
            aligner: Aligner::with_max_cells(600),
            ..Workspace::default()
        };
        let params = PoaParams::default();
        let segments = [overlong.clone(), deletion.clone(), backbone.clone()];
        // Had the skipped segment counted as a third vote against the deleted base, the
        // consensus would drop it; with two votes the backbone keeps it.
        assert_eq!(workspace.consensus(&segments, &params), backbone);
        // With room for every segment, the overlong one votes against the deleted base.
        let mut workspace = Workspace::default();
        assert_eq!(workspace.consensus(&segments, &params), deletion);
    }

    #[test]
    fn estimates_memory_per_thread() {
        assert_eq!(PoaBackend::MATRIX_BYTES_PER_THREAD, 48 << 20);
        assert_eq!(
            PoaBackend::peak_memory_per_thread(1000),
            (48 << 20) + 1000 * PoaBackend::BYTES_PER_WINDOW_BASE
        );
        assert_eq!(PoaBackend::peak_memory_per_thread(usize::MAX), usize::MAX);
    }

    #[test]
    fn rejects_invalid_params() {
        let params = PoaParams {
            gap_extend: 0,
            ..PoaParams::default()
        };
        assert_eq!(
            PoaBackend::new(params),
            Err(InvalidPoaParamsError::GapExtend)
        );
        let backend = PoaBackend::new(PoaParams::default()).expect("defaults are valid");
        assert_eq!(backend, PoaBackend::default());
        assert_eq!(backend.params(), &PoaParams::default());
    }
}
