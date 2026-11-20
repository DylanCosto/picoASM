//! Banded, affine-gap global alignment of a segment to a partial-order graph.
//!
//! Rows are graph nodes in topological order (row 0 is the virtual source) and columns are
//! segment positions `0..=m`. Three matrices are kept, as in Gotoh's algorithm: `h` (best
//! score), `e` (best score ending in a deletion, which consumes a node) and `f` (best score
//! ending in an insertion, which consumes a segment base).
//!
//! Each row is computed only on a band of columns around the best-scoring columns of its
//! predecessor rows, in the spirit of abPOA's adaptive band. If the traceback touches its row's
//! band edge, passes a cell with a candidate outside a predecessor row's band, or ends while a
//! sink predecessor's band misses the last column, the band may have cut off a better path, so
//! the alignment is repeated with a doubled band, until the band covers every column or the
//! matrices would exceed [`MAX_DP_CELLS`]. The predecessor and sink checks count only when an
//! upper bound on the cut-off cell's score could beat the traced path. These checks only see the traced path, so they are a
//! heuristic: an optimum that leaves the band away from the traced path can still be missed.

use crate::consensus::PoaParams;
use crate::consensus::graph::{AlignedPair, PoaGraph, SINK, bases_match};

/// The score of an unreachable cell. Stored scores are clamped to it, so a transition from any
/// stored score can subtract a penalty without overflowing.
const NEG_INF: i32 = i32::MIN / 4;

/// Scores at or below this are treated as unreachable when choosing where an alignment ends.
const UNREACHABLE: i32 = NEG_INF / 2;

/// The most banded matrix cells (per matrix) one alignment may use, 48 MiB for the three
/// matrices.
///
/// It bounds memory per thread on pathological windows (the matrices never reserve more than
/// this), and it bounds positive scores: with match scores of at most 255, no score can exceed
/// `MAX_DP_CELLS * 255`, well within `i32`. Raising it costs memory on huge windows only;
/// lowering it makes the band stop widening sooner.
pub(crate) const MAX_DP_CELLS: usize = 1 << 22;

/// The bytes of the three alignment matrices at [`MAX_DP_CELLS`] cells each.
pub(crate) const MAX_DP_BYTES: usize = 3 * MAX_DP_CELLS * size_of::<i32>();

/// Scores as signed values, ready for the recurrences.
#[derive(Clone, Copy, Debug)]
struct Scoring {
    matched: i32,
    mismatch: i32,
    open_extend: i32,
    extend: i32,
}

impl Scoring {
    fn new(params: &PoaParams) -> Self {
        Self {
            matched: i32::from(params.match_score),
            mismatch: i32::from(params.mismatch_penalty),
            open_extend: i32::from(params.gap_open) + i32::from(params.gap_extend),
            extend: i32::from(params.gap_extend),
        }
    }

    fn substitution(self, node_base: u8, seq_base: u8) -> i32 {
        if bases_match(node_base, seq_base) {
            self.matched
        } else {
            -self.mismatch
        }
    }
}

/// The matrix a traceback step is in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum State {
    Best,
    Deletion,
    Insertion,
}

/// The outcome of one traceback.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Trace {
    /// The path stays inside the band.
    Clean,
    /// The path touches a band edge; a wider band might find a better one.
    TouchesBand,
    /// No path reaches the end, or the matrices are inconsistent.
    Failed,
}

/// Reusable dynamic-programming buffers for aligning segments to a [`PoaGraph`].
#[derive(Debug)]
pub(crate) struct Aligner {
    /// The most cells per matrix; [`MAX_DP_CELLS`] except in tests.
    max_cells: usize,
    h: Vec<i32>,
    e: Vec<i32>,
    f: Vec<i32>,
    /// First column of each row's band.
    lo: Vec<usize>,
    /// Last column (inclusive) of each row's band.
    hi: Vec<usize>,
    /// Offset of each row's band in the flat matrices.
    offset: Vec<usize>,
    /// Column of each row's best `h` score, which centres the bands of its successors.
    best_col: Vec<usize>,
    /// The most nodes on any path from the source to each row's node (0 for the source).
    depth: Vec<usize>,
    path: Vec<AlignedPair>,
    best_path: Vec<AlignedPair>,
}

impl Default for Aligner {
    fn default() -> Self {
        Self::with_max_cells(MAX_DP_CELLS)
    }
}

impl Aligner {
    /// Creates an aligner whose matrices hold at most `max_cells` cells each.
    pub(crate) fn with_max_cells(max_cells: usize) -> Self {
        Self {
            max_cells,
            h: Vec::new(),
            e: Vec::new(),
            f: Vec::new(),
            lo: Vec::new(),
            hi: Vec::new(),
            offset: Vec::new(),
            best_col: Vec::new(),
            depth: Vec::new(),
            path: Vec::new(),
            best_path: Vec::new(),
        }
    }

    /// Globally aligns `seq` to `graph` (from source to sink), widening the band as needed.
    ///
    /// `graph` must be topologically sorted. If a wider band does not fit in the matrices
    /// ([`MAX_DP_CELLS`]), the best alignment from a narrower band is kept even though it touched
    /// a band edge. Returns `None` if no alignment fits; the caller should then leave the
    /// segment out.
    pub(crate) fn align(
        &mut self,
        graph: &PoaGraph,
        seq: &[u8],
        params: &PoaParams,
    ) -> Option<&[AlignedPair]> {
        let scoring = Scoring::new(params);
        let mut width = params.band_width.max(1);
        let mut found = false;
        loop {
            if !self.fill(graph, seq, width, scoring) {
                break;
            }
            match self.trace(graph, seq, scoring) {
                Trace::Clean => {
                    std::mem::swap(&mut self.path, &mut self.best_path);
                    found = true;
                    break;
                }
                Trace::TouchesBand => {
                    std::mem::swap(&mut self.path, &mut self.best_path);
                    found = true;
                }
                Trace::Failed => {}
            }
            if width > seq.len() {
                // The band already covered every column.
                break;
            }
            width = width.saturating_mul(2);
        }
        found.then_some(self.best_path.as_slice())
    }

    /// Fills the banded matrices. Returns `false` if they would exceed the cell limit.
    fn fill(&mut self, graph: &PoaGraph, seq: &[u8], width: usize, scoring: Scoring) -> bool {
        let m = seq.len();
        let rows = graph.rows();
        let n_rows = rows.len();
        for buf in [
            &mut self.lo,
            &mut self.hi,
            &mut self.offset,
            &mut self.best_col,
            &mut self.depth,
        ] {
            buf.clear();
            buf.resize(n_rows, 0);
        }
        self.h.clear();
        self.e.clear();
        self.f.clear();

        // Row 0, the source: only leading insertions.
        let hi0 = m.min(width);
        if !self.reserve_cells(hi0 + 1) {
            return false;
        }
        self.hi[0] = hi0;
        self.h.push(0);
        self.e.push(NEG_INF);
        self.f.push(NEG_INF);
        let mut gap = NEG_INF;
        for j in 1..=hi0 {
            gap = if j == 1 {
                -scoring.open_extend
            } else {
                (gap - scoring.extend).max(NEG_INF)
            };
            self.h.push(gap);
            self.e.push(NEG_INF);
            self.f.push(gap);
        }

        for (r, &node) in rows.iter().enumerate().skip(1) {
            let preds = graph.preds(node);
            let (mut min_col, mut max_col, mut depth) = (usize::MAX, 0, 0);
            for edge in preds {
                let p = graph.rank(edge.from);
                min_col = min_col.min(self.best_col[p]);
                max_col = max_col.max(self.best_col[p]);
                depth = depth.max(self.depth[p]);
            }
            self.depth[r] = depth + 1;
            debug_assert!(
                !preds.is_empty(),
                "only the source, row 0, has no predecessors"
            );
            let lo = min_col.saturating_add(1).saturating_sub(width).min(m);
            let hi = max_col.saturating_add(1).saturating_add(width).min(m);
            let len = hi - lo + 1;
            let offset = self.h.len();
            if !self.reserve_cells(offset + len) {
                return false;
            }
            self.lo[r] = lo;
            self.hi[r] = hi;
            self.offset[r] = offset;
            self.h.resize(offset + len, NEG_INF);
            self.e.resize(offset + len, NEG_INF);
            self.f.resize(offset + len, NEG_INF);

            let base = graph.base(node);
            let (done_h, row_h) = self.h.split_at_mut(offset);
            let (done_e, row_e) = self.e.split_at_mut(offset);
            for edge in preds {
                let p = graph.rank(edge.from);
                let (plo, phi, poff) = (self.lo[p], self.hi[p], self.offset[p]);
                // Deletion of this node: same column in the predecessor row.
                for j in lo.max(plo)..=hi.min(phi) {
                    let k = poff + j - plo;
                    let cand = (done_h[k] - scoring.open_extend).max(done_e[k] - scoring.extend);
                    let cell = &mut row_e[j - lo];
                    *cell = (*cell).max(cand);
                }
                // Match or mismatch: previous column in the predecessor row.
                for j in lo.max(plo + 1)..=hi.min(phi + 1) {
                    let cand = done_h[poff + j - 1 - plo] + scoring.substitution(base, seq[j - 1]);
                    let cell = &mut row_h[j - lo];
                    *cell = (*cell).max(cand);
                }
            }
            let row_f = &mut self.f[offset..];
            self.best_col[r] = lo + finish_row(row_h, row_e, row_f, scoring);
        }
        true
    }

    /// Makes room for `needed` cells in each matrix, or returns `false` if that exceeds the
    /// cell limit.
    ///
    /// Capacity grows geometrically but never past the limit, so the matrices' memory stays
    /// within it across windows.
    fn reserve_cells(&mut self, needed: usize) -> bool {
        if needed > self.max_cells {
            return false;
        }
        for matrix in [&mut self.h, &mut self.e, &mut self.f] {
            if needed > matrix.capacity() {
                let target = needed
                    .max(matrix.capacity().saturating_mul(2))
                    .min(self.max_cells);
                matrix.reserve_exact(target - matrix.len());
            }
        }
        true
    }

    /// Traces the best path back from the sink into `self.path`, in segment order.
    fn trace(&mut self, graph: &PoaGraph, seq: &[u8], scoring: Scoring) -> Trace {
        let m = seq.len();
        let mut end: Option<(i32, usize)> = None;
        for edge in graph.preds(SINK) {
            let r = graph.rank(edge.from);
            if self.hi[r] == m {
                let score = self.h[self.offset[r] + m - self.lo[r]];
                if end.is_none_or(|(best, _)| score > best) {
                    end = Some((score, r));
                }
            }
        }
        let Some((score, mut r)) = end else {
            return Trace::Failed;
        };
        if score <= UNREACHABLE {
            return Trace::Failed;
        }
        // A sink predecessor whose band stops short of the last column might hold a better end,
        // unless even its score bound cannot beat the end found.
        let mut touched = graph.preds(SINK).iter().any(|edge| {
            let p = graph.rank(edge.from);
            self.hi[p] < m && self.score_bound(p, m, scoring) > i64::from(score)
        });

        self.path.clear();
        let mut j = m;
        let mut state = State::Best;
        loop {
            let (lo, hi) = (self.lo[r], self.hi[r]);
            touched |= (j == lo && lo > 0) || (j == hi && hi < m);
            if r == 0 {
                // The source row: whatever is left of the segment is a leading insertion.
                self.path
                    .extend((0..j).rev().map(|pos| AlignedPair::Insertion { pos }));
                break;
            }
            let node = graph.rows()[r];
            let idx = self.offset[r] + j - lo;
            match state {
                State::Best => {
                    let score = self.h[idx];
                    // A predecessor band that misses column `j - 1` cut off a diagonal
                    // candidate worth at most its bound plus a match.
                    touched |= j > 0
                        && self.clipped(
                            graph,
                            node,
                            j - 1,
                            i64::from(score) - i64::from(scoring.matched),
                            scoring,
                        );
                    let diagonal = if j > 0 {
                        let sub = scoring.substitution(graph.base(node), seq[j - 1]);
                        graph.preds(node).iter().find_map(|edge| {
                            let p = graph.rank(edge.from);
                            let (plo, phi) = (self.lo[p], self.hi[p]);
                            let inside = (plo..=phi).contains(&(j - 1));
                            (inside && self.h[self.offset[p] + j - 1 - plo] + sub == score)
                                .then_some(p)
                        })
                    } else {
                        None
                    };
                    if let Some(p) = diagonal {
                        self.path.push(AlignedPair::Aligned { node, pos: j - 1 });
                        r = p;
                        j -= 1;
                    } else if score == self.e[idx] {
                        state = State::Deletion;
                    } else if score == self.f[idx] {
                        state = State::Insertion;
                    } else {
                        return Trace::Failed;
                    }
                }
                State::Deletion => {
                    let score = self.e[idx];
                    // A predecessor band that misses column `j` cut off a deletion candidate
                    // worth at most its bound minus a gap extension.
                    touched |= self.clipped(
                        graph,
                        node,
                        j,
                        i64::from(score) + i64::from(scoring.extend),
                        scoring,
                    );
                    let step = graph.preds(node).iter().find_map(|edge| {
                        let p = graph.rank(edge.from);
                        let (plo, phi) = (self.lo[p], self.hi[p]);
                        if !(plo..=phi).contains(&j) {
                            return None;
                        }
                        let k = self.offset[p] + j - plo;
                        if self.h[k] - scoring.open_extend == score {
                            Some((p, State::Best))
                        } else if self.e[k] - scoring.extend == score {
                            Some((p, State::Deletion))
                        } else {
                            None
                        }
                    });
                    let Some((p, next)) = step else {
                        return Trace::Failed;
                    };
                    self.path.push(AlignedPair::Deletion { node });
                    r = p;
                    state = next;
                }
                State::Insertion => {
                    if j == lo {
                        return Trace::Failed;
                    }
                    let score = self.f[idx];
                    self.path.push(AlignedPair::Insertion { pos: j - 1 });
                    state = if self.h[idx - 1] - scoring.open_extend == score {
                        State::Best
                    } else if self.f[idx - 1] - scoring.extend == score {
                        State::Insertion
                    } else {
                        return Trace::Failed;
                    };
                    j -= 1;
                }
            }
        }
        self.path.reverse();
        if touched {
            Trace::TouchesBand
        } else {
            Trace::Clean
        }
    }

    /// Returns whether some predecessor row of `node` has a band that misses column `j` while
    /// its [score bound](Self::score_bound) there exceeds `beats`: a candidate the band cut off
    /// that might have won.
    fn clipped(
        &self,
        graph: &PoaGraph,
        node: usize,
        j: usize,
        beats: i64,
        scoring: Scoring,
    ) -> bool {
        graph.preds(node).iter().any(|edge| {
            let p = graph.rank(edge.from);
            !(self.lo[p]..=self.hi[p]).contains(&j) && self.score_bound(p, j, scoring) > beats
        })
    }

    /// Returns an upper bound on the score of any alignment of the first `col` segment bases
    /// that ends on row `r`, banded or not.
    ///
    /// A path to the row has at most `depth[r]` nodes, so at most that many bases match and
    /// the rest of the `col` bases are insertions, which cost at least one gap.
    fn score_bound(&self, r: usize, col: usize, scoring: Scoring) -> i64 {
        let to_i64 = |n: usize| i64::try_from(n).unwrap_or(i64::MAX);
        let depth = self.depth[r];
        let matched = to_i64(depth.min(col)).saturating_mul(i64::from(scoring.matched));
        let inserted = col.saturating_sub(depth);
        if inserted == 0 {
            return matched;
        }
        let gap_open = i64::from(scoring.open_extend - scoring.extend);
        let gap = to_i64(inserted)
            .saturating_mul(i64::from(scoring.extend))
            .saturating_add(gap_open);
        matched.saturating_sub(gap)
    }
}

/// Completes one row of the matrices: `h` holds the best match or mismatch score of each cell
/// and `e` its best deletion score (both unclamped); fills in `f` (insertions, along the row),
/// clamps all three to [`NEG_INF`], leaves the best of them in `h`, and returns the first index
/// of the row's best `h`.
///
/// An insertion run continues from the cell before as
/// `f[j] = max(h[j - 1] - open_extend, f[j - 1] - extend)`, where `h[j - 1]` is itself at least
/// `f[j - 1]`; since `open_extend >= extend`, that is
/// `max(g[j - 1] - open_extend, f[j - 1] - extend)` with `g` the best of match and deletion
/// alone. So only the short `f` chain runs cell by cell, and the rest runs on whole slices.
fn finish_row(h: &mut [i32], e: &mut [i32], f: &mut [i32], scoring: Scoring) -> usize {
    debug_assert!(
        scoring.open_extend >= scoring.extend,
        "a gap opening is never cheaper than extending one"
    );
    debug_assert!(h.len() == e.len() && h.len() == f.len(), "one row each");
    // g: the best of a match, a mismatch and a deletion, clamped.
    for (h, e) in h.iter_mut().zip(e.iter_mut()) {
        *e = (*e).max(NEG_INF);
        *h = (*h).max(*e);
    }
    // Insertions opened after each cell, then carried along the row.
    f[0] = NEG_INF;
    for (f, &g) in f[1..].iter_mut().zip(h.iter()) {
        *f = (g - scoring.open_extend).max(NEG_INF);
    }
    for j in 1..f.len() {
        f[j] = f[j].max(f[j - 1] - scoring.extend);
    }
    for (h, &f) in h.iter_mut().zip(f.iter()) {
        *h = (*h).max(f);
    }
    let best = h.iter().copied().max().unwrap_or(NEG_INF);
    h.iter().position(|&s| s == best).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The row recurrences cell by cell, as [`finish_row`] must compute them.
    fn finish_row_by_cell(h: &mut [i32], e: &mut [i32], f: &mut [i32], scoring: Scoring) -> usize {
        let mut best = (NEG_INF, 0);
        for idx in 0..h.len() {
            let ins = if idx == 0 {
                NEG_INF
            } else {
                (h[idx - 1] - scoring.open_extend).max(f[idx - 1] - scoring.extend)
            };
            f[idx] = ins.max(NEG_INF);
            e[idx] = e[idx].max(NEG_INF);
            let score = h[idx].max(e[idx]).max(f[idx]).max(NEG_INF);
            h[idx] = score;
            if score > best.0 {
                best = (score, idx);
            }
        }
        best.1
    }

    #[test]
    fn a_row_finishes_as_cell_by_cell() {
        let mut state = 7_u64;
        let mut next = move |n: u64| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % n
        };
        for case in 0..2_000 {
            let scoring = Scoring::new(
                &PoaParams::new()
                    .with_gap_open(u8::try_from(next(8)).expect("small"))
                    .with_gap_extend(u8::try_from(1 + next(6)).expect("small")),
            );
            let len = usize::try_from(1 + next(40)).expect("small");
            // Scores near each other, some unreachable and some below NEG_INF, as the
            // match and deletion candidates leave them.
            let mut cell = || match next(6) {
                0 => NEG_INF,
                1 => NEG_INF - i32::try_from(next(10)).expect("small"),
                _ => i32::try_from(next(60)).expect("small") - 30,
            };
            let h: Vec<i32> = (0..len).map(|_| cell()).collect();
            let e: Vec<i32> = (0..len).map(|_| cell()).collect();
            let (mut h1, mut e1, mut f1) = (h.clone(), e.clone(), vec![0; len]);
            let (mut h2, mut e2, mut f2) = (h, e, vec![0; len]);
            let best1 = finish_row(&mut h1, &mut e1, &mut f1, scoring);
            let best2 = finish_row_by_cell(&mut h2, &mut e2, &mut f2, scoring);
            assert_eq!(
                (best1, &h1, &e1, &f1),
                (best2, &h2, &e2, &f2),
                "case {case}"
            );
        }
    }

    fn aligned(graph_seq: &[u8], seq: &[u8], params: &PoaParams) -> Vec<AlignedPair> {
        let mut graph = PoaGraph::default();
        graph.reset(graph_seq);
        assert!(graph.topo_sort());
        let mut aligner = Aligner::default();
        aligner
            .align(&graph, seq, params)
            .expect("a toy alignment fits in the matrices")
            .to_vec()
    }

    /// Counts (matches, mismatches, deletions, insertions) in an alignment to a linear graph.
    fn summary(graph_seq: &[u8], seq: &[u8], path: &[AlignedPair]) -> (usize, usize, usize, usize) {
        let mut counts = (0, 0, 0, 0);
        for pair in path {
            match *pair {
                // Linear graph: node ids are 2.. in sequence order.
                AlignedPair::Aligned { node, pos } if graph_seq[node - 2] == seq[pos] => {
                    counts.0 += 1;
                }
                AlignedPair::Aligned { .. } => counts.1 += 1,
                AlignedPair::Deletion { .. } => counts.2 += 1,
                AlignedPair::Insertion { .. } => counts.3 += 1,
            }
        }
        counts
    }

    /// A tiny deterministic generator (splitmix64).
    struct SplitMix(u64);

    impl SplitMix {
        fn below(&mut self, n: u64) -> usize {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut x = self.0;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            usize::try_from((x ^ (x >> 31)) % n).expect("a small bound fits in usize")
        }

        fn bases(&mut self, len: usize) -> Vec<u8> {
            (0..len).map(|_| b"ACGT"[self.below(4)]).collect()
        }

        /// Copies `truth` with a few substitutions and indels of up to `max_indel` bases.
        fn mutate(&mut self, truth: &[u8], max_indel: u64) -> Vec<u8> {
            let mut out = Vec::new();
            let mut i = 0;
            while i < truth.len() {
                match self.below(40) {
                    0 => out.push(b"ACGT"[self.below(4)]),
                    1 => {
                        let len = 1 + self.below(max_indel);
                        out.extend(self.bases(len));
                    }
                    2 => i += self.below(max_indel),
                    _ => {}
                }
                if let Some(&base) = truth.get(i) {
                    out.push(base);
                }
                i += 1;
            }
            out
        }
    }

    /// Recomputes the affine-gap score of an alignment path.
    fn path_score(graph: &PoaGraph, seq: &[u8], path: &[AlignedPair], params: &PoaParams) -> i32 {
        let scoring = Scoring::new(params);
        let mut score = 0;
        let mut prev: Option<State> = None;
        for &pair in path {
            let (state, step) = match pair {
                AlignedPair::Aligned { node, pos } => (
                    State::Best,
                    scoring.substitution(graph.base(node), seq[pos]),
                ),
                AlignedPair::Deletion { .. } => (State::Deletion, 0),
                AlignedPair::Insertion { .. } => (State::Insertion, 0),
            };
            let gap = match state {
                State::Best => 0,
                _ if prev == Some(state) => -scoring.extend,
                _ => -scoring.open_extend,
            };
            score += step + gap;
            prev = Some(state);
        }
        score
    }

    /// Builds a graph from `segments`, aligning each with an unbanded alignment.
    fn build_graph(segments: &[Vec<u8>]) -> PoaGraph {
        let wide = PoaParams {
            band_width: 1 << 20,
            ..PoaParams::default()
        };
        let mut graph = PoaGraph::default();
        graph.reset(&segments[0]);
        let mut aligner = Aligner::default();
        for segment in &segments[1..] {
            assert!(graph.topo_sort());
            let path = aligner
                .align(&graph, segment, &wide)
                .expect("a toy alignment fits")
                .to_vec();
            graph.add_alignment(segment, &path);
        }
        assert!(graph.topo_sort());
        graph
    }

    /// Asserts that a band of width 1 (widened as needed) finds an optimal alignment.
    fn assert_band_is_optimal(graph: &PoaGraph, seq: &[u8], context: &str) {
        let narrow = PoaParams {
            band_width: 1,
            ..PoaParams::default()
        };
        let wide = PoaParams {
            band_width: 1 << 20,
            ..PoaParams::default()
        };
        let mut aligner = Aligner::default();
        let best = aligner.align(graph, seq, &wide).expect("fits").to_vec();
        let banded = aligner.align(graph, seq, &narrow).expect("fits").to_vec();
        assert_eq!(
            path_score(graph, seq, &banded, &narrow),
            path_score(graph, seq, &best, &wide),
            "{context}: the banded alignment is not optimal"
        );
    }

    /// Aligns a mutated copy of a random truth to a graph of four other copies, with indels of
    /// up to 8 bases, and checks that a narrow band finds the optimum.
    fn assert_branchy_case_is_optimal(seed: u64) {
        let mut rng = SplitMix(seed);
        let truth = rng.bases(60);
        let segments: Vec<Vec<u8>> = (0..4).map(|_| rng.mutate(&truth, 8)).collect();
        let graph = build_graph(&segments);
        let seq = rng.mutate(&truth, 8);
        assert_band_is_optimal(&graph, &seq, &format!("seed {seed}"));
    }

    #[test]
    fn widens_when_a_predecessor_band_clips_the_path() {
        // On this seed the traced path passes cells whose diagonal or deletion candidate lay
        // outside a predecessor row's band, deep inside the cell's own row band.
        assert_branchy_case_is_optimal(48);
    }

    #[test]
    fn widens_when_the_best_end_lies_outside_a_band() {
        // On this seed the best alignment ends on a sink predecessor whose narrow band never
        // reaches the last column.
        assert_branchy_case_is_optimal(441);
    }

    #[test]
    fn truncated_segments_do_not_force_a_wider_band() {
        // A prefix-truncated segment adds an edge to the sink from mid-graph, and a
        // suffix-truncated one an edge from the source to mid-graph. Neither can hold the best
        // alignment of a full-length segment, so they must not trigger widening.
        let mut rng = SplitMix(9);
        let truth = rng.bases(300);
        let graph = build_graph(&[truth.clone(), truth[..150].to_vec(), truth[150..].to_vec()]);
        let seq = rng.mutate(&truth, 2);
        let scoring = Scoring::new(&PoaParams::default());
        let mut aligner = Aligner::default();
        assert!(aligner.fill(&graph, &seq, PoaParams::DEFAULT_BAND_WIDTH, scoring));
        assert_eq!(aligner.trace(&graph, &seq, scoring), Trace::Clean);
    }

    #[test]
    fn gives_up_when_nothing_fits() {
        let graph_seq = b"ACGTTGCATGCAAGTCCGAT";
        let graph = build_graph(&[graph_seq.to_vec()]);
        let mut aligner = Aligner::with_max_cells(10);
        assert!(
            aligner
                .align(&graph, graph_seq, &PoaParams::default())
                .is_none()
        );
        assert!(
            aligner.h.capacity() <= 10,
            "capacity {}",
            aligner.h.capacity()
        );
    }

    #[test]
    fn keeps_a_narrow_alignment_when_a_wider_band_does_not_fit() {
        let graph_seq = b"ACGTTGCATGCAAGTCCGAT";
        let seq = b"ACGTTGCATGGGCAAGTCCGAT";
        let graph = build_graph(&[graph_seq.to_vec()]);
        let scoring = Scoring::new(&PoaParams::default());
        // A 2-base insertion touches the edge of a band of width 2. Size the matrices to
        // exactly what that band needs, so the doubled band does not fit.
        let width = 2;
        let mut probe = Aligner::default();
        assert!(probe.fill(&graph, seq, width, scoring));
        assert_eq!(probe.trace(&graph, seq, scoring), Trace::TouchesBand);
        let narrow_path = probe.path.clone();
        let cells = probe.h.len();

        let params = PoaParams {
            band_width: width,
            ..PoaParams::default()
        };
        let mut aligner = Aligner::with_max_cells(cells);
        let path = aligner
            .align(&graph, seq, &params)
            .expect("the width-2 alignment fits");
        assert_eq!(path, narrow_path.as_slice());
        for matrix in [&aligner.h, &aligner.e, &aligner.f] {
            assert!(matrix.capacity() <= cells, "capacity {}", matrix.capacity());
        }
    }

    #[test]
    fn fails_on_a_narrow_band_then_recovers() {
        // A segment far longer than the graph: no narrow band reaches the last column.
        let graph = build_graph(&[b"ACGT".to_vec()]);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let scoring = Scoring::new(&PoaParams::default());
        let mut aligner = Aligner::default();
        assert!(aligner.fill(&graph, seq, 1, scoring));
        assert_eq!(aligner.trace(&graph, seq, scoring), Trace::Failed);
        let params = PoaParams {
            band_width: 1,
            ..PoaParams::default()
        };
        let path = aligner
            .align(&graph, seq, &params)
            .expect("a wider band fits");
        assert_eq!(summary(b"ACGT", seq, path), (4, 0, 0, 28));
    }

    #[test]
    fn aligns_with_a_gap() {
        let params = PoaParams::default();
        let path = aligned(b"ACGTACGTAC", b"ACGTCGTAC", &params);
        assert_eq!(summary(b"ACGTACGTAC", b"ACGTCGTAC", &path), (9, 0, 1, 0));
        let path = aligned(b"ACGTACGTAC", b"ACGTTACGTAC", &params);
        assert_eq!(summary(b"ACGTACGTAC", b"ACGTTACGTAC", &path), (10, 0, 0, 1));
    }

    #[test]
    fn widens_a_narrow_band_for_a_long_insertion() {
        let params = PoaParams {
            band_width: 1,
            ..PoaParams::default()
        };
        let graph_seq = b"ACGTTGCATGCAAGTCCGAT";
        let seq = b"ACGTTGCATGGGGGGGGGGGGCAAGTCCGAT";
        let path = aligned(graph_seq, seq, &params);
        assert_eq!(summary(graph_seq, seq, &path), (20, 0, 0, 11));
    }

    #[test]
    fn aligns_a_truncated_segment_with_end_gaps() {
        let params = PoaParams::default();
        let graph_seq = b"ACGTTGCATGCAAGTCCGATCGATTACG";
        let path = aligned(graph_seq, &graph_seq[10..], &params);
        assert_eq!(summary(graph_seq, &graph_seq[10..], &path), (18, 0, 10, 0));
    }
}
