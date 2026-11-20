//! Haplotype-aware read correction, hifiasm's rule: each read is corrected only from the reads
//! that agree with it where copies differ.
//!
//! For each read, the reads overlapping it are stacked on its bases. A column is a **site**
//! where at least [`CorrectionParams::min_alt`] stacked reads carry one base other than the
//! read's, at least as many (the read itself included) carry the read's own, and no run of four
//! of a base lies within three columns (runs are where HiFi reads err). Two copies of a repeat,
//! or two haplotypes, differ there. A stacked read with another base at a site comes from
//! another copy and is dropped from the stack. Every other column takes the majority of the
//! read and its remaining stack (a base, a deletion, or an insertion before it), where the stack
//! covers at least [`CorrectionParams::min_depth`] reads and the majority holds more than half
//! of them.
//!
//! # Limits
//!
//! - A read's whole stack is held at once, two bytes per column of each overlap: a 20 kb read
//!   with 3,000 overlaps (a young transposon) takes some 120 MB, on each thread.
//! - A read the overlaps miss (one with more errors than the overlap step allows) is left as
//!   it is.

use std::ops::Range;

use crate::seq::complement;
use rayon::prelude::*;

use crate::overlap::align::{Edit, EditOp};
use crate::overlap::overlap::Overlap;
use crate::overlap::reads::{ReadIndex, Reads};

/// Settings for [`correct`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CorrectionParams {
    /// A site needs at least this many reads with each of two bases.
    pub min_alt: u16,
    /// A column is corrected only where at least this many reads cover it.
    pub min_depth: u16,
}

impl Default for CorrectionParams {
    fn default() -> Self {
        Self {
            min_alt: 3,
            min_depth: 3,
        }
    }
}

/// The reads after one round of correction, and which overlaps join reads of one copy.
#[derive(Clone, Debug)]
pub struct Correction {
    /// The corrected reads, in the input order.
    pub reads: Reads,
    /// For each input overlap, whether neither read has the other's base at one of its sites.
    pub consistent: Vec<bool>,
    /// Sites found, over all reads.
    pub sites: u64,
    /// Bases changed, inserted or deleted, over all reads.
    pub changes: u64,
    /// For each read, whether its bases changed.
    pub changed: Vec<bool>,
}

/// Corrects every read in `reads` once from `overlaps` (found on these reads).
///
/// The result is the same for any number of threads; work runs on the caller's rayon pool.
///
/// # Panics
///
/// Panics if an overlap names a read not in `reads`.
#[must_use]
pub fn correct(reads: &Reads, overlaps: &[Overlap], params: &CorrectionParams) -> Correction {
    // For each read, the overlaps it is in, and whether it is their target.
    let mut by_read: Vec<Vec<(u32, Role)>> = vec![Vec::new(); reads.len()];
    for (i, o) in overlaps.iter().enumerate() {
        let i = u32::try_from(i).expect("fewer than u32::MAX overlaps");
        by_read[o.target.index()].push((i, Role::Target));
        by_read[o.query.index()].push((i, Role::Query));
    }
    let mut consistent = vec![true; overlaps.len()];
    let mut changed = Vec::with_capacity(reads.len());
    let (mut sites, mut changes) = (0, 0);
    // Corrections change a read's length by a few bases; the room spare saves regrowing.
    let mut corrected = Reads::with_capacity(reads.len(), reads.total_bases() * 101 / 100);
    let indices: Vec<ReadIndex> = reads.indices().collect();
    // A chunk of reads at a time, so the corrected reads are held twice only a chunk at once.
    for chunk in indices.chunks(CHUNK) {
        let per_read: Vec<ReadResult> = chunk
            .par_iter()
            .map(|&read| correct_one(reads, overlaps, &by_read[read.index()], read, params))
            .collect();
        for (&read, result) in chunk.iter().zip(per_read) {
            changed.push(result.bases != reads.get(read));
            for i in result.disagreeing {
                consistent[usize::try_from(i).expect("overlap indices fit")] = false;
            }
            sites += result.sites;
            changes += result.changes;
            corrected.push(&result.bases);
        }
    }
    Correction {
        reads: corrected,
        consistent,
        sites,
        changes,
        changed,
    }
}

/// Reads corrected in parallel at once: enough to keep every thread busy.
const CHUNK: usize = 4_096;

/// Which read of an overlap is being corrected.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Role {
    Target,
    Query,
}

/// Another read laid on the read being corrected, as runs of alignment columns of one op in
/// the read's order, with each indel moved as far left as it goes (so that reads aligned from
/// either strand, or either end, vote for the same column).
struct View {
    /// The read position of the first column.
    on_read: usize,
    /// The columns, as maximal runs of one op.
    runs: Vec<(EditOp, usize)>,
    /// The other read's base at a match, mismatch or insertion; the read's at a deletion.
    bases: Vec<u8>,
    /// The read positions the columns cover.
    covered: usize,
}

impl View {
    fn new(reads: &Reads, own: &[u8], o: &Overlap, role: Role) -> Self {
        let (target, query) = (reads.get(o.target), reads.get(o.query));
        let reverse = o.strand.is_reverse();
        // The other read as aligned to this one, read in place: reverse-complemented on the
        // reverse strand. Copying it for each overlap took a tenth of the correction's time.
        let (other, on_read, other_start) = match role {
            Role::Target => {
                let start = if reverse {
                    query.len() - pos(o.query_end)
                } else {
                    pos(o.query_start)
                };
                (query, pos(o.target_start), start)
            }
            Role::Query => {
                let start = if reverse {
                    target.len() - pos(o.target_end)
                } else {
                    pos(o.target_start)
                };
                (target, pos(o.query_start), start)
            }
        };
        // The other read's bases `from..from + len` as aligned, appended to `bases`.
        let push_other = |bases: &mut Vec<u8>, from: usize, len: usize| {
            if reverse {
                let end = other.len() - from;
                bases.extend(other[end - len..end].iter().rev().map(|&b| complement(b)));
            } else {
                bases.extend_from_slice(&other[from..from + len]);
            }
        };
        // Laid on the query, the target's alignment swaps insertions and deletions, and on the
        // reverse strand runs backwards.
        let edit_op = |e: &Edit| match (role, e.op()) {
            (Role::Query, EditOp::Insertion) => EditOp::Deletion,
            (Role::Query, EditOp::Deletion) => EditOp::Insertion,
            (_, op) => op,
        };
        let edits: Box<dyn Iterator<Item = &Edit>> = match role {
            Role::Query if reverse => Box::new(o.edits.iter().rev()),
            _ => Box::new(o.edits.iter()),
        };
        let columns = o.edits.iter().map(|e| pos(e.len())).sum();
        let (mut ops, mut bases) = (Vec::with_capacity(columns), Vec::with_capacity(columns));
        // The indels' columns, for the left alignment.
        let mut gaps: Vec<Range<usize>> = Vec::new();
        let (mut i, mut j) = (on_read, other_start);
        for e in edits {
            let op = edit_op(e);
            let len = pos(e.len());
            if matches!(op, EditOp::Insertion | EditOp::Deletion) {
                gaps.push(ops.len()..ops.len() + len);
            }
            ops.extend(std::iter::repeat_n(op, len));
            match op {
                EditOp::Match | EditOp::Mismatch => {
                    push_other(&mut bases, j, len);
                    i += len;
                    j += len;
                }
                EditOp::Deletion => {
                    bases.extend_from_slice(&own[i..i + len]);
                    i += len;
                }
                EditOp::Insertion => {
                    push_other(&mut bases, j, len);
                    j += len;
                }
            }
        }
        left_align(&mut ops, &bases, &gaps);
        Self::from_columns(on_read, &ops, bases)
    }

    /// The view of the columns `ops`, whose bases are `bases`, from read position `on_read`.
    fn from_columns(on_read: usize, ops: &[EditOp], bases: Vec<u8>) -> Self {
        let mut runs: Vec<(EditOp, usize)> = Vec::new();
        let mut c = 0;
        while c < ops.len() {
            let op = ops[c];
            let len = ops[c..].iter().take_while(|&&o| o == op).count();
            runs.push((op, len));
            c += len;
        }
        let covered = runs
            .iter()
            .filter(|&&(op, _)| op != EditOp::Insertion)
            .map(|&(_, len)| len)
            .sum();
        Self {
            on_read,
            runs,
            bases,
            covered,
        }
    }

    /// Calls `matched` for each run of the read's positions where this view has a match (its
    /// base is the read's own there), `column` for each other position it covers, with the
    /// view's base there (`None` for a deletion), and `insertion` for each run of bases it has
    /// between two of the read's columns it covers (before the given one).
    fn walk(
        &self,
        mut matched: impl FnMut(Range<usize>),
        mut column: impl FnMut(usize, Option<u8>),
        mut insertion: impl FnMut(usize, &[u8]),
    ) {
        let (mut i, mut c) = (self.on_read, 0);
        for &(op, len) in &self.runs {
            match op {
                EditOp::Match => {
                    matched(i..i + len);
                    i += len;
                }
                EditOp::Mismatch => {
                    for (k, &b) in self.bases[c..c + len].iter().enumerate() {
                        column(i + k, Some(b));
                    }
                    i += len;
                }
                EditOp::Deletion => {
                    for k in 0..len {
                        column(i + k, None);
                    }
                    i += len;
                }
                EditOp::Insertion => {
                    // Only between two of the read's columns the view covers: a run at either
                    // end of the view lies outside it.
                    if c > 0 && c + len < self.bases.len() {
                        insertion(i, &self.bases[c..c + len]);
                    }
                }
            }
            c += len;
        }
    }

    /// Returns the read positions this view covers.
    fn span(&self) -> (usize, usize) {
        (self.on_read, self.on_read + self.covered)
    }
}

/// Moves every run of insertions or deletions left past the matches before it while the base
/// the run would give up equals the one it would take on. In a homopolymer or a short tandem
/// repeat an indel can sit at any of several columns; this puts it at the first.
///
/// `gaps` are the columns of the indels as aligned (in order): the runs are found there rather
/// than by a scan of every column. A run's move only changes columns before its end, and the
/// scan goes on past it, so it finds the indels after it where they were aligned.
fn left_align(ops: &mut [EditOp], bases: &[u8], gaps: &[Range<usize>]) {
    // Where the scan goes on: every column before it is done with.
    let mut c = 0;
    for gap_run in gaps {
        if gap_run.end <= c {
            continue;
        }
        let mut start = gap_run.start.max(c);
        let gap = ops[start];
        debug_assert!(
            matches!(gap, EditOp::Insertion | EditOp::Deletion),
            "an indel's column"
        );
        let mut end = start;
        while end < ops.len() && ops[end] == gap {
            end += 1;
        }
        // Columns keep their bases; only the ops move. Shifting the run [start, end) one
        // column left swaps the match before it to its last column, which holds the same base.
        while start > 0 && ops[start - 1] == EditOp::Match && bases[start - 1] == bases[end - 1] {
            ops[start - 1..end].rotate_left(1);
            start -= 1;
            end -= 1;
        }
        c = end + 1;
    }
}

/// One read's correction.
struct ReadResult {
    bases: Vec<u8>,
    /// The overlaps whose other read has another base at one of this read's sites.
    disagreeing: Vec<u32>,
    sites: u64,
    changes: u64,
}

/// The index of a base among A, C, G, T, or 4 for a deletion.
const DELETION: usize = 4;

fn class(base: Option<u8>) -> Option<usize> {
    match base {
        None => Some(DELETION),
        Some(b'A') => Some(0),
        Some(b'C') => Some(1),
        Some(b'G') => Some(2),
        Some(b'T') => Some(3),
        Some(_) => None,
    }
}

fn correct_one(
    reads: &Reads,
    overlaps: &[Overlap],
    mine: &[(u32, Role)],
    read: ReadIndex,
    params: &CorrectionParams,
) -> ReadResult {
    let own = reads.get(read);
    let views: Vec<View> = mine
        .iter()
        .map(|&(i, role)| View::new(reads, own, &overlaps[pos(i)], role))
        .collect();

    // Sites, from every stacked read; and each stacked read's insertions, for the majority.
    // u32 counts do not overflow; the u16 counts before them saturated only past 65,535 reads
    // stacked on one column.
    let mut counts = vec![[0_u32; 5]; own.len()];
    let mut all_insertions: Vec<(usize, usize, Vec<u8>)> = Vec::new();
    // Matches, which carry the read's own base, are counted as runs: +1 where one starts, -1
    // past its end.
    let mut matching = vec![0_i32; own.len() + 1];
    for (k, v) in views.iter().enumerate() {
        v.walk(
            |run| {
                matching[run.start] += 1;
                matching[run.end] -= 1;
            },
            |p, b| {
                if let Some(c) = class(b) {
                    counts[p][c] += 1;
                }
            },
            |p, bases| all_insertions.push((k, p, bases.to_vec())),
        );
    }
    let mut matched = 0_i32;
    for ((column, &b), &step) in counts.iter_mut().zip(own).zip(&matching) {
        matched += step;
        if let Some(c) = class(Some(b)) {
            column[c] += u32::try_from(matched).expect("runs end after they start");
        }
    }
    let is_site: Vec<bool> = (0..own.len())
        .map(|p| site(own, &counts[p], p, params.min_alt))
        .collect();
    let sites = is_site.iter().filter(|&&s| s).count() as u64;

    // Drop the stacked reads with another base at a site, and their bases from the counts.
    let mut disagreeing = Vec::new();
    let mut agrees = vec![true; views.len()];
    if sites > 0 {
        for ((v, &(i, _)), agree) in views.iter().zip(mine).zip(&mut agrees) {
            // A match has the read's own base: only the other columns can disagree.
            v.walk(
                |_| {},
                |p, b| {
                    if is_site[p] && b.is_some_and(|b| b != own[p] && class(Some(b)).is_some()) {
                        *agree = false;
                    }
                },
                |_, _| {},
            );
            if !*agree {
                disagreeing.push(i);
                // A rare path: the matches column by column, then the other columns.
                v.walk(
                    |run| {
                        for p in run {
                            if let Some(c) = class(Some(own[p])) {
                                counts[p][c] -= 1;
                            }
                        }
                    },
                    |_, _| {},
                    |_, _| {},
                );
                v.walk(
                    |_| {},
                    |p, b| {
                        if let Some(c) = class(b) {
                            counts[p][c] -= 1;
                        }
                    },
                    |_, _| {},
                );
            }
        }
    }

    // The majority of the read and the rest of its stack.
    // Reads (the read included) covering each column, and the gap before it (both neighbours).
    let mut depth = vec![1_u32; own.len()];
    let mut gap_depth = vec![1_u32; own.len()];
    for (v, _) in views.iter().zip(&agrees).filter(|(_, a)| **a) {
        let (from, to) = v.span();
        for d in &mut depth[from..to] {
            *d += 1;
        }
        for d in &mut gap_depth[(from + 1).min(to)..to] {
            *d += 1;
        }
    }
    let mut insertions: Vec<(usize, Vec<u8>)> = all_insertions
        .into_iter()
        .filter(|(k, ..)| agrees[*k])
        .map(|(_, p, run)| (p, run))
        .collect();
    for (p, &b) in own.iter().enumerate() {
        if let Some(c) = class(Some(b)) {
            counts[p][c] += 1;
        }
    }
    insertions.sort_unstable();

    let mut bases = Vec::with_capacity(own.len());
    let mut changes = 0;
    let mut next = 0;
    for p in 0..own.len() {
        // The most common insertion before column p, if most reads spanning the gap carry it.
        let from = next;
        while next < insertions.len() && insertions[next].0 == p {
            next += 1;
        }
        if p > 0
            && gap_depth[p] >= u32::from(params.min_depth)
            && let Some((run, n)) = most_common(&insertions[from..next])
            && 2 * n > gap_depth[p]
        {
            bases.extend_from_slice(run);
            changes += run.len() as u64;
        }
        let column = &counts[p];
        let best = (0..5)
            .max_by_key(|&c| (column[c], std::cmp::Reverse(c)))
            .expect("five classes");
        let keep_own = depth[p] < u32::from(params.min_depth) || 2 * column[best] <= depth[p];
        match (keep_own, best) {
            (true, _) => bases.push(own[p]),
            (false, DELETION) => changes += 1,
            (false, c) => {
                let base = b"ACGT"[c];
                changes += u64::from(base != own[p]);
                bases.push(base);
            }
        }
    }
    ReadResult {
        bases,
        disagreeing,
        sites,
        changes,
    }
}

/// Returns `true` if column `p` of `own` is a site: at least `min_alt` stacked reads carry one
/// base other than the read's, at least `min_alt` reads (the read included) carry the read's,
/// and no run of four of a base lies within three columns.
fn site(own: &[u8], counts: &[u32; 5], p: usize, min_alt: u16) -> bool {
    let min_alt = u32::from(min_alt);
    let Some(mine) = class(Some(own[p])).filter(|&c| c != DELETION) else {
        return false;
    };
    if counts[mine] + 1 < min_alt || !(0..4).any(|b| b != mine && counts[b] >= min_alt) {
        return false;
    }
    let around = &own[p.saturating_sub(3)..(p + 4).min(own.len())];
    !around.windows(4).any(|w| w.iter().all(|&b| b == w[0]))
}

/// Returns the most common run among `runs` (all at one position) and how many carry it.
fn most_common(runs: &[(usize, Vec<u8>)]) -> Option<(&[u8], u32)> {
    let mut best: Option<(&[u8], u32)> = None;
    let mut i = 0;
    while i < runs.len() {
        let mut j = i;
        while j < runs.len() && runs[j].1 == runs[i].1 {
            j += 1;
        }
        let n = u32::try_from(j - i).unwrap_or(u32::MAX);
        if best.is_none_or(|(_, b)| n > b) {
            best = Some((&runs[i].1, n));
        }
        i = j;
    }
    best
}

fn pos(x: u32) -> usize {
    usize::try_from(x).unwrap_or_else(|_| panic!("position {x} does not fit in usize"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::overlap::{OverlapParams, find_overlaps};

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

    fn other(b: u8) -> u8 {
        if b == b'A' { b'C' } else { b'A' }
    }

    fn params() -> OverlapParams {
        OverlapParams {
            density: 0.2,
            min_overlap: 500,
            ..OverlapParams::default()
        }
    }

    /// Six reads tiling `g` from `offset`, every 400 bases, each 3,000 long.
    fn tiling(g: &[u8], count: usize) -> Vec<Vec<u8>> {
        (0..count)
            .map(|i| g[i * 400..i * 400 + 3_000].to_vec())
            .collect()
    }

    /// The indels' columns of `ops`, one range per run of one op.
    fn gap_runs(ops: &[EditOp]) -> Vec<Range<usize>> {
        let mut gaps = Vec::new();
        let mut c = 0;
        while c < ops.len() {
            let op = ops[c];
            let len = ops[c..].iter().take_while(|&&o| o == op).count();
            if matches!(op, EditOp::Insertion | EditOp::Deletion) {
                gaps.push(c..c + len);
            }
            c += len;
        }
        gaps
    }

    /// [`left_align`] as it was, scanning every column.
    fn left_align_scanning(ops: &mut [EditOp], bases: &[u8]) {
        let mut c = 0;
        while c < ops.len() {
            let gap = ops[c];
            if !matches!(gap, EditOp::Insertion | EditOp::Deletion) {
                c += 1;
                continue;
            }
            let mut end = c;
            while end < ops.len() && ops[end] == gap {
                end += 1;
            }
            let mut start = c;
            while start > 0 && ops[start - 1] == EditOp::Match && bases[start - 1] == bases[end - 1]
            {
                ops[start - 1..end].rotate_left(1);
                start -= 1;
                end -= 1;
            }
            c = end + 1;
        }
    }

    #[test]
    fn indels_move_as_the_scan_of_every_column_moved_them() {
        use EditOp::{Deletion as D, Insertion as I, Match as M, Mismatch as X};
        let mut state = 21_u64;
        let mut next = move |n: u64| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            usize::try_from((state >> 33) % n).expect("small")
        };
        for case in 0..5_000 {
            // Runs of every op side by side over a two-letter alphabet, so that indels move far.
            let mut ops = Vec::new();
            while ops.len() < 60 {
                let op = [M, M, M, X, I, D][next(6)];
                ops.extend(std::iter::repeat_n(op, 1 + next(4)));
            }
            let bases: Vec<u8> = (0..ops.len()).map(|_| b"AC"[next(2)]).collect();
            let gaps = gap_runs(&ops);
            let mut by_scan = ops.clone();
            left_align_scanning(&mut by_scan, &bases);
            left_align(&mut ops, &bases, &gaps);
            assert_eq!(ops, by_scan, "case {case}");
        }
    }

    #[test]
    fn indels_move_to_the_first_column_of_a_run() {
        use EditOp::{Deletion as D, Insertion as I, Match as M};
        // A deletion of one A from the read's run AAA, aligned at its last A.
        let mut ops = [M, M, M, D, M];
        let gaps = gap_runs(&ops);
        left_align(&mut ops, b"CAAAG", &gaps);
        assert_eq!(ops, [M, D, M, M, M]);
        // An insertion of CA into the read's CACA, aligned after it.
        let mut ops = [M, M, M, M, M, I, I, M];
        let gaps = gap_runs(&ops);
        left_align(&mut ops, b"GCACACAT", &gaps);
        assert_eq!(ops, [M, I, I, M, M, M, M, M]);
        // A mismatch stops the move.
        let mut ops = [M, EditOp::Mismatch, M, D];
        let gaps = gap_runs(&ops);
        left_align(&mut ops, b"AAAA", &gaps);
        assert_eq!(ops, [M, EditOp::Mismatch, D, M]);
    }

    #[test]
    fn insertions_at_a_views_ends_are_not_voted() {
        use EditOp::{Insertion as I, Match as M};
        let view = View::from_columns(5, &[I, M, I, I, M, I], b"AACGTT".to_vec());
        let mut seen = Vec::new();
        view.walk(|_| {}, |_, _| {}, |p, run| seen.push((p, run.to_vec())));
        assert_eq!(seen, [(6, b"CG".to_vec())]);
        assert_eq!(view.span(), (5, 7));
    }

    /// What a view shows, column by column: each covered read position with its base (`None`
    /// for a deletion), and each insertion between covered columns.
    type Shown = (Vec<(usize, Option<u8>)>, Vec<(usize, Vec<u8>)>);

    /// [`View::new`]'s view, column by column, as the walk shows it (matches with the read's
    /// own base).
    fn shown(view: &View, own: &[u8]) -> Shown {
        let columns = std::cell::RefCell::new(Vec::new());
        let mut insertions = Vec::new();
        view.walk(
            |run| columns.borrow_mut().extend(run.map(|p| (p, Some(own[p])))),
            |p, b| columns.borrow_mut().push((p, b)),
            |p, run| insertions.push((p, run.to_vec())),
        );
        (columns.into_inner(), insertions)
    }

    /// The view of overlap `o` built and walked one column at a time, as it was before views
    /// were runs.
    fn shown_by_column(reads: &Reads, own: &[u8], o: &Overlap, role: Role) -> Shown {
        let (target, query) = (reads.get(o.target), reads.get(o.query));
        let reverse = o.strand.is_reverse();
        let (other, on_read, other_start) = match role {
            Role::Target => {
                let start = if reverse {
                    query.len() - pos(o.query_end)
                } else {
                    pos(o.query_start)
                };
                (query, pos(o.target_start), start)
            }
            Role::Query => {
                let start = if reverse {
                    target.len() - pos(o.target_end)
                } else {
                    pos(o.target_start)
                };
                (target, pos(o.query_start), start)
            }
        };
        let other_at = |j: usize| {
            if reverse {
                complement(other[other.len() - 1 - j])
            } else {
                other[j]
            }
        };
        let edits: Vec<Edit> = match role {
            Role::Query if reverse => o.edits.iter().rev().copied().collect(),
            _ => o.edits.clone(),
        };
        let (mut ops, mut bases) = (Vec::new(), Vec::new());
        let (mut i, mut j) = (on_read, other_start);
        for e in &edits {
            let op = match (role, e.op()) {
                (Role::Query, EditOp::Insertion) => EditOp::Deletion,
                (Role::Query, EditOp::Deletion) => EditOp::Insertion,
                (_, op) => op,
            };
            for _ in 0..e.len() {
                ops.push(op);
                match op {
                    EditOp::Match | EditOp::Mismatch => {
                        bases.push(other_at(j));
                        i += 1;
                        j += 1;
                    }
                    EditOp::Deletion => {
                        bases.push(own[i]);
                        i += 1;
                    }
                    EditOp::Insertion => {
                        bases.push(other_at(j));
                        j += 1;
                    }
                }
            }
        }
        left_align_scanning(&mut ops, &bases);
        let (mut columns, mut insertions) = (Vec::new(), Vec::new());
        let (mut i, mut c) = (on_read, 0);
        while c < ops.len() {
            match ops[c] {
                EditOp::Match | EditOp::Mismatch => {
                    columns.push((i, Some(bases[c])));
                    i += 1;
                    c += 1;
                }
                EditOp::Deletion => {
                    columns.push((i, None));
                    i += 1;
                    c += 1;
                }
                EditOp::Insertion => {
                    let from = c;
                    while c < ops.len() && ops[c] == EditOp::Insertion {
                        c += 1;
                    }
                    if from > 0 && c < ops.len() {
                        insertions.push((i, bases[from..c].to_vec()));
                    }
                }
            }
        }
        (columns, insertions)
    }

    #[test]
    fn views_show_what_column_by_column_views_did() {
        // Noisy reads of a genome with homopolymers and tandem repeats, on both strands, so
        // that indels move left and runs of every op meet.
        let mut g = genome(12_000, 11);
        for at in (500..11_000).step_by(700) {
            g[at..at + 6].fill(b'T');
            let unit = g[at + 20..at + 23].to_vec();
            for k in 1..5 {
                g[at + 20 + 3 * k..at + 23 + 3 * k].copy_from_slice(&unit);
            }
        }
        let mut state = 13_u64;
        let mut next = move |n: u64| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % n
        };
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for i in 0..22 {
            let start = i * 400;
            let mut read = Vec::new();
            for &b in &g[start..start + 3_000] {
                match next(150) {
                    0 => read.push(other(b)),
                    1 => {}
                    2 => read.extend_from_slice(&[b, b]),
                    _ => read.push(b),
                }
            }
            if i % 2 == 1 {
                read = crate::seq::reverse_complement(&read);
            }
            reads.push(read);
        }
        let set: Reads = reads.iter().map(Vec::as_slice).collect();
        let params = OverlapParams {
            max_divergence: 0.05,
            ..params()
        };
        let overlaps = find_overlaps(&set, &params).expect("valid");
        assert!(overlaps.len() > 40, "{} overlaps", overlaps.len());
        for o in &overlaps {
            for (role, read) in [(Role::Target, o.target), (Role::Query, o.query)] {
                let own = set.get(read);
                let view = View::new(&set, own, o, role);
                assert_eq!(
                    shown(&view, own),
                    shown_by_column(&set, own, o, role),
                    "{role:?} {} {}",
                    o.query.0,
                    o.target.0
                );
            }
        }
    }

    #[test]
    fn errors_are_corrected_by_the_majority() {
        let g = genome(6_000, 3);
        let mut reads = tiling(&g, 6);
        // An error in read 2 at genome position 2,000, and an insertion in read 3 at 2,500.
        reads[2][1_200] = other(reads[2][1_200]);
        reads[3].insert(1_300, b'G');
        let set: Reads = reads.iter().map(Vec::as_slice).collect();
        let overlaps = find_overlaps(&set, &params()).expect("valid");
        let c = correct(&set, &overlaps, &CorrectionParams::default());
        assert_eq!(c.reads.get(ReadIndex(2)), &g[800..3_800]);
        assert_eq!(c.reads.get(ReadIndex(3)), &g[1_200..4_200]);
        assert_eq!(c.sites, 0);
        assert!(c.consistent.iter().all(|&x| x));
    }

    #[test]
    fn reads_of_another_copy_are_dropped_and_not_used() {
        // Two copies differing at one base, four reads of each, all over the same stretch.
        let a = genome(4_000, 5);
        let mut b = a.clone();
        b[2_000] = other(b[2_000]);
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for i in 0..4 {
            reads.push(a[i * 100..i * 100 + 3_000].to_vec());
        }
        for i in 0..4 {
            reads.push(b[i * 100 + 50..i * 100 + 3_050].to_vec());
        }
        let set: Reads = reads.iter().map(Vec::as_slice).collect();
        let overlaps = find_overlaps(&set, &params()).expect("valid");
        let c = correct(&set, &overlaps, &CorrectionParams::default());
        // Every read keeps its own copy's base.
        for (i, r) in reads.iter().enumerate() {
            assert_eq!(
                c.reads.get(ReadIndex(u32::try_from(i).expect("small"))),
                &r[..]
            );
        }
        // Overlaps within a copy agree, overlaps across copies do not.
        for (o, &ok) in overlaps.iter().zip(&c.consistent) {
            let same = (o.query.0 < 4) == (o.target.0 < 4);
            assert_eq!(ok, same, "{} {}", o.query.0, o.target.0);
        }
    }
}
