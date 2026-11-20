//! Base-level alignment of two similar sequences in O(ND) time (Landau and Vishkin 1989).
//!
//! For each edit count `d` and diagonal `k` (query position minus target position), the aligner
//! keeps the furthest target position reachable with `d` edits, then slides along matching
//! bases. Unit costs: a mismatch, an insertion and a deletion each cost one edit. HiFi reads
//! differ at about one base in a thousand, so `d` stays small and the work with it. The furthest
//! positions of every `d` are kept for the traceback, so memory grows as `d` squared; alignments
//! past [`DiffAligner::max_edits`] are given up, and so, with a [`Stall`] set, are alignments
//! whose front has stopped moving (one sequence has run into bases the other lacks). `N`
//! matches nothing, itself included.

/// One kind of alignment column, in target order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EditOp {
    /// A target base and a query base that are equal.
    Match,
    /// A target base and a query base that differ.
    Mismatch,
    /// A query base the target lacks.
    Insertion,
    /// A target base the query lacks.
    Deletion,
}

/// A run of columns of one [`EditOp`], in four bytes (an overlap holds tens of runs, and an
/// assembly millions of overlaps).
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Edit(u32);

impl Edit {
    /// The most columns one run holds; a longer stretch of one op takes several runs.
    pub const MAX_LEN: u32 = u32::MAX >> 2;

    /// Returns a run of `len` columns of `op`.
    ///
    /// # Panics
    ///
    /// Panics if `len` is zero or more than [`Edit::MAX_LEN`].
    #[must_use]
    pub fn new(op: EditOp, len: u32) -> Self {
        assert!(
            (1..=Self::MAX_LEN).contains(&len),
            "a run of {len} columns is out of range"
        );
        let code = match op {
            EditOp::Match => 0,
            EditOp::Mismatch => 1,
            EditOp::Insertion => 2,
            EditOp::Deletion => 3,
        };
        Self(len << 2 | code)
    }

    /// Returns the kind of the columns.
    #[must_use]
    pub const fn op(self) -> EditOp {
        match self.0 & 3 {
            0 => EditOp::Match,
            1 => EditOp::Mismatch,
            2 => EditOp::Insertion,
            _ => EditOp::Deletion,
        }
    }

    /// Returns the number of columns, at least one.
    #[must_use]
    #[expect(clippy::len_without_is_empty, reason = "a run is never empty")]
    pub const fn len(self) -> u32 {
        self.0 >> 2
    }
}

impl std::fmt::Debug for Edit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Edit")
            .field("op", &self.op())
            .field("len", &self.len())
            .finish()
    }
}

/// Where an extension alignment stopped: the bases of each sequence it covers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub struct Extent {
    /// Target bases covered, from the start.
    pub target: usize,
    /// Query bases covered, from the start.
    pub query: usize,
}

/// When an alignment is given up before [`DiffAligner::max_edits`]: its last `edits` edits moved
/// its front (the furthest target plus query position reached) fewer than `advance` bases.
///
/// Two unrelated sequences move the front a few bases per edit, two copies of one sequence
/// hundreds; an indel of `n` bases takes `n` edits that hardly move it, so `edits` must be more
/// than the longest indel the alignments should cross.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Stall {
    /// The edits looked back over.
    pub edits: usize,
    /// The bases, target and query together, the front must move over them.
    pub advance: usize,
}

/// Reusable buffers for [`DiffAligner::global`] and [`DiffAligner::extend`].
#[derive(Clone, Debug)]
pub struct DiffAligner {
    max_edits: usize,
    stall: Option<Stall>,
    max_lag: Option<usize>,
    /// Furthest target positions: the diagonals `bands[d]` (first and last) of edit count `d`
    /// from `offsets[d]`, each edit count's row with at least [`PAD`] unreached diagonals
    /// either side.
    furthest: Vec<i64>,
    offsets: Vec<usize>,
    bands: Vec<(i64, i64)>,
    /// The front after each edit count: the most target plus query bases any diagonal covers.
    fronts: Vec<i64>,
    /// The traceback's runs, last first.
    runs: Vec<Edit>,
}

/// How an alignment ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    /// Both sequences to their ends.
    Global,
    /// Either sequence to its end, from the start of both.
    Extend,
}

impl DiffAligner {
    /// Returns an aligner that gives up past `max_edits` edits.
    #[must_use]
    pub fn new(max_edits: usize) -> Self {
        Self {
            max_edits,
            stall: None,
            max_lag: None,
            furthest: Vec::new(),
            offsets: Vec::new(),
            bands: Vec::new(),
            fronts: Vec::new(),
            runs: Vec::new(),
        }
    }

    /// Sets the edit count past which alignments are given up.
    pub fn set_max_edits(&mut self, max_edits: usize) {
        self.max_edits = max_edits;
    }

    /// Returns the edit count past which alignments are given up.
    #[must_use]
    pub const fn max_edits(&self) -> usize {
        self.max_edits
    }

    /// Sets when alignments are given up before [`DiffAligner::max_edits`] (`None`: never).
    pub fn set_stall(&mut self, stall: Option<Stall>) {
        self.stall = stall;
    }

    /// Returns when alignments are given up before [`DiffAligner::max_edits`].
    #[must_use]
    pub const fn stall(&self) -> Option<Stall> {
        self.stall
    }

    /// Sets how far (target plus query bases) a diagonal may lag the furthest one before it is
    /// dropped (`None`: none is, and the alignment found has the fewest edits).
    pub fn set_max_lag(&mut self, max_lag: Option<usize>) {
        self.max_lag = max_lag;
    }

    /// Aligns all of `query` to all of `target` and appends the alignment to `edits`.
    ///
    /// Returns the edit distance, or `None` (leaving `edits` as it was) if it is more than
    /// [`DiffAligner::max_edits`] or the alignment stalls.
    pub fn global(&mut self, target: &[u8], query: &[u8], edits: &mut Vec<Edit>) -> Option<usize> {
        let (d, end) = self.forward(target, query, Mode::Global)?;
        self.traceback(target, query, d, end, edits);
        Some(d)
    }

    /// Aligns `query` to `target` from the start of both until either ends, and appends the
    /// alignment to `edits`.
    ///
    /// Returns the edit count and the bases of each sequence covered, or `None` (leaving `edits`
    /// as it was) if neither end is reached within [`DiffAligner::max_edits`] edits or the
    /// alignment stalls first. Of the diagonals that reach an end with the fewest edits, the one
    /// nearest the main diagonal wins.
    pub fn extend(
        &mut self,
        target: &[u8],
        query: &[u8],
        edits: &mut Vec<Edit>,
    ) -> Option<(usize, Extent)> {
        let (d, end) = self.forward(target, query, Mode::Extend)?;
        self.traceback(target, query, d, end, edits);
        Some((d, end))
    }

    /// Fills the furthest positions up to the first edit count that ends, and returns it with
    /// the end.
    fn forward(&mut self, target: &[u8], query: &[u8], mode: Mode) -> Option<(usize, Extent)> {
        // With at most `max_edits` edits an extension covers at most that many target bases
        // more than query bases: past one more, which it never reaches, a long target (the rest
        // of a contig) is never read, and is left out of the scan for N below.
        let target = match mode {
            Mode::Global => target,
            Mode::Extend => &target[..target.len().min(query.len() + self.max_edits + 1)],
        };
        // Only a target's N stops a slide that matching bases would not (a query's N differs
        // from the target's base anyway); most sequences have none.
        let has_n = target.contains(&b'N');
        let (n, m) = (signed(target.len()), signed(query.len()));
        self.furthest.clear();
        self.offsets.clear();
        self.bands.clear();
        self.fronts.clear();
        for d in 0..=self.max_edits {
            // Each edit count's row is written in place, after the last, between two unreached
            // diagonals either side ([`PAD`]), so that the next reads each diagonal's three
            // neighbours without a range check and without a copy.
            let base = self.furthest.len();
            let mut front = 0;
            let lo;
            if d == 0 {
                let i = slide(target, query, 0, 0, has_n);
                self.furthest.extend_from_slice(&[UNREACHED; PAD]);
                self.furthest.push(i);
                self.furthest.extend_from_slice(&[UNREACHED; PAD]);
                front = 2 * i;
                lo = 0;
            } else {
                // The diagonals of edit count d - 1, one more on each side. The best of the
                // steps that stay inside both sequences, as `step` takes them; an unreached
                // diagonal offers a negative position.
                let (before_lo, before_hi) = self.bands[d - 1];
                let start = self.offsets[d - 1];
                let width = index(before_hi - before_lo + 1);
                lo = before_lo - 1;
                self.furthest
                    .resize(base + PAD + width + 2 + PAD, UNREACHED);
                let (done, new) = self.furthest.split_at_mut(base);
                let before = &done[start - PAD..start + width + PAD];
                let row = &mut new[PAD..PAD + width + 2];
                for ((k, around), cell) in (lo..).zip(before.windows(3)).zip(row) {
                    let (below, same, above) = (around[0], around[1], around[2]);
                    let mismatch = if same < n && same + k < m {
                        same + 1
                    } else {
                        same
                    };
                    let deletion = if above >= 0 && above < n {
                        above + 1
                    } else {
                        UNREACHED
                    };
                    let insertion = if below + k <= m { below } else { UNREACHED };
                    let i = mismatch.max(deletion).max(insertion);
                    if i >= 0 {
                        let i = slide(target, query, i, k, has_n);
                        front = front.max(2 * i + k);
                        *cell = i;
                    }
                }
            }
            let row_start = base + PAD;
            let row_len = self.furthest.len() - PAD - row_start;
            if let Some(end) = finished(
                &self.furthest[row_start..row_start + row_len],
                lo,
                n,
                m,
                mode,
            ) {
                self.offsets.push(row_start);
                self.bands.push((lo, lo + signed(row_len) - 1));
                return Some((d, end));
            }
            // Diagonals that lag the front by more than `max_lag` are dropped from the band's
            // edges: the alignment hardly ever comes back through them.
            let max_lag = self.max_lag;
            let lagging = |i: i64, k: i64| {
                i < 0 || max_lag.is_some_and(|lag| 2 * i + k + signed(lag) < front)
            };
            let row = &mut self.furthest[row_start..row_start + row_len];
            let mut first = 0;
            let mut last = row.len();
            while last - first > 1 && lagging(row[first], lo + signed(first)) {
                first += 1;
            }
            while last - first > 1 && lagging(row[last - 1], lo + signed(last) - 1) {
                last -= 1;
            }
            // The dropped diagonals are unreached to the next edit count.
            row[..first].fill(UNREACHED);
            row[last..].fill(UNREACHED);
            let none_reached = row[first..last].iter().all(|&i| i < 0);
            self.offsets.push(row_start + first);
            self.bands.push((lo + signed(first), lo + signed(last) - 1));
            self.fronts.push(front);
            if let Some(stall) = self.stall
                && d >= stall.edits
                && front - self.fronts[d - stall.edits] < signed(stall.advance)
            {
                return None;
            }
            if none_reached {
                return None;
            }
        }
        None
    }

    /// Returns the furthest target position on diagonal `k` with edit count `d` before sliding
    /// along matches, and the step from edit count `d - 1` that reaches it, if any does. Only
    /// steps that stay inside both sequences count; ties go to the first of no edit, a mismatch,
    /// a deletion and an insertion.
    fn step(&self, d: usize, k: i64, n: i64, m: i64) -> Option<(i64, Step)> {
        let mut best: Option<(i64, Step)> = None;
        let mut offer = |i: i64, step: Step| {
            if best.is_none_or(|(b, _)| i > b) {
                best = Some((i, step));
            }
        };
        if let Some(f) = self.at(d - 1, k) {
            offer(f, Step::Same);
            if f < n && f + k < m {
                offer(f + 1, Step::Mismatch);
            }
        }
        if let Some(f) = self.at(d - 1, k + 1)
            && f < n
        {
            offer(f + 1, Step::Deletion);
        }
        if let Some(f) = self.at(d - 1, k - 1)
            && f + k <= m
        {
            offer(f, Step::Insertion);
        }
        best
    }

    /// Returns the furthest target position on diagonal `k` at edit count `d`, if it was
    /// reached.
    fn at(&self, d: usize, k: i64) -> Option<i64> {
        let (lo, hi) = self.bands[d];
        if k < lo || k > hi {
            return None;
        }
        let i = self.furthest[self.offsets[d] + index(k - lo)];
        (i >= 0).then_some(i)
    }

    /// Appends the columns of the alignment ending at `end` with `d` edits to `edits`.
    fn traceback(
        &mut self,
        target: &[u8],
        query: &[u8],
        d: usize,
        end: Extent,
        edits: &mut Vec<Edit>,
    ) {
        let (n, m) = (signed(target.len()), signed(query.len()));
        let mut runs = std::mem::take(&mut self.runs);
        runs.clear();
        let mut k = signed(end.query) - signed(end.target);
        let mut i = signed(end.target);
        for e in (1..=d).rev() {
            let (start, step) = self
                .step(e, k, n, m)
                .expect("a reached diagonal has a step");
            push(&mut runs, EditOp::Match, columns(i - start));
            i = start;
            match step {
                Step::Same => {}
                Step::Mismatch => {
                    push(&mut runs, EditOp::Mismatch, 1);
                    i -= 1;
                }
                Step::Deletion => {
                    push(&mut runs, EditOp::Deletion, 1);
                    i -= 1;
                    k += 1;
                }
                Step::Insertion => {
                    push(&mut runs, EditOp::Insertion, 1);
                    k -= 1;
                }
            }
        }
        debug_assert_eq!(k, 0, "the traceback ends off the main diagonal");
        push(&mut runs, EditOp::Match, columns(i));
        for e in runs.iter().rev() {
            push(edits, e.op(), e.len());
        }
        self.runs = runs;
    }
}

/// A diagonal no step reached: far enough below zero that a step from it stays negative.
const UNREACHED: i64 = i64::MIN / 4;

/// Unreached diagonals kept either side of each edit count's row in [`DiffAligner`]'s furthest
/// positions: the next row reads one beyond each of its own ends.
const PAD: usize = 2;

/// Returns where the alignment ends if the edit count just filled into `row`, whose first
/// diagonal is `lo`, ends it (`n` and `m` the target's and the query's lengths).
fn finished(row: &[i64], lo: i64, n: i64, m: i64, mode: Mode) -> Option<Extent> {
    match mode {
        Mode::Global => {
            let k = m - n;
            (k >= lo && k < lo + signed(row.len()) && row[index(k - lo)] == n).then(|| Extent {
                target: index(n),
                query: index(m),
            })
        }
        Mode::Extend => row
            .iter()
            .zip(lo..)
            .filter(|&(&i, k)| i >= 0 && (i == n || i + k == m))
            // Nearest the main diagonal first, the negative diagonal before the positive.
            .min_by_key(|&(_, k)| (k.abs(), k > 0))
            .map(|(&i, k)| Extent {
                target: index(i),
                query: index(i + k),
            }),
    }
}

/// The step that reaches a diagonal's furthest position from the edit count before.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Step {
    /// The same diagonal, with no further edit.
    Same,
    /// The same diagonal, one mismatch on.
    Mismatch,
    /// The diagonal above, one target base on.
    Deletion,
    /// The diagonal below, one query base on.
    Insertion,
}

/// Appends `len` columns of `op` to `edits`, extending the last run when it has the same op (up
/// to [`Edit::MAX_LEN`] columns a run).
pub fn push(edits: &mut Vec<Edit>, op: EditOp, mut len: u32) {
    if let Some(last) = edits.last_mut()
        && last.op() == op
    {
        let more = len.min(Edit::MAX_LEN - last.len());
        *last = Edit::new(op, last.len() + more);
        len -= more;
    }
    while len > 0 {
        let run = len.min(Edit::MAX_LEN);
        edits.push(Edit::new(op, run));
        len -= run;
    }
}

/// Slides from target position `i` along diagonal `k` while the bases match, eight at a time;
/// `has_n` says whether `target` holds an N (which matches nothing, itself included).
fn slide(target: &[u8], query: &[u8], i: i64, k: i64, has_n: bool) -> i64 {
    /// The top bit of each byte.
    const HIGH: u64 = 0x8080_8080_8080_8080;
    /// The other bits of each byte.
    const LOW: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    const N: u64 = u64::from_le_bytes([b'N'; 8]);
    /// The top bit of each byte of `v` that is not zero.
    fn nonzero(v: u64) -> u64 {
        (((v & LOW) + LOW) | v) & HIGH
    }
    let word =
        |s: &[u8], at: usize| u64::from_le_bytes(s[at..at + 8].try_into().expect("eight bytes"));
    let (mut t, mut q) = (index(i), index(i + k));
    while t + 8 <= target.len() && q + 8 <= query.len() {
        let a = word(target, t);
        // A base that differs, or an N, stops the slide.
        let mut stop = nonzero(a ^ word(query, q));
        if has_n {
            stop |= !nonzero(a ^ N) & HIGH;
        }
        if stop != 0 {
            return signed(t + usize::try_from(stop.trailing_zeros() / 8).expect("below 8"));
        }
        t += 8;
        q += 8;
    }
    while t < target.len()
        && q < query.len()
        && target[t] == query[q]
        && (!has_n || target[t] != b'N')
    {
        t += 1;
        q += 1;
    }
    signed(t)
}

/// Returns a run length of `x` columns.
fn columns(x: i64) -> u32 {
    u32::try_from(x).expect("a run within one sequence")
}

fn signed(x: usize) -> i64 {
    i64::try_from(x).expect("sequence lengths fit in i64")
}

fn index(x: i64) -> usize {
    usize::try_from(x).expect("positions are non-negative")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Returns the target and query bases the alignment covers, and its edit count.
    fn apply(target: &[u8], query: &[u8], edits: &[Edit]) -> (usize, usize, usize) {
        let (mut t, mut q, mut d) = (0, 0, 0);
        for e in edits {
            let len = e.len() as usize;
            match e.op() {
                EditOp::Match => {
                    assert_eq!(target[t..t + len], query[q..q + len], "match at {t}/{q}");
                    t += len;
                    q += len;
                }
                EditOp::Mismatch => {
                    for x in 0..len {
                        assert_ne!(target[t + x], query[q + x], "mismatch at {t}/{q}");
                    }
                    t += len;
                    q += len;
                    d += len;
                }
                EditOp::Insertion => {
                    q += len;
                    d += len;
                }
                EditOp::Deletion => {
                    t += len;
                    d += len;
                }
            }
        }
        (t, q, d)
    }

    /// Edit distance by the quadratic dynamic programme, for checking.
    fn distance(a: &[u8], b: &[u8]) -> usize {
        let mut row: Vec<usize> = (0..=b.len()).collect();
        for (i, &x) in a.iter().enumerate() {
            let mut prev = row[0];
            row[0] = i + 1;
            for (j, &y) in b.iter().enumerate() {
                let cur = row[j + 1];
                row[j + 1] = (prev + usize::from(x != y)).min(row[j] + 1).min(cur + 1);
                prev = cur;
            }
        }
        row[b.len()]
    }

    #[test]
    fn global_finds_the_edit_distance() {
        let cases: [(&[u8], &[u8]); 7] = [
            (b"ACGTACGT", b"ACGTACGT"),
            (b"ACGTACGT", b"ACGAACGT"),
            (b"ACGTACGT", b"ACGTTACGT"),
            (b"ACGTACGT", b"ACGACGT"),
            (b"", b"ACG"),
            (b"ACG", b""),
            (b"GATTACAGATTACA", b"GATACAGGATTTACA"),
        ];
        let mut aligner = DiffAligner::new(10);
        for (target, query) in cases {
            let mut edits = Vec::new();
            let d = aligner
                .global(target, query, &mut edits)
                .expect("within the cap");
            assert_eq!(d, distance(target, query), "{target:?} {query:?}");
            assert_eq!(
                apply(target, query, &edits),
                (target.len(), query.len(), d),
                "{target:?} {query:?} {edits:?}"
            );
        }
    }

    #[test]
    fn global_matches_the_quadratic_distance_on_random_pairs() {
        // A fixed linear congruential sequence: no randomness crate needed.
        let mut state = 12345_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            state >> 33
        };
        let mut aligner = DiffAligner::new(64);
        for _ in 0..200 {
            let len = 20 + usize::try_from(next() % 60).expect("small");
            let target: Vec<u8> = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let mut query = target.clone();
            for _ in 0..next() % 6 {
                let at = usize::try_from(next()).expect("small") % (query.len() + 1);
                match next() % 3 {
                    0 if at < query.len() => query[at] = b"ACGT"[(next() % 4) as usize],
                    1 => query.insert(at, b"ACGT"[(next() % 4) as usize]),
                    _ if at < query.len() => {
                        query.remove(at);
                    }
                    _ => {}
                }
            }
            let mut edits = Vec::new();
            let d = aligner
                .global(&target, &query, &mut edits)
                .expect("within the cap");
            assert_eq!(d, distance(&target, &query));
            assert_eq!(
                apply(&target, &query, &edits),
                (target.len(), query.len(), d)
            );
        }
    }

    #[test]
    fn n_matches_nothing() {
        let mut aligner = DiffAligner::new(10);
        let mut edits = Vec::new();
        assert_eq!(aligner.global(b"ANA", b"ANA", &mut edits), Some(1));
        assert_eq!(
            edits,
            [
                Edit::new(EditOp::Match, 1),
                Edit::new(EditOp::Mismatch, 1),
                Edit::new(EditOp::Match, 1),
            ]
        );
    }

    #[test]
    fn global_gives_up_past_the_cap() {
        let mut aligner = DiffAligner::new(2);
        let mut edits = Vec::new();
        assert_eq!(aligner.global(b"AAAAAA", b"CCCCCC", &mut edits), None);
        assert!(edits.is_empty());
    }

    #[test]
    fn extend_stops_at_the_shorter_end() {
        let mut aligner = DiffAligner::new(10);
        let mut edits = Vec::new();
        // The query ends first: everything after it in the target is left out.
        let (d, end) = aligner
            .extend(b"ACGTTGCAGGT", b"ACGATGC", &mut edits)
            .expect("within the cap");
        assert_eq!(
            (d, end),
            (
                1,
                Extent {
                    target: 7,
                    query: 7
                }
            )
        );
        assert_eq!(apply(b"ACGTTGCAGGT", b"ACGATGC", &edits), (7, 7, 1));
    }

    #[test]
    fn extend_takes_an_indel_near_the_end() {
        let mut aligner = DiffAligner::new(10);
        let mut edits = Vec::new();
        let (d, end) = aligner
            .extend(b"ACGTACGTAC", b"ACGTTACGTACGG", &mut edits)
            .expect("within the cap");
        assert_eq!(d, 1);
        assert_eq!(
            end,
            Extent {
                target: 10,
                query: 11
            }
        );
    }

    /// `len` random bases from a fixed linear congruential sequence.
    /// `seq` with an error at about one base in `every` (a substitution, an insertion or a
    /// deletion in turn).
    fn noisy(seq: &[u8], every: u64, seed: u64) -> Vec<u8> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            state >> 33
        };
        let mut out = Vec::with_capacity(seq.len() + seq.len() / 4);
        for &b in seq {
            if next() % every != 0 {
                out.push(b);
                continue;
            }
            let other = b"ACGT"[usize::try_from(next() % 4).expect("small")];
            match next() % 3 {
                0 => out.push(if other == b {
                    b"CGTA"[usize::try_from(next() % 4).expect("small")]
                } else {
                    other
                }),
                1 => out.extend_from_slice(&[b, other]),
                _ => {}
            }
        }
        out
    }

    /// The edit counts of `query` aligned globally to `target` with no band and with `lag`.
    fn banded_and_not(target: &[u8], query: &[u8], lag: usize) -> (usize, usize) {
        let mut aligner = DiffAligner::new(2_000);
        let mut edits = Vec::new();
        let full = aligner
            .global(target, query, &mut edits)
            .expect("within the cap");
        aligner.set_max_lag(Some(lag));
        edits.clear();
        let banded = aligner
            .global(target, query, &mut edits)
            .expect("within the cap");
        assert_eq!(
            apply(target, query, &edits),
            (target.len(), query.len(), banded)
        );
        (full, banded)
    }

    #[test]
    fn a_band_keeps_noisy_alignments_near_the_fewest_edits() {
        // Two reads of one stretch, each with an error at one base in six.
        let genome = random(3_000, 71);
        let (a, b) = (noisy(&genome, 6, 72), noisy(&genome, 6, 73));
        let (full, banded) = banded_and_not(&a, &b, 150);
        assert!(full > 600, "{full} edits");
        assert!(
            banded >= full && banded <= full + full / 50,
            "{banded} against {full}"
        );
    }

    #[test]
    fn a_band_crosses_an_indel_of_up_to_half_its_width() {
        // 60 bases one read lacks crossed as the fewest edits do, inside a band of 150; 200
        // bases are not (653 edits against 518), but the reads are still aligned.
        let genome = random(4_000, 81);
        let a = noisy(&genome, 20, 82);
        for (gap, exact) in [(60, true), (200, false)] {
            let lacking: Vec<u8> = genome[..2_000]
                .iter()
                .chain(&genome[2_000 + gap..])
                .copied()
                .collect();
            let b = noisy(&lacking, 20, 83);
            let (full, banded) = banded_and_not(&a, &b, 150);
            assert_eq!(banded == full, exact, "gap {gap}: {banded} against {full}");
        }
    }

    fn random(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
            })
            .collect()
    }

    const STALL: Stall = Stall {
        edits: 64,
        advance: 512,
    };

    #[test]
    fn extend_never_reaches_past_the_query_by_more_than_its_cap() {
        // A query of a stretch of a much longer target, deletions and noise and all: the
        // target cut `max_edits + 1` bases past the query's length aligns the same (mapping
        // relies on it to copy only that much of a contig).
        let target = random(20_000, 91);
        for (seed, max_edits, band) in [(92, 400, None), (93, 150, Some(64)), (94, 40, None)] {
            let query = noisy(&target[..1_500], 12, seed);
            let mut aligner = DiffAligner::new(max_edits);
            aligner.set_max_lag(band);
            aligner.set_stall(Some(STALL));
            let (mut whole, mut cut) = (Vec::new(), Vec::new());
            let on_whole = aligner.extend(&target, &query, &mut whole);
            let reach = query.len() + max_edits + 1;
            let on_cut = aligner.extend(&target[..reach], &query, &mut cut);
            assert_eq!(on_whole, on_cut, "seed {seed}");
            assert_eq!(whole, cut, "seed {seed}");
        }
    }

    #[test]
    fn a_stall_gives_up_soon_after_the_sequences_part() {
        // 1 kb in common, then unrelated: the front stops moving at 2 kb (target plus query).
        let shared = random(1_000, 1);
        let mut target = shared.clone();
        target.extend(random(10_000, 2));
        let mut query = shared;
        query.extend(random(10_000, 3));
        let mut aligner = DiffAligner::new(1_000);
        aligner.set_stall(Some(STALL));
        let mut edits = Vec::new();
        assert_eq!(aligner.extend(&target, &query, &mut edits), None);
        assert!(edits.is_empty());
        assert_eq!(aligner.fronts.len(), STALL.edits + 1);
    }

    #[test]
    fn a_stall_still_crosses_an_indel_shorter_than_its_edits() {
        let (left, right) = (random(2_000, 4), random(2_000, 5));
        let mut target = left.clone();
        target.extend_from_slice(&right);
        let mut query = left;
        query.extend(random(50, 6));
        query.extend_from_slice(&right);
        let mut aligner = DiffAligner::new(1_000);
        aligner.set_stall(Some(STALL));
        let mut edits = Vec::new();
        let (d, end) = aligner
            .extend(&target, &query, &mut edits)
            .expect("crosses the insertion");
        assert_eq!(d, 50);
        assert_eq!(
            end,
            Extent {
                target: 4_000,
                query: 4_050
            }
        );
    }

    #[test]
    fn slide_stops_where_base_by_base_would() {
        let slow = |target: &[u8], query: &[u8], i: usize, k: usize| {
            let (mut t, mut q) = (i, i + k);
            while t < target.len() && q < query.len() && target[t] == query[q] && target[t] != b'N'
            {
                t += 1;
                q += 1;
            }
            t
        };
        let base = random(40, 7);
        for at in 0..base.len() {
            for change in *b"NAC" {
                let mut target = base.clone();
                let mut query = base.clone();
                // A mismatch, or an N in either sequence, at every offset.
                if change == b'N' {
                    target[at] = b'N';
                } else if query[at] == change {
                    query[at] = b'G';
                } else {
                    query[at] = change;
                }
                for from in [0, 3, 9] {
                    for cut in [0, 5] {
                        let q = &query[..query.len() - cut];
                        assert_eq!(
                            index(slide(&target, q, signed(from), 0, target.contains(&b'N'))),
                            slow(&target, q, from, 0),
                            "at {at}, {change}, from {from}, cut {cut}"
                        );
                    }
                }
            }
        }
        let query: Vec<u8> = [&b"TT"[..], &base].concat();
        assert_eq!(
            index(slide(&base, &query, 0, 2, false)),
            slow(&base, &query, 0, 2)
        );
    }

    #[test]
    fn a_run_keeps_its_op_and_length_in_four_bytes() {
        assert_eq!(size_of::<Edit>(), 4);
        for op in [
            EditOp::Match,
            EditOp::Mismatch,
            EditOp::Insertion,
            EditOp::Deletion,
        ] {
            for len in [1, 2, 1_000, Edit::MAX_LEN] {
                let e = Edit::new(op, len);
                assert_eq!((e.op(), e.len()), (op, len));
            }
        }
    }

    #[test]
    fn push_starts_another_run_past_the_longest() {
        let mut edits = Vec::new();
        push(&mut edits, EditOp::Match, Edit::MAX_LEN - 1);
        push(&mut edits, EditOp::Match, 3);
        assert_eq!(
            edits,
            [
                Edit::new(EditOp::Match, Edit::MAX_LEN),
                Edit::new(EditOp::Match, 2)
            ]
        );
    }

    #[test]
    fn push_merges_runs() {
        let mut edits = Vec::new();
        push(&mut edits, EditOp::Match, 3);
        push(&mut edits, EditOp::Match, 2);
        push(&mut edits, EditOp::Insertion, 0);
        push(&mut edits, EditOp::Deletion, 1);
        assert_eq!(
            edits,
            [Edit::new(EditOp::Match, 5), Edit::new(EditOp::Deletion, 1),]
        );
    }
}
