//! A vote of the reads at each base of the polished contigs, per strand.
//!
//! Partial-order consensus leaves some single-base errors where the reads mostly agree with
//! each other: where most reads have the right base but the rest disagree in several ways at
//! once (another base, a deletion), the graph's heaviest path can follow a minority. And at a
//! base the strain modifies (a methylated motif), the reads of one strand often call another
//! base, so a plain majority over both strands follows whichever strand has more reads. Here
//! the reads are mapped onto the polished contigs, and at each base each strand's reads are
//! counted by what they have there (each base, or none), and between bases by the single base
//! they insert. A base is replaced, removed or added when the reads of the two strands, each
//! strand's share averaged, favour that by a margin ([`MARGIN`]) over what the contig has. Bases
//! in runs of three are voted on only for another base (where the reads put an indel in a run is
//! arbitrary), and runs of four or more not at all: their lengths are called by
//! `homopolymer.rs`. Counts are 16-bit and saturate: a column under more than 65,535 reads of one
//! strand (only a very high-copy plasmid's) is voted on as if it had that many.

use crate::overlap::{EditOp, Overlap, ReadIndex, Reads};
use crate::seq::{Strand, complement};
use rayon::prelude::*;

use crate::motif::{Modified, Motifs};
use crate::polish::Draft;

/// A change to a draft's bases: `len` bases from `start` become `bases`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Change {
    pub(crate) start: u32,
    pub(crate) len: u32,
    pub(crate) bases: Vec<u8>,
}

/// A change is made when the reads favour it by this much (in each strand's share of its
/// reads, averaged over the strands) over what the contig has (the default
/// [`crate::BactParams::vote_margin`]). At a modified base the right base can have as little as
/// 40% of the reads, the rest split between another base and none: at 0.05 the votes fixed 6
/// of K. pneumoniae's 29 errors, at 0.15 4, and on four other genomes neither changed any.
pub(crate) const MARGIN: f64 = 0.05;

/// A strand's reads count only where at least this many cover the base.
const MIN_STRAND_DEPTH: u32 = 3;

/// What each strand's reads have at each base: A, C, G, T, nothing; and the single base they
/// insert after it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Column {
    /// `at[strand][kind]`, kind 0-3 for A, C, G, T and 4 for a deletion.
    pub(crate) at: [[u16; 5]; 2],
    /// `after[strand][base]`, base 0-3 for A, C, G, T.
    pub(crate) after: [[u16; 4]; 2],
}

/// Returns, per draft, the changes the reads vote for (see the module docs), from their
/// `mappings` onto the drafts ([`crate::overlap::map_reads`]); with [`Modified::Learn`], also
/// where one strand misreads a learned motif (see `motif.rs`), and how many motifs were
/// learned.
pub(crate) fn changes(
    drafts: &[Draft],
    reads: &Reads,
    mappings: &[Option<Overlap>],
    margin: f64,
    modified: Modified,
) -> (Vec<Vec<Change>>, usize) {
    let columns = columns(drafts, reads, mappings);
    let motifs = if modified == Modified::Learn {
        let seqs: Vec<&[u8]> = drafts.iter().map(|d| d.seq.as_slice()).collect();
        Motifs::learn(&seqs, &columns)
    } else {
        Motifs::default()
    };
    let changes = drafts
        .iter()
        .zip(&columns)
        .map(|(d, c)| vote(&d.seq, c, margin, &motifs))
        .collect();
    (changes, motifs.len())
}

/// Returns, per draft, what each strand's reads have at each of its bases, from their
/// `mappings` onto the drafts.
pub(crate) fn columns(
    drafts: &[Draft],
    reads: &Reads,
    mappings: &[Option<Overlap>],
) -> Vec<Vec<Column>> {
    // Each draft's mappings, in read order. Counts saturate, so they come out the same in any
    // order, and each stretch of a draft is counted on its own.
    let mut on: Vec<Vec<(ReadIndex, &Overlap)>> = vec![Vec::new(); drafts.len()];
    for (read, mapping) in reads.indices().zip(mappings) {
        if let Some(m) = mapping {
            on[m.target.index()].push((read, m));
        }
    }
    drafts
        .iter()
        .zip(&on)
        .map(|(d, on)| {
            let mut columns = vec![Column::default(); d.seq.len()];
            columns
                .par_chunks_mut(COLUMN_CHUNK)
                .enumerate()
                .for_each(|(k, chunk)| {
                    let from = k * COLUMN_CHUNK;
                    for &(read, m) in on {
                        count(chunk, from, reads.get(read), m);
                    }
                });
            columns
        })
        .collect()
}

/// Bases of a draft whose columns one task counts.
const COLUMN_CHUNK: usize = 1 << 16;

/// Adds what `read` (whose bases are `bases`) has at each base, by its mapping `m`, to
/// `columns`, the draft's columns from `from` on.
fn count(columns: &mut [Column], from: usize, bases: &[u8], m: &Overlap) {
    let to = from + columns.len();
    // An insertion before the mapping's first base counts at the base before it.
    if pos(m.target_end) <= from || pos(m.target_start) > to {
        return;
    }
    let strand = usize::from(m.strand.is_reverse());
    // The read's base at `q` as aligned.
    let at = |q: u32| match m.strand {
        Strand::Forward => bases[pos(q)],
        Strand::Reverse => complement(bases[bases.len() - 1 - pos(q)]),
    };
    let mut q = match m.strand {
        Strand::Forward => m.query_start,
        Strand::Reverse => {
            u32::try_from(bases.len()).expect("read lengths fit in u32") - m.query_end
        }
    };
    let mut t = m.target_start;
    for e in &m.edits {
        if pos(t) > to {
            break;
        }
        let n = e.len();
        // The run's columns inside this stretch.
        let inside = |t: u32| pos(t).max(from)..pos(t + n).min(to);
        match e.op() {
            EditOp::Match | EditOp::Mismatch => {
                for c in inside(t) {
                    if let Some(b) = code(at(q + u32::try_from(c - pos(t)).expect("within a run")))
                    {
                        bump(&mut columns[c - from].at[strand][b]);
                    }
                }
                t += n;
                q += n;
            }
            EditOp::Deletion => {
                for c in inside(t) {
                    bump(&mut columns[c - from].at[strand][4]);
                }
                t += n;
            }
            EditOp::Insertion => {
                if n == 1
                    && t > 0
                    && (from..to).contains(&(pos(t) - 1))
                    && let Some(b) = code(at(q))
                {
                    bump(&mut columns[pos(t) - 1 - from].after[strand][b]);
                }
                q += n;
            }
        }
    }
}

/// Returns the changes the columns of `seq` vote for, a strand's misreads of `motifs` first.
fn vote(seq: &[u8], columns: &[Column], margin: f64, motifs: &Motifs) -> Vec<Change> {
    let mut out = Vec::new();
    for p in 1..seq.len().saturating_sub(1) {
        let (before, here, next) = (seq[p - 1], seq[p], seq[p + 1]);
        let Some(current) = code(here) else { continue };
        let column = &columns[p];
        // A substitution is placed as well in a run as anywhere (an error that made the run,
        // AGCT read as AACT or AAGCT as AAACT, is one), so every base is voted on for another
        // base. A deletion only outside runs of three or more, which homopolymer calling
        // settles; in a run of two, a deletion split between the two places loses the vote.
        let in_triple = (p >= 2 && seq[p - 2] == here && before == here)
            || (p + 2 < seq.len() && seq[p + 2] == here && next == here)
            || (before == here && here == next);
        // Runs of homopolymer::MIN_RUN or more are homopolymer calling's alone.
        if run_len(seq, p) >= crate::homopolymer::MIN_RUN {
            continue;
        }
        if let Some(y) = motifs.corrected(seq, p, column) {
            out.push(Change {
                start: u32::try_from(p).expect("contigs are shorter than 4 Gb"),
                len: 1,
                bases: vec![BASES[y]],
            });
            continue;
        }
        let kinds = if in_triple { 4 } else { 5 };
        let share = |kind: usize| strand_share(column, |s| f64::from(column.at[s][kind]));
        if let (Some(have), Some((best, best_share))) = (
            share(current),
            (0..kinds)
                .filter_map(|k| share(k).map(|v| (k, v)))
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0))),
        ) && best != current
            && best_share >= have + margin
        {
            out.push(Change {
                start: u32::try_from(p).expect("contigs are shorter than 4 Gb"),
                len: 1,
                bases: if best == 4 {
                    Vec::new()
                } else {
                    vec![BASES[best]]
                },
            });
            continue;
        }
        // A base the reads insert after `p`, unlike both its neighbours (else it lengthens a
        // run, for homopolymer calling).
        for (b, &base) in BASES.iter().enumerate() {
            if base == here || base == next {
                continue;
            }
            let inserted = strand_share(column, |s| f64::from(column.after[s][b]));
            if inserted.is_some_and(|share| share >= 0.5 + margin / 2.0) {
                out.push(Change {
                    start: u32::try_from(p + 1).expect("contigs are shorter than 4 Gb"),
                    len: 0,
                    bases: vec![base],
                });
                break;
            }
        }
    }
    out
}

/// Returns the length of the run of one base that holds `seq[p]`.
fn run_len(seq: &[u8], p: usize) -> usize {
    let base = seq[p];
    let before = seq[..p].iter().rev().take_while(|&&b| b == base).count();
    let after = seq[p + 1..].iter().take_while(|&&b| b == base).count();
    before + 1 + after
}

/// Returns the share of each strand's reads at a column that `count` counts, averaged over
/// the strands with at least [`MIN_STRAND_DEPTH`] reads there; `None` if neither has.
fn strand_share(column: &Column, count: impl Fn(usize) -> f64) -> Option<f64> {
    let mut sum = 0.0;
    let mut strands = 0_u32;
    for (s, counts) in column.at.iter().enumerate() {
        let depth: u32 = counts.iter().map(|&n| u32::from(n)).sum();
        if depth >= MIN_STRAND_DEPTH {
            sum += count(s) / f64::from(depth);
            strands += 1;
        }
    }
    (strands > 0).then(|| sum / f64::from(strands))
}

/// Applies `changes` (per draft, sorted by start) to `drafts`; returns those applied, per draft.
/// A change that starts inside the one before is skipped.
pub(crate) fn apply(drafts: &mut [Draft], changes: &[Vec<Change>]) -> Vec<Vec<Applied>> {
    let mut applied = Vec::with_capacity(drafts.len());
    for (draft, changes) in drafts.iter_mut().zip(changes) {
        let mut made = Vec::new();
        if changes.is_empty() {
            applied.push(made);
            continue;
        }
        let mut seq = Vec::with_capacity(draft.seq.len() + changes.len());
        let mut at = 0_usize;
        for c in changes {
            let start = pos(c.start);
            if start < at {
                continue;
            }
            seq.extend_from_slice(&draft.seq[at..start]);
            seq.extend_from_slice(&c.bases);
            at = start + pos(c.len);
            made.push(Applied {
                start: c.start,
                len: c.len,
                new_len: u32::try_from(c.bases.len()).expect("a change is shorter than 4 Gb"),
            });
        }
        seq.extend_from_slice(&draft.seq[at..]);
        draft.seq = seq;
        applied.push(made);
    }
    applied
}

/// A change [`apply`] made: `len` bases from `start` became `new_len` bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Applied {
    pub(crate) start: u32,
    pub(crate) len: u32,
    pub(crate) new_len: u32,
}

/// Where the places of a draft (its bases, and its end) are after the changes [`apply`] made to
/// it: bases inserted at or before a place move it on, and a place inside a changed stretch
/// keeps its offset into the new bases, up to their end.
pub(crate) struct Moves<'a> {
    /// The changes, in order (as [`apply`] returns them: none overlaps another, and at one
    /// start an insertion comes before a change of bases).
    applied: &'a [Applied],
    /// Per change, how many bases the changes before it added (less those they took out); then
    /// all of them.
    shift: Vec<i64>,
}

impl<'a> Moves<'a> {
    pub(crate) fn new(applied: &'a [Applied]) -> Self {
        let mut shift = Vec::with_capacity(applied.len() + 1);
        let mut sum = 0_i64;
        shift.push(sum);
        for a in applied {
            sum += i64::from(a.new_len) - i64::from(a.len);
            shift.push(sum);
        }
        Self { applied, shift }
    }

    /// Where place `at` (a base, or the end) is now.
    pub(crate) fn at(&self, at: u32) -> u32 {
        let place = |x: i64| u32::try_from(x).expect("a place on a draft is under 4 Gb");
        // The changes before `at`: those starting before it, and insertions just before it.
        let before = self
            .applied
            .partition_point(|a| a.start < at || (a.start == at && a.len == 0));
        if let Some(a) = before.checked_sub(1).map(|i| self.applied[i])
            && at < a.start + a.len
        {
            let into = i64::from((at - a.start).min(a.new_len));
            return place(i64::from(a.start) + self.shift[before - 1] + into);
        }
        place(i64::from(at) + self.shift[before])
    }
}

pub(crate) const BASES: [u8; 4] = *b"ACGT";

pub(crate) fn code(base: u8) -> Option<usize> {
    BASES.iter().position(|&b| b == base)
}

fn bump(n: &mut u16) {
    *n = n.saturating_add(1);
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::{Edit, map_reads};
    use crate::seq::reverse_complement;

    use crate::params::{BactParams, ReadKind};

    /// [`columns`] one mapping at a time over each whole draft.
    fn columns_in_one_pass(
        drafts: &[Draft],
        reads: &Reads,
        mappings: &[Option<Overlap>],
    ) -> Vec<Vec<Column>> {
        let mut columns: Vec<Vec<Column>> = drafts
            .iter()
            .map(|d| vec![Column::default(); d.seq.len()])
            .collect();
        for (read, mapping) in reads.indices().zip(mappings) {
            let Some(m) = mapping else { continue };
            let strand = usize::from(m.strand.is_reverse());
            let (oriented, mut q) = match m.strand {
                Strand::Forward => (reads.get(read).to_vec(), m.query_start),
                Strand::Reverse => (
                    reverse_complement(reads.get(read)),
                    reads.read_len(read) - m.query_end,
                ),
            };
            let column = &mut columns[m.target.index()];
            let mut t = m.target_start;
            for e in &m.edits {
                let n = e.len();
                match e.op() {
                    EditOp::Match | EditOp::Mismatch => {
                        for k in 0..n {
                            if let Some(b) = code(oriented[pos(q + k)]) {
                                bump(&mut column[pos(t + k)].at[strand][b]);
                            }
                        }
                        t += n;
                        q += n;
                    }
                    EditOp::Deletion => {
                        for k in 0..n {
                            bump(&mut column[pos(t + k)].at[strand][4]);
                        }
                        t += n;
                    }
                    EditOp::Insertion => {
                        if n == 1
                            && t > 0
                            && let Some(b) = code(oriented[pos(q)])
                        {
                            bump(&mut column[pos(t - 1)].after[strand][b]);
                        }
                        q += n;
                    }
                }
            }
        }
        columns
    }

    #[test]
    fn columns_counted_by_stretches_are_those_counted_in_one_pass() {
        let mut state = 17_u64;
        let mut next = move |n: u64| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % n
        };
        let base = |r: u64| b"ACGT"[usize::try_from(r).expect("below 4")];
        // Two drafts, one over three stretches long; noisy reads of both strands across the
        // stretches' ends.
        let drafts: Vec<Draft> = [150_000, 9_000]
            .iter()
            .map(|&len| Draft::new((0..len).map(|_| base(next(4))).collect(), false))
            .collect();
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for (d, starts) in [(0, 0..60), (1, 0..4)] {
            for i in starts {
                let seq = &drafts[d].seq;
                let start = (i * 2_400 + 300).min(seq.len() - 5_000);
                let mut read = Vec::new();
                for &b in &seq[start..start + 5_000] {
                    match next(40) {
                        0 => read.push(base(next(4))),
                        1 => {}
                        2 => read.extend_from_slice(&[b, base(next(4))]),
                        _ => read.push(b),
                    }
                }
                if i % 7 == 3 {
                    // Bases of nowhere at both ends: a mapping clipped at both.
                    let junk = |next: &mut dyn FnMut(u64) -> u64| -> Vec<u8> {
                        (0..400).map(|_| base(next(4))).collect()
                    };
                    read = [junk(&mut next), read, junk(&mut next)].concat();
                }
                if i % 2 == 1 {
                    read = reverse_complement(&read);
                }
                reads.push(read);
            }
        }
        let reads: Reads = reads.iter().map(Vec::as_slice).collect();
        let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
        let mut mappings = map_reads(&targets, &reads, &BactParams::new(ReadKind::Ont).map)
            .expect("valid settings");
        assert!(mappings.iter().filter(|m| m.is_some()).count() > 50);
        // And one mapping that starts and ends with a one-base insertion.
        let clipped = mappings
            .iter_mut()
            .zip(reads.indices())
            .find_map(|(m, r)| {
                m.as_mut()
                    .filter(|m| m.query_start > 0 && m.query_end < reads.read_len(r))
            })
            .expect("a mapping clipped at both ends");
        let mut edits = vec![Edit::new(EditOp::Insertion, 1)];
        edits.extend_from_slice(&clipped.edits);
        edits.push(Edit::new(EditOp::Insertion, 1));
        clipped.edits = edits;
        clipped.query_start -= 1;
        clipped.query_end += 1;
        assert_eq!(
            columns(&drafts, &reads, &mappings),
            columns_in_one_pass(&drafts, &reads, &mappings)
        );
    }

    fn column(fwd: [u16; 5], rev: [u16; 5]) -> Column {
        Column {
            at: [fwd, rev],
            after: [[0; 4]; 2],
        }
    }

    #[test]
    fn a_base_most_reads_of_both_strands_disagree_with_is_replaced() {
        let seq = b"TAGCT";
        let mut columns = vec![column([0; 5], [0; 5]); seq.len()];
        // Base 2 is G; forward reads say G 21, A 2; reverse reads say G 6, A 18.
        columns[2] = column([0, 0, 21, 0, 0], [2, 0, 6, 0, 0]);
        assert!(
            vote(seq, &columns, MARGIN, &Motifs::default()).is_empty(),
            "G is the strands' choice"
        );
        // Contig has A there instead: the strands' average favours G.
        let seq = b"TAACT";
        columns[2] = column([2, 0, 21, 0, 0], [18, 0, 6, 0, 0]);
        // Base 2 in a run of two As is voted on too.
        assert_eq!(
            vote(seq, &columns, MARGIN, &Motifs::default()),
            [Change {
                start: 2,
                len: 1,
                bases: b"G".to_vec()
            }]
        );
        // In a run of three, too: a substitution's place is not in doubt.
        let seq = b"TAAAC";
        assert_eq!(
            vote(seq, &columns, MARGIN, &Motifs::default()),
            [Change {
                start: 2,
                len: 1,
                bases: b"G".to_vec()
            }]
        );
        // But not a deletion.
        let mut deleted = columns.clone();
        deleted[2] = column([2, 0, 0, 0, 21], [6, 0, 0, 0, 18]);
        assert!(vote(seq, &deleted, MARGIN, &Motifs::default()).is_empty());
        let seq = b"TCACT";
        assert_eq!(
            vote(seq, &columns, MARGIN, &Motifs::default()),
            [Change {
                start: 2,
                len: 1,
                bases: b"G".to_vec()
            }]
        );
    }

    #[test]
    fn deletions_and_insertions_follow_the_reads() {
        let seq = b"ACGTA";
        let mut columns = vec![column([10, 10, 10, 10, 0], [10, 10, 10, 10, 0]); seq.len()];
        columns[2] = column([1, 0, 2, 0, 9], [0, 0, 3, 0, 9]);
        columns[3].after = [[0, 8, 0, 0], [0, 7, 0, 0]];
        columns[3].at = [[0, 0, 0, 10, 0], [0, 0, 0, 10, 0]];
        let c = vote(seq, &columns, MARGIN, &Motifs::default());
        assert_eq!(
            c,
            [
                Change {
                    start: 2,
                    len: 1,
                    bases: Vec::new()
                },
                Change {
                    start: 4,
                    len: 0,
                    bases: b"C".to_vec()
                },
            ]
        );
        let mut drafts = vec![Draft::new(seq.to_vec(), false)];
        let applied = apply(&mut drafts, &[c]);
        assert_eq!(applied[0].len(), 2);
        assert_eq!(drafts[0].seq, b"ACTCA");
        // ACGTA became ACTCA: G deleted, C inserted before the last A.
        let moves = Moves::new(&applied[0]);
        let now: Vec<u32> = (0..=5).map(|at| moves.at(at)).collect();
        assert_eq!(now, [0, 1, 2, 2, 4, 5]);
    }

    /// Applies `changes` (start, bases taken out, bases put in) to `seq` and checks that every
    /// base no change touched is found where [`Moves`] says it went.
    fn moves_follow(seq: &[u8], changes: &[(u32, u32, &[u8])]) -> (Vec<u8>, Vec<u32>) {
        let changes: Vec<Change> = changes
            .iter()
            .map(|&(start, len, bases)| Change {
                start,
                len,
                bases: bases.to_vec(),
            })
            .collect();
        let mut drafts = vec![Draft::new(seq.to_vec(), false)];
        let applied = apply(&mut drafts, &[changes]);
        let moves = Moves::new(&applied[0]);
        let now: Vec<u32> = (0..=u32::try_from(seq.len()).expect("short"))
            .map(|at| moves.at(at))
            .collect();
        let touched = |at: u32| {
            applied[0]
                .iter()
                .any(|a| a.start <= at && at < a.start + a.len)
        };
        for (at, &to) in now.iter().enumerate().take(seq.len()) {
            let at32 = u32::try_from(at).expect("short");
            if !touched(at32) {
                assert_eq!(
                    drafts[0].seq[to as usize], seq[at],
                    "base {at} moved to {to}"
                );
            }
        }
        assert_eq!(
            *now.last().expect("the end"),
            u32::try_from(drafts[0].seq.len()).expect("short"),
            "the end moves to the new end"
        );
        assert!(
            now.windows(2).all(|w| w[0] <= w[1]),
            "places keep their order: {now:?}"
        );
        (drafts[0].seq.clone(), now)
    }

    #[test]
    fn places_move_with_every_kind_of_change() {
        let seq = b"ACGTACGTAC";
        // At the first base, a base made three, three made one, an insertion at the end.
        let (new, now) = moves_follow(
            seq,
            &[(0, 1, b"TT"), (3, 1, b"GGG"), (5, 3, b"A"), (10, 0, b"CC")],
        );
        assert_eq!(new, b"TTCGGGGAAACCC");
        assert_eq!(now, [0, 2, 3, 4, 7, 8, 9, 9, 9, 10, 13]);
        // A place inside a stretch made shorter stops at its new end; one inside a stretch
        // made longer keeps its offset.
        assert_eq!(now[6], 9);
    }

    #[test]
    fn an_insertion_before_a_change_at_one_place_moves_it_past_the_insertion() {
        let (new, now) = moves_follow(b"AAAA", &[(2, 0, b"G"), (2, 1, b"T")]);
        assert_eq!(new, b"AAGTA");
        assert_eq!(now, [0, 1, 3, 4, 5]);
    }

    #[test]
    fn no_changes_move_nothing() {
        let (_, now) = moves_follow(b"ACGT", &[]);
        assert_eq!(now, [0, 1, 2, 3, 4]);
    }
}
