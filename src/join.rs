//! Joining linear drafts whose ends overlap, and closing those whose end runs into their start.
//!
//! The string graph stops a contig wherever its last read has no overlap on: an overlap missed (a
//! read's stretch of errors), or cut as transitive or short beside a wrong one. Where that happens
//! the reads on both sides still overlap, so the two contigs (or one contig's two ends) share tens
//! of kilobases. Here the drafts' ends are overlapped as if they were reads: each linear draft
//! gives its first and its last [`END_LEN`] bases (or halves), and a dovetail between the outer
//! ends of two of these pieces joins the two drafts' ends. A join is made when its two ends have no
//! other, or when its overlap is [`JOIN_DOMINANCE`] times as long as each other join of either end
//! (two ends that stop in copies of one repeat overlap along the repeat only; the true overlap runs
//! on into unique sequence): the drafts are merged along the overlap (each draft's bases up to
//! where the next one's start lies), and a chain of drafts that comes back to its first is a
//! circle. The few bases the overlap's coordinates are off by are mended by polishing. First, a
//! short linear draft that holds a small replicon two or more times over is cut to one circle of it
//! ([`collapse_tandem`]). Drafts far shallower than the longest ([`MIN_JOIN_DEPTH_SHARE`]) join
//! no other (they may still close on themselves): a few reads of another sample, or a chimera's,
//! overlapping a chromosome's ends would make its joins ambiguous.

use std::collections::{BTreeMap, BTreeSet};

use crate::consensus::PoaBackend;
use crate::overlap::{
    EditOp, InvalidOverlapParamsError, Overlap, OverlapParams, ReadIndex, Reads, find_overlaps,
    map_reads,
};
use crate::seq::{Strand, reverse_complement};

use crate::hang;
use crate::polish::Draft;

/// A linear draft whose mean depth is below this share of the longest draft's is not joined
/// (it is left out as shallow later; see [`crate::MIN_LINEAR_DEPTH_SHARE`]).
const MIN_JOIN_DEPTH_SHARE: f64 = 0.25;

/// A join is made despite rival joins at its ends when its overlap is at least this many times
/// as long as each rival's (as a ratio, numerator over denominator).
const JOIN_DOMINANCE: (u32, u32) = (5, 4);

/// Each linear draft's ends are overlapped over this many bases.
pub(crate) const END_LEN: usize = 150_000;

/// One end of a draft: its index, and whether it is the draft's end (its last base) rather
/// than its start.
type End = (usize, Side);

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Side {
    Start,
    End,
}

/// A join of two draft ends: reading on out of `from`, the sequence runs through `between`
/// into `to`, the first `to_overlap` bases of which it already holds.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Join {
    from: End,
    to: End,
    /// The bases of `to`'s draft, counted from its joined end, that the overlap covers.
    to_overlap: u32,
    /// Bases between the two ends (a gap reads crossed), read out of `from`.
    between: Vec<u8>,
}

/// Joins `drafts` along overlaps of their ends (see the module docs), one round; returns the
/// drafts after joining and how many joins were made. The bridging comes next, then
/// [`join_more`].
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
pub(crate) fn join_ends(
    drafts: Vec<Draft>,
    params: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    join_ends_once(drafts, params, Ends::Raw)
}

/// Joins the drafts' ends again, round after round until no join is left: a join's rival may
/// be gone once another join took its end (two plasmids sharing a long repeat each overlap the
/// other's ends along it, and only the one whose own overlap is clearly longer closes at
/// first; pieces of a chromosome wait on each other likewise). Run after the bridging, which
/// first closes a draft through another that holds a variant of the stretch between its ends
/// and leaves that one out: joined to one of the ends, it took the chromosome's circle
/// (*E. hormaechei* 4). A round that finds no join with the drafts' raw ends is tried again
/// without them ([`Draft::shallow_ends`]): a raw end runs on past the overlap with a read's
/// unpolished bases, over a kilobase of them too many to snap (a *Burkholderia* chromosome left
/// linear across a 13 kb overlap; a *N. meningitidis* chromosome left in three pieces whose
/// ends each overlapped one other by 5-8 kb). Returns the drafts after joining and how many
/// joins were made.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
pub(crate) fn join_more(
    drafts: Vec<Draft>,
    params: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let mut drafts = drafts;
    let mut made = 0;
    loop {
        let mut joined;
        (drafts, joined) = join_ends_once(drafts, params, Ends::Raw)?;
        if joined == 0 {
            (drafts, joined) = join_ends_once(drafts, params, Ends::Polished)?;
        }
        if joined == 0 {
            return Ok((drafts, made));
        }
        made += joined;
    }
}

/// Unique bases either side of a repeat that joins are weighed over: the joined sequence a read
/// is mapped onto holds this many before the overlap and as many after it.
const SPAN_FLANK: usize = 2_000;

/// A read spans a join's overlap if it aligns over at least this many of the unique bases on
/// each side of it.
pub(crate) const MIN_SPAN_FLANK: usize = 500;

/// Two ends that overlap only along a repeat (each overlapping others as long) are joined if at
/// least this many reads span the overlap and the unique bases either side, and...
/// One is enough where no read spans a rival: HiFi reads of 5.5 kb on an *X. oryzae* chromosome
/// left in three pieces ending in 5 kb repeats, one read spanned one of the four rival joins
/// and none the others; with it the chromosome closed as the reference has it.
pub(crate) const MIN_SPANNING_READS: usize = 1;

/// ...this many times as many as span any rival join's.
pub(crate) const SPANNING_DOMINANCE: usize = 2;

/// Joins the linear `drafts`' ends that overlap only along copies of one repeat, by the reads
/// that span it; returns the drafts after joining and how many joins were made.
///
/// A chromosome whose reads are shorter than a few of its repeats (HiFi reads of 5-8 kb; a
/// 10 kb repeat in three copies) lays out as pieces that each end in a copy: every end overlaps
/// each of the others along the whole repeat, no join dominates ([`JOIN_DOMINANCE`]), and the
/// reads mapped onto the ends' last bases ([`join_by_reads`]) fall on any copy. Here each such
/// join is spelled out with [`SPAN_FLANK`] unique bases either side of its overlap and the
/// reads are mapped onto all of them: a read that runs from one side's unique bases across the
/// repeat into the other's ([`MIN_SPAN_FLANK`]) aligns that far only on the join it came from.
/// A join with [`MIN_SPANNING_READS`] such reads and [`SPANNING_DOMINANCE`] times as many as
/// each rival is made; the ends it leaves may then have no rival ([`join_more`]).
///
/// # Errors
///
/// Returns an error if `params` or `map` cannot be used.
pub(crate) fn join_across_repeats(
    drafts: Vec<Draft>,
    reads: &Reads,
    params: &OverlapParams,
    map: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let drafts: Vec<Draft> = drafts.into_iter().map(collapse_tandem).collect();
    let joins = end_joins(&drafts, params, Ends::Raw)?;
    let settled = dominant(&joins);
    let shares_end = |a: &Join, b: &Join| [b.from, b.to].iter().any(|&e| e == a.from || e == a.to);
    // The joins with a rival, other than those a dominant join leaves no room for.
    let tied: Vec<&Join> = joins
        .iter()
        .filter(|j| !settled.contains(j) && !settled.iter().any(|s| shares_end(j, s)))
        .collect();
    if tied.is_empty() {
        return Ok((drafts, 0));
    }
    // Each join spelled out: its flank before the overlap, the overlap, and its flank after.
    let spelled: Vec<(Vec<u8>, usize)> = tied.iter().map(|j| spell(&drafts, j)).collect();
    let windows: Reads = spelled.iter().map(|(seq, _)| seq.as_slice()).collect();
    let mut spanning = vec![0_usize; tied.len()];
    for m in map_reads(&windows, reads, map)?.into_iter().flatten() {
        let k = m.target.index();
        let (seq, before) = (&spelled[k].0, spelled[k].1);
        let after = before + pos(tied[k].to_overlap);
        if pos(m.target_start) + MIN_SPAN_FLANK <= before
            && pos(m.target_end) >= after + MIN_SPAN_FLANK
            && after + MIN_SPAN_FLANK <= seq.len()
        {
            spanning[k] += 1;
        }
    }
    let chosen: Vec<Join> = tied
        .iter()
        .enumerate()
        .filter(|&(k, j)| {
            log::debug!(
                "bacteria: ends {:?} and {:?} overlap by {} along a repeat; {} reads span it",
                j.from,
                j.to,
                j.to_overlap,
                spanning[k]
            );
            spanning[k] >= MIN_SPANNING_READS
                && tied.iter().enumerate().all(|(r, rival)| {
                    r == k
                        || !shares_end(j, rival)
                        || spanning[r] * SPANNING_DOMINANCE <= spanning[k]
                })
        })
        .map(|(_, &j)| j.clone())
        .collect();
    let made = chosen.len();
    Ok((merge(drafts, &chosen), made))
}

/// Returns `join` spelled out (see [`join_across_repeats`]): up to [`SPAN_FLANK`] bases of its
/// first draft before the overlap, the overlap, and up to as many of its second draft's after;
/// and how many bases lie before the overlap.
fn spell(drafts: &[Draft], join: &Join) -> (Vec<u8>, usize) {
    debug_assert!(join.between.is_empty(), "the ends of an end join overlap");
    // Out of `from` (its joined end last), then into `to` (its joined end first); the overlap is
    // the last bases of the one and the first of the other.
    let out = oriented(&drafts[join.from.0].seq, join.from.1 == Side::Start);
    let into = oriented(&drafts[join.to.0].seq, join.to.1 == Side::End);
    let overlap = (join.to_overlap as usize).min(into.len()).min(out.len());
    let before = SPAN_FLANK.min(out.len() - overlap);
    let after = SPAN_FLANK.min(into.len() - overlap);
    let mut seq = out[out.len() - overlap - before..].to_vec();
    seq.extend_from_slice(&into[overlap..overlap + after]);
    (seq, before)
}

/// Closes the linear drafts whose ends overlap with a stretch in one of them only (see
/// [`closed_by_anchors`]); returns the drafts and how many were closed. Run after the other
/// joins and the bridging: a draft that holds a variant of another's stretch can have ends
/// that meet too, and is the one the other is closed through (*E. hormaechei* 4).
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
pub(crate) fn close_by_anchors(
    drafts: Vec<Draft>,
    params: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let mut drafts = drafts;
    let mut closed = 0;
    for i in 0..drafts.len() {
        if drafts[i].circular {
            continue;
        }
        let Some(seq) = closed_by_anchors(drafts[i].deep_part(), params)? else {
            continue;
        };
        // Ends in copies of a repeat close the draft through whichever copy the anchor maps to.
        if !(crate::bridge::end_unique(&drafts, i, false, params)?
            && crate::bridge::end_unique(&drafts, i, true, params)?)
        {
            log::debug!(
                "bacteria: a draft of {} bases not closed through its end anchors: an end lies \
                 in a repeat",
                drafts[i].seq.len()
            );
            continue;
        }
        log::debug!(
            "bacteria: a draft of {} bases closed through its end anchors into {}",
            drafts[i].seq.len(),
            seq.len()
        );
        let depth = drafts[i].depth;
        drafts[i] = Draft::new(seq, true);
        drafts[i].depth = depth;
        closed += 1;
    }
    Ok((drafts, closed))
}

/// Two drafts' ends overlap only if their alignment has substitutions at no more than this share of the
/// columns where both have polished bases (before its overhangs are carried): polished drafts
/// of one replicon agree to well under a percent, while copies of a family of plasmids
/// (*Borrelia*'s cp32s, 4-15% apart) overlap along kilobases too, longer than a plasmid's own
/// closing overlap, and took its join. The raw reads past a dead end are not held to it: one
/// with a poor stretch (often why the layout stopped there) differed at 4-9% from the end it
/// overlapped, over 4-26 kb, in a *B. arboris* chromosome whose joins were all right.
///
/// Only substitutions are counted: the cp32s differ by them (2.3-12%), while drafts polished
/// from short reads (a *V. parahaemolyticus* with reads of 1.2 kb) differ by 1-2% in indels
/// on top of about 1% in substitutions, and their joins were lost.
const MAX_JOIN_DIVERGENCE: f64 = 0.02;

/// The cap is put to the polished columns alone where there are at least this many; to all
/// the columns otherwise (a 1.9 kb overlap almost all of raw read closed a *Borrelia* cp32 with
/// 5 kb of another member in it).
const MIN_POLISHED_COLUMNS: u64 = 1_000;

/// How a join's overlap is judged against [`MAX_JOIN_DIVERGENCE`].
#[derive(Clone, Copy, Debug, PartialEq)]
enum Divergence {
    /// Its share of substitutions where both ends have polished bases.
    Polished(f64),
    /// Its share of substitutions in all its columns: too few are polished to tell (raw reads
    /// past dead ends overlap).
    Raw(f64),
}

impl Divergence {
    fn share(self) -> f64 {
        match self {
            Self::Polished(x) | Self::Raw(x) => x,
        }
    }
}

/// The substitutions of `o`, an overlap between two of `pieces`, in its columns where both
/// pieces' bases lie in their `polished` ranges (in all its columns, if there are too few such
/// to tell).
fn divergence(o: &Overlap, pieces: &Reads, polished: &[(usize, usize)]) -> Divergence {
    let (q_from, q_to) = polished[o.query.index()];
    let (t_from, t_to) = polished[o.target.index()];
    let q_len = pos(pieces.read_len(o.query));
    // The query is walked along the target's strand; its bases are counted back to its own.
    let q_own = |x: usize| match o.strand {
        Strand::Forward => x,
        Strand::Reverse => q_len - 1 - x,
    };
    let mut q = match o.strand {
        Strand::Forward => pos(o.query_start),
        Strand::Reverse => q_len - pos(o.query_end),
    };
    let mut t = pos(o.target_start);
    let (mut columns, mut substitutions) = (0_u64, 0_u64);
    let mut all_substitutions = 0_u64;
    for e in &o.edits {
        if e.op() == EditOp::Mismatch {
            all_substitutions += u64::from(e.len());
        }
        for _ in 0..e.len() {
            let qq = q_own(q.min(q_len - 1));
            if (q_from..q_to).contains(&qq) && (t_from..t_to).contains(&t) {
                columns += 1;
                substitutions += u64::from(e.op() == EditOp::Mismatch);
            }
            match e.op() {
                EditOp::Match | EditOp::Mismatch => (q, t) = (q + 1, t + 1),
                EditOp::Insertion => q += 1,
                EditOp::Deletion => t += 1,
            }
        }
    }
    #[expect(clippy::cast_precision_loss, reason = "a share, not an exact count")]
    if columns < MIN_POLISHED_COLUMNS {
        Divergence::Raw(all_substitutions as f64 / o.columns().max(1) as f64)
    } else {
        Divergence::Polished(substitutions as f64 / columns as f64)
    }
}

/// Which of `overlaps` (between `pieces`, with their `polished` ranges) may be joins: those
/// within [`MAX_JOIN_DIVERGENCE`], and those over it that are mostly raw reads and the least
/// divergent overlap at both their pieces. A raw read's substitutions alone can pass the cap
/// (reads of a *V. parahaemolyticus* differ at 3-7% over 1.6-2 kb overlaps of three joins it
/// needs, and no other overlap competes); where a relative competes (*Borrelia* cp32s, at
/// 2-10% over as short overlaps), the plasmid's own overlap is the closer.
fn within_cap(overlaps: &[Overlap], pieces: &Reads, polished: &[(usize, usize)]) -> Vec<bool> {
    let shares: Vec<Divergence> = overlaps
        .iter()
        .map(|o| divergence(o, pieces, polished))
        .collect();
    // The least divergence of an overlap at each piece.
    let mut least = vec![f64::INFINITY; pieces.len()];
    for (o, d) in overlaps.iter().zip(&shares) {
        for p in [o.query.index(), o.target.index()] {
            least[p] = least[p].min(d.share());
        }
    }
    overlaps
        .iter()
        .zip(&shares)
        .map(|(o, &d)| match d {
            Divergence::Polished(x) => x <= MAX_JOIN_DIVERGENCE,
            Divergence::Raw(x) => {
                x <= MAX_JOIN_DIVERGENCE
                    || [o.query.index(), o.target.index()]
                        .iter()
                        .all(|&p| x <= least[p])
            }
        })
        .collect()
}

/// A draft's end anchors are this long: longer than an rRNA operon or an IS element, so that
/// an anchor maps whole only where the draft's other end runs on.
const ANCHOR: usize = 10_000;

/// A draft's end anchor maps whole if at most this many of its bases are left out.
const ANCHOR_CLIP: u32 = 500;

/// End anchors are taken this far in from the ends, the first that maps whole: bases of
/// another replicon that chimeric reads put on chromosomes' ends kept anchors there from
/// mapping (1.9 kb on a *Burkholderia* sp. chromosome's start; 5 kb and 21 kb on a
/// *B. arboris* chromosome's two ends, 7.9 kb on another's, whose ends overlap by 13 kb).
const ANCHOR_SKIPS: [usize; 11] = [
    0, 2_500, 5_000, 7_500, 10_000, 12_500, 15_000, 17_500, 20_000, 22_500, 25_000,
];

/// Returns `deep` (a linear draft's polished bases) closed into a circle, if its ends overlap
/// but differ inside the overlap, as when one holds an element that the other lacks (part of
/// the cells have it): see [`anchored_cut`], with the draft as both sides. The circle starts
/// past where the end's last base lies in the start, so it keeps the end's version of the
/// overlap; the majority's is set when loci are repaired.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
fn closed_by_anchors(
    deep: &[u8],
    params: &OverlapParams,
) -> Result<Option<Vec<u8>>, InvalidOverlapParamsError> {
    Ok(anchored_cut(deep, deep, params)?
        .filter(|&(cut, trim)| cut < deep.len() - trim)
        .map(|(cut, trim)| deep[cut..deep.len() - trim].to_vec()))
}

/// Where in `right` the last base of `left` lies, and how many of `left`'s last bases are
/// left out, if `left`'s end overlaps `right`'s start through end anchors (both polished
/// bases, at least four anchors long): `left`'s bases up to that end then `right`'s past that
/// place are the two joined. An [`ANCHOR`] of `left`'s last bases is mapped onto the first
/// [`END_LEN`] of `right` (or half), and one of `right`'s first onto the last of `left`, each
/// at the end or, past bases a chimeric read put there, [`ANCHOR_SKIPS`] in (the first that
/// maps whole): both must map whole, forward, and the overlaps they imply (`right`'s bases up
/// to where `left`'s anchored end lies, and `left`'s from where `right`'s first one lies) must
/// each be at least an anchor long and within half of each other.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
fn anchored_cut(
    left: &[u8],
    right: &[u8],
    params: &OverlapParams,
) -> Result<Option<(usize, usize)>, InvalidOverlapParamsError> {
    let (left_len, right_len) = (left.len(), right.len());
    if left_len.min(right_len) < 4 * ANCHOR {
        return Ok(None);
    }
    let (left_n, right_n) = ((left_len / 2).min(END_LEN), (right_len / 2).min(END_LEN));
    // Each anchor onto the other side's end only: its own end holds it too.
    let onto =
        |region: &[u8], anchor: &[u8]| -> Result<Option<Overlap>, InvalidOverlapParamsError> {
            let regions: Reads = [region].into_iter().collect();
            let anchors: Reads = [anchor].into_iter().collect();
            Ok(map_reads(&regions, &anchors, params)?.pop().flatten())
        };
    let anchor = u32::try_from(ANCHOR).expect("anchors are short");
    let whole =
        |m: &Overlap| !m.strand.is_reverse() && m.query_end - m.query_start + ANCHOR_CLIP >= anchor;
    let mut end = None;
    for skip in ANCHOR_SKIPS
        .into_iter()
        .filter(|&skip| left_len >= 4 * ANCHOR + skip)
    {
        let from = left_len - ANCHOR - skip;
        let found = onto(&right[..right_n], &left[from..from + ANCHOR])?;
        if let Some(m) = found.filter(|m| whole(m)) {
            end = Some((m, skip));
            break;
        }
    }
    let Some((end, end_skip)) = end else {
        return Ok(None);
    };
    let mut start = None;
    for skip in ANCHOR_SKIPS
        .into_iter()
        .filter(|&skip| right_len >= 4 * ANCHOR + skip)
    {
        let found = onto(&left[left_len - left_n..], &right[skip..skip + ANCHOR])?;
        if let Some(m) = found.filter(|m| whole(m)) {
            start = Some((m, skip));
            break;
        }
    }
    let Some((start, start_skip)) = start else {
        return Ok(None);
    };
    // Where `left`'s anchored end lies in `right`, and `right`'s first base (before its
    // anchor) in `left`'s anchored end.
    let cut = pos(end.target_end + (anchor - end.query_end));
    let in_right = cut;
    let in_left = (left_n - end_skip).saturating_sub(
        pos(start.target_start.saturating_sub(start.query_start)).saturating_sub(start_skip),
    );
    let shorter = in_right.min(in_left);
    let ok = shorter >= ANCHOR && in_right.abs_diff(in_left) * 2 <= in_right.max(in_left);
    Ok(ok.then_some((cut, end_skip)))
}

/// Joins linear drafts whose ends overlap through end anchors ([`anchored_cut`]), where no
/// overlap between the ends was found: kilobases of another replicon that a chimeric read put
/// on an end (3.3 kb on one of a *B. arboris* chromosome's two pieces, 1.4 kb on the other)
/// leave the overlap with an overhang at the wrong side. An end joins only if it has one such
/// partner. Returns the drafts and how many joins were made.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
pub(crate) fn join_by_anchors(
    drafts: Vec<Draft>,
    params: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let mut drafts = drafts;
    let mut joined = 0;
    loop {
        let longest = drafts.iter().map(|d| d.seq.len()).max().unwrap_or(0);
        let least = drafts
            .iter()
            .find(|d| d.seq.len() == longest)
            .map_or(0.0, |d| MIN_JOIN_DEPTH_SHARE * d.depth);
        let lines: Vec<usize> = (0..drafts.len())
            .filter(|&i| {
                let d = &drafts[i];
                !d.circular
                    && d.hairpins == [false; 2]
                    && d.depth >= least
                    && d.deep_part().len() >= 4 * ANCHOR
            })
            .collect();
        // Each line forward or turned round: (draft, turned) of the left side's end and of the
        // right side's start, and where the left's anchored end lies in the right and how many
        // bases past it the left loses.
        type Found = ((usize, bool), (usize, bool), (usize, usize));
        let mut found: Vec<Found> = Vec::new();
        for (k, &i) in lines.iter().enumerate() {
            for &j in &lines[k + 1..] {
                for (ti, tj) in [(false, false), (false, true), (true, false), (true, true)] {
                    let left = oriented(drafts[i].deep_part(), ti);
                    let right = oriented(drafts[j].deep_part(), tj);
                    if let Some(cut) = anchored_cut(&left, &right, params)? {
                        found.push(((i, ti), (j, tj), cut));
                    }
                }
            }
        }
        // The ends each join takes: a line's end is its start if it is turned round.
        let ends = |&((i, ti), (j, tj), _): &Found| [(i, !ti), (j, tj)];
        // Ends in copies of a repeat join whichever draft holds the copy the anchor maps to.
        let mut unique: BTreeMap<(usize, bool), bool> = BTreeMap::new();
        for f in &found {
            for (d, at_end) in ends(f) {
                if let std::collections::btree_map::Entry::Vacant(v) = unique.entry((d, at_end)) {
                    v.insert(crate::bridge::end_unique(&drafts, d, at_end, params)?);
                }
            }
        }
        found.retain(|f| ends(f).iter().all(|e| unique[e]));
        let mut uses: BTreeMap<(usize, bool), usize> = BTreeMap::new();
        for f in &found {
            for e in ends(f) {
                *uses.entry(e).or_default() += 1;
            }
        }
        let Some(&((i, ti), (j, tj), (cut, trim))) =
            found.iter().find(|f| ends(f).iter().all(|e| uses[e] == 1))
        else {
            return Ok((drafts, joined));
        };
        let left = oriented(drafts[i].deep_part(), ti);
        let right = oriented(drafts[j].deep_part(), tj);
        log::debug!(
            "bacteria: drafts of {} and {} bases joined through end anchors",
            left.len(),
            right.len()
        );
        let mut seq = left;
        seq.truncate(seq.len() - trim);
        seq.extend_from_slice(&right[cut..]);
        #[expect(clippy::cast_precision_loss, reason = "a mean depth")]
        let depth = (drafts[i].depth * drafts[i].seq.len() as f64
            + drafts[j].depth * drafts[j].seq.len() as f64)
            / (drafts[i].seq.len() + drafts[j].seq.len()) as f64;
        let mut merged = Draft::new(seq, false);
        merged.depth = depth;
        drafts = drafts
            .into_iter()
            .enumerate()
            .filter(|&(k, _)| k != i && k != j)
            .map(|(_, d)| d)
            .collect();
        drafts.push(merged);
        joined += 1;
    }
}

/// `seq`, turned round if `turned`.
fn oriented(seq: &[u8], turned: bool) -> Vec<u8> {
    if turned {
        reverse_complement(seq)
    } else {
        seq.to_vec()
    }
}

/// Which of a draft's ends [`join_ends_once`] overlaps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Ends {
    /// Its bases to the last, with the raw reads past a dead end.
    Raw,
    /// Its polished bases only ([`Draft::deep_part`]).
    Polished,
}

/// One round of [`join_ends`].
fn join_ends_once(
    drafts: Vec<Draft>,
    params: &OverlapParams,
    ends: Ends,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let drafts: Vec<Draft> = drafts.into_iter().map(collapse_tandem).collect();
    let joins = end_joins(&drafts, params, ends)?;
    let mut unique = dominant(&joins);
    unique.sort_unstable_by_key(|j| (j.from, j.to));
    for j in &unique {
        log::debug!(
            "bacteria: ends {:?} and {:?} of drafts of {} and {} bases overlap by {}",
            j.from,
            j.to,
            drafts[j.from.0].seq.len(),
            drafts[j.to.0].seq.len(),
            j.to_overlap
        );
    }
    let mut drafts = drafts;
    if ends == Ends::Polished {
        drop_shallow_ends(&mut drafts, &unique);
    }
    let made = unique.len();
    Ok((merge(drafts, &unique), made))
}

/// Every join the overlaps of the linear `drafts`' `ends` make, rivals and all (see the module
/// docs).
fn end_joins(
    drafts: &[Draft],
    params: &OverlapParams,
    ends: Ends,
) -> Result<Vec<Join>, InvalidOverlapParamsError> {
    let longest = (0..drafts.len()).max_by_key(|&i| (drafts[i].seq.len(), std::cmp::Reverse(i)));
    let least = longest.map_or(0.0, |i| MIN_JOIN_DEPTH_SHARE * drafts[i].depth);
    let deep = |i: usize| drafts[i].depth >= least;
    let linear: Vec<usize> = (0..drafts.len()).filter(|&i| !drafts[i].circular).collect();
    if linear.is_empty() {
        return Ok(Vec::new());
    }
    // Pieces 2k and 2k + 1 are the start and the end of linear draft k.
    let mut pieces = Reads::new();
    for &i in &linear {
        let seq = match ends {
            Ends::Raw => drafts[i].seq.as_slice(),
            Ends::Polished => drafts[i].deep_part(),
        };
        // A short draft's two halves, so that its end and its start are two pieces.
        let n = (seq.len() / 2).min(END_LEN);
        pieces.push(&seq[..n]);
        pieces.push(&seq[seq.len() - n..]);
    }
    let overlaps = find_overlaps(&pieces, params)?;
    // Where each piece's polished bases lie in it: the raw reads past a dead end may differ
    // from the other end at several percent (a read with a poor stretch, which is why the
    // layout stopped there).
    let polished: Vec<(usize, usize)> = linear
        .iter()
        .flat_map(|&i| {
            let d = &drafts[i];
            let len = d.seq.len();
            let n = (len / 2).min(END_LEN);
            let [s0, s1] = match ends {
                Ends::Raw if d.shallow_ends[0] + d.shallow_ends[1] < len => d.shallow_ends,
                _ => [0; 2],
            };
            [(s0.min(n), n), (0, n.saturating_sub(s1))]
        })
        .collect();
    let capped = within_cap(&overlaps, &pieces, &polished);
    let joins: Vec<Join> = overlaps
        .iter()
        .zip(capped)
        .filter(|&(o, ok)| {
            log::debug!(
                "bacteria: pieces {} {}..{} and {} {}..{} ({:?}) overlap: {} columns, {} \
                 differences{}",
                o.query.index(),
                o.query_start,
                o.query_end,
                o.target.index(),
                o.target_start,
                o.target_end,
                o.strand,
                o.columns(),
                o.differences(),
                if ok { "" } else { ", too divergent" }
            );
            ok
        })
        .map(|(o, _)| o)
        .filter_map(|o| hang::snap(&pieces, o))
        .filter_map(|o| join_of(&o, &pieces, &linear))
        // A shallow draft may close on itself (a plasmid of few copies), not join another.
        .filter(|j| j.from.0 == j.to.0 || (deep(j.from.0) && deep(j.to.0)))
        .collect();
    Ok(joins)
}

/// Takes the shallow bases off the ends `joins` join, which were measured from the polished
/// ones.
fn drop_shallow_ends(drafts: &mut [Draft], joins: &[Join]) {
    for j in joins {
        for (d, side) in [j.from, j.to] {
            let draft = &mut drafts[d];
            let n = draft.shallow_ends[usize::from(side == Side::End)];
            match side {
                Side::Start => drop(draft.seq.drain(..n.min(draft.seq.len()))),
                Side::End => draft.seq.truncate(draft.seq.len().saturating_sub(n)),
            }
            draft.shallow_ends[usize::from(side == Side::End)] = 0;
        }
    }
}

/// Reads that bridge two draft ends are found in the last this many polished bases of each
/// (or half of it).
const BRIDGE_PIECE: usize = 5_000;

/// A read runs on past a draft's end if at least this many of its bases lie beyond it.
const MIN_OUTWARD: u32 = 300;

/// A read's alignment reaches a draft's end if it stops at most this many bases short of it.
const END_SLACK: u32 = 100;

/// Two ends are joined through the reads that cross the gap between them if at least this
/// many do...
const MIN_BRIDGE_READS: usize = 2;

/// ...and each other end either is reached by at most half as many.
const BRIDGE_DOMINANCE: usize = 2;

/// One read alone joins two ends (a dip at 12-15x leaves a single read across, where the
/// chromosome broke) only if no other read comes into either end and the reads that run on
/// past either end without crossing agree with the bases the join puts there: at least this
/// many of them, none against, and between them reaching across the gap but for at most
/// [`MAX_GAP_UNVOUCHED`] bases. A read alone once joined across 4 kb it lacked (an *E. hormaechei*
/// chromosome came out 4 kb short); three reads agreed with its first kilobases, but none
/// reached the 4 kb.
const MIN_AGREEING: usize = 1;

/// At most this many of a lone crossing read's gap bases lie beyond the agreeing reads from
/// either side (whose ends fray): a deletion there could be no longer.
const MAX_GAP_UNVOUCHED: usize = 200;

/// A read running on past an end is weighed against a lone crossing read only if at least
/// this many of its bases lie beyond the end: fewer hardly tell a gap's bases from another's.
const MIN_CHECK_OUTWARD: usize = 500;

/// Such a read agrees if its outward bases align from their first to their last (or to the
/// last of the join's) within this many bases.
const CHECK_SLACK: u32 = 200;

/// A gap the reads cross holds at most this many bases.
const MAX_BRIDGE_GAP: usize = 10_000;

/// The bases between two ends are the consensus of at most this many reads'.
const MAX_BRIDGE_READS: usize = 20;

/// A read's gap is taken only if it is within this share (one over it) of the reads' median
/// length...
const GAP_SLACK_DIVISOR: i64 = 5;

/// ...or within this many bases of it.
const MIN_GAP_SLACK: i64 = 200;

/// What one read has between two ends: the gap's bases, or how far the ends overlap.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Crossing {
    gap: Vec<u8>,
    overlap: u32,
    /// The read that crosses.
    read: ReadIndex,
}

impl Crossing {
    /// The gap's length, or minus the overlap.
    fn signed_len(&self) -> i64 {
        self.gap.len() as i64 - i64::from(self.overlap)
    }
}

/// Joins linear drafts' ends that reads run between, through those reads' consensus bases;
/// returns the drafts after joining and how many joins were made.
///
/// Where reads are scarce over a stretch (a few reads at 20x), the string graph ends its contigs
/// either side, and their ends do not overlap: the reads across were left out as a shallow
/// contig's. Here the reads are mapped onto the last [`BRIDGE_PIECE`] polished bases at each
/// end of each linear draft (past its [`Draft::shallow_ends`]), onto the starts and the ends
/// apart, so that a read the length of a short draft comes off both its ends; those that align
/// up to the end and run on past it ([`MIN_OUTWARD`]) are mapped again with only those outward
/// bases, and where these come into another end, from its outer side, the bases before it are
/// the gap. A read found crossing from both ends counts once. Two ends at least
/// [`MIN_BRIDGE_READS`] reads cross between (counting those whose gap is
/// within [`GAP_SLACK_DIVISOR`] of the median), and [`BRIDGE_DOMINANCE`] times as many as to
/// any other end, are joined through the consensus of those bases. Drafts far shallower than
/// the longest ([`MIN_JOIN_DEPTH_SHARE`]) join no other.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn join_by_reads(
    drafts: Vec<Draft>,
    reads: &Reads,
    map: &OverlapParams,
    backend: &PoaBackend,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let linear: Vec<usize> = (0..drafts.len()).filter(|&i| !drafts[i].circular).collect();
    if linear.is_empty() {
        return Ok((drafts, 0));
    }
    // Pieces 2k and 2k + 1 are the start and the end of linear draft k's polished bases.
    let mut pieces = Reads::new();
    for &i in &linear {
        let deep = drafts[i].deep_part();
        let n = (deep.len() / 2).min(BRIDGE_PIECE);
        pieces.push(&deep[..n]);
        pieces.push(&deep[deep.len() - n..]);
    }
    let end_of = |piece: ReadIndex| -> End {
        let side = if piece.index() % 2 == 1 {
            Side::End
        } else {
            Side::Start
        };
        (linear[piece.index() / 2], side)
    };
    // The bases of each read past an end it aligns up to, read outward from the end. The reads
    // are mapped onto the starts and the ends apart: one that runs the length of a short draft
    // comes off both its ends.
    let mut outward: Vec<(ReadIndex, End, Vec<u8>)> = Vec::new();
    for side in [Side::Start, Side::End] {
        let first = usize::from(side == Side::End);
        let sided: Reads = (0..linear.len())
            .map(|k| pieces.get(ReadIndex(piece_index(2 * k + first))))
            .collect();
        for (read, m) in reads.indices().zip(map_reads(&sided, reads, map)?) {
            let Some(m) = m else { continue };
            let end = (linear[m.target.index()], side);
            let piece_len = sided.read_len(m.target);
            let len = reads.read_len(read);
            let bases = reads.get(read);
            // The read in the draft's orientation, and its aligned stretch there.
            let (oriented, from, to) = if m.strand.is_reverse() {
                (
                    reverse_complement(bases),
                    len - m.query_end,
                    len - m.query_start,
                )
            } else {
                (bases.to_vec(), m.query_start, m.query_end)
            };
            let out = match side {
                Side::End if m.target_end + END_SLACK >= piece_len && len - to >= MIN_OUTWARD => {
                    // From the base past the draft's last one.
                    let past = to + (piece_len - m.target_end);
                    oriented[(past.min(len)) as usize..].to_vec()
                }
                Side::Start if m.target_start <= END_SLACK && from >= MIN_OUTWARD => {
                    let before = from.saturating_sub(m.target_start);
                    reverse_complement(&oriented[..before as usize])
                }
                _ => continue,
            };
            outward.push((read, end, out));
        }
    }
    // A draft far shallower than the longest joins no other (another sample's reads, or a
    // chimera's, laid out as a draft of their own).
    let longest = (0..drafts.len()).max_by_key(|&i| (drafts[i].seq.len(), std::cmp::Reverse(i)));
    let least = longest.map_or(0.0, |i| MIN_JOIN_DEPTH_SHARE * drafts[i].depth);
    let deep = |d: usize| drafts[d].depth >= least;
    let queries: Reads = outward
        .iter()
        .map(|(_, _, o)| &o[..o.len().min(MAX_BRIDGE_GAP + BRIDGE_PIECE)])
        .collect();
    // Each read that comes into another end: the gap's bases (or how far the ends overlap),
    // keyed by the two ends in order and read out of the first.
    let mut crossings: BTreeMap<(End, End), Vec<Crossing>> = BTreeMap::new();
    // A read that crosses between two ends is found from both: it counts once.
    let mut counted: BTreeSet<(End, End, ReadIndex)> = BTreeSet::new();
    for ((read, from, out), m) in outward.iter().zip(map_reads(&pieces, &queries, map)?) {
        let Some(m) = m else { continue };
        let to = end_of(m.target);
        let piece_len = pieces.read_len(m.target);
        // Into a start from its first base on, or into an end from its last base back.
        let short_of_end = match (to.1, m.strand.is_reverse()) {
            (Side::Start, false) => m.target_start,
            (Side::End, true) => piece_len - m.target_end,
            _ => continue,
        };
        if to == *from
            || short_of_end > END_SLACK
            || (from.0 != to.0 && !(deep(from.0) && deep(to.0)))
        {
            continue;
        }
        // The read's bases up to where the other end's outermost base lies in it; if they
        // start inside that end, the two ends overlap by the difference.
        let crossing = if m.query_start >= short_of_end {
            let gap = &out[..(m.query_start - short_of_end) as usize];
            if gap.len() > MAX_BRIDGE_GAP {
                continue;
            }
            Crossing {
                gap: gap.to_vec(),
                overlap: 0,
                read: *read,
            }
        } else {
            Crossing {
                gap: Vec::new(),
                overlap: short_of_end - m.query_start,
                read: *read,
            }
        };
        if !counted.insert((*from.min(&to), *from.max(&to), *read)) {
            continue;
        }
        if *from < to {
            crossings.entry((*from, to)).or_default().push(crossing);
        } else {
            crossings.entry((to, *from)).or_default().push(Crossing {
                gap: reverse_complement(&crossing.gap),
                ..crossing
            });
        }
    }
    // Only the reads whose gap agrees with the others': a read whose outward bases land on a
    // repeat's copy at the far end gives another length.
    for found in crossings.values_mut() {
        let mut lens: Vec<i64> = found.iter().map(Crossing::signed_len).collect();
        lens.sort_unstable();
        let median = lens[lens.len() / 2];
        let slack = (median.abs() / GAP_SLACK_DIVISOR).max(MIN_GAP_SLACK);
        found.retain(|c| (c.signed_len() - median).abs() <= slack);
    }
    // The bases of each end's piece read inward from its outer end.
    let inward = |(d, side): End| -> Vec<u8> {
        let k = linear
            .iter()
            .position(|&i| i == d)
            .expect("ends are of linear drafts");
        match side {
            Side::Start => pieces.get(ReadIndex(piece_index(2 * k))).to_vec(),
            Side::End => reverse_complement(pieces.get(ReadIndex(piece_index(2 * k + 1)))),
        }
    };
    let rivals = |a: End, b: End| {
        crossings
            .iter()
            .filter(move |&(&(c, d), _)| {
                (c, d) != (a, b) && [c, d].iter().any(|&e| e == a || e == b)
            })
            .map(|(_, others)| others.len())
    };
    // The ends one read alone crosses between that pass the check of [`MIN_AGREEING`].
    let mut agreed_alone: BTreeSet<(End, End)> = BTreeSet::new();
    for (&(a, b), found) in &crossings {
        if let [crossing] = found.as_slice()
            && rivals(a, b).sum::<usize>() == 0
            && agreed(a, b, crossing, &outward, &inward, map)?
        {
            agreed_alone.insert((a, b));
        }
    }
    let joins: Vec<Join> = crossings
        .iter()
        .filter(|&(&(a, b), found)| {
            (found.len() >= MIN_BRIDGE_READS
                && rivals(a, b).all(|n| n * BRIDGE_DOMINANCE <= found.len()))
                || agreed_alone.contains(&(a, b))
        })
        .map(|(&(from, to), found)| {
            log::debug!(
                "bacteria: ends {from:?} and {to:?} crossed by reads with gaps {:?}",
                found.iter().map(Crossing::signed_len).collect::<Vec<_>>()
            );
            let mut overlaps: Vec<u32> = found.iter().map(|c| c.overlap).collect();
            overlaps.sort_unstable();
            let to_overlap = overlaps[overlaps.len() / 2];
            let gaps: Vec<&[u8]> = found
                .iter()
                .filter(|c| c.overlap == 0)
                .take(MAX_BRIDGE_READS)
                .map(|c| c.gap.as_slice())
                .collect();
            let between = if to_overlap > 0 {
                Vec::new()
            } else {
                backend.consensus_segments(&gaps)
            };
            Join {
                from,
                to,
                to_overlap,
                between,
            }
        })
        .collect();
    if joins.is_empty() {
        return Ok((drafts, 0));
    }
    let mut drafts = drafts;
    drop_shallow_ends(&mut drafts, &joins);
    let made = joins.len();
    Ok((merge(drafts, &joins), made))
}

/// Whether the reads in `outward` that run on past end `a` or `b` (other than `crossing`'s,
/// the one read across between them, its gap read out of `a`) agree with the bases the join
/// puts past each end: see [`MIN_AGREEING`]. `inward` gives an end's bases read inward.
fn agreed(
    a: End,
    b: End,
    crossing: &Crossing,
    outward: &[(ReadIndex, End, Vec<u8>)],
    inward: &impl Fn(End) -> Vec<u8>,
    map: &OverlapParams,
) -> Result<bool, InvalidOverlapParamsError> {
    let overlap = crossing.overlap as usize;
    let (mut agree, mut against) = (0, 0);
    // How far past each end the agreeing reads reach.
    let mut reach = [0_usize; 2];
    for (side, (end, other, gap)) in [
        (a, b, crossing.gap.clone()),
        (b, a, reverse_complement(&crossing.gap)),
    ]
    .into_iter()
    .enumerate()
    {
        // The bases past `end` once joined: the gap, then the other end's read inward.
        let mut past = gap;
        let into = inward(other);
        past.extend_from_slice(&into[overlap.min(into.len())..]);
        let queries: Reads = outward
            .iter()
            .filter(|(read, e, out)| {
                *e == end && *read != crossing.read && out.len() >= MIN_CHECK_OUTWARD
            })
            .map(|(_, _, out)| &out[..out.len().min(past.len())])
            .collect();
        if queries.is_empty() {
            continue;
        }
        let target: Reads = std::iter::once(past.as_slice()).collect();
        for (q, m) in queries.indices().zip(map_reads(&target, &queries, map)?) {
            let q_len = queries.read_len(q);
            let fits = m.as_ref().is_some_and(|m| {
                !m.strand.is_reverse()
                    && m.query_start <= CHECK_SLACK
                    && m.target_start <= CHECK_SLACK
                    && (q_len - m.query_end <= CHECK_SLACK
                        || target.read_len(m.target) - m.target_end <= CHECK_SLACK)
            });
            if let Some(m) = m.filter(|_| fits) {
                agree += 1;
                reach[side] = reach[side].max(m.target_end as usize);
            } else {
                against += 1;
            }
        }
    }
    let gap = crossing.gap.len();
    log::debug!(
        "bacteria: ends {a:?} and {b:?} crossed by one read (gap {}): {agree} reads past them \
         agree, reaching {reach:?}, {against} against",
        crossing.signed_len()
    );
    Ok(agree >= MIN_AGREEING && against == 0 && reach[0] + reach[1] + MAX_GAP_UNVOUCHED >= gap)
}

/// Returns the joins whose overlap is at least [`JOIN_DOMINANCE`] times as long as that of every
/// rival (a join that shares an end with it); two of them never share an end, as each would
/// have to be the longer.
fn dominant(joins: &[Join]) -> Vec<Join> {
    let (num, den) = JOIN_DOMINANCE;
    joins
        .iter()
        .enumerate()
        .filter(|&(i, j)| {
            joins.iter().enumerate().all(|(k, r)| {
                let rival = k != i && [r.from, r.to].iter().any(|&e| e == j.from || e == j.to);
                !rival
                    || u64::from(j.to_overlap) * u64::from(den)
                        >= u64::from(r.to_overlap) * u64::from(num)
            })
        })
        .map(|(_, j)| j.clone())
        .collect()
}

/// Returns the join an overlap between two pieces makes, if it is a dovetail between their
/// outer ends.
fn join_of(o: &Overlap, pieces: &Reads, linear: &[usize]) -> Option<Join> {
    let q_len = pieces.read_len(o.query);
    let t_len = pieces.read_len(o.target);
    let (q_start, q_end) = if o.strand.is_reverse() {
        (q_len - o.query_end, q_len - o.query_start)
    } else {
        (o.query_start, o.query_end)
    };
    let end_of = |piece: ReadIndex, side_of_piece: Side| -> Option<End> {
        let draft = linear[piece.index() / 2];
        let is_end_piece = piece.index() % 2 == 1;
        // A start piece's outer end is its first base; an end piece's, its last.
        match (is_end_piece, side_of_piece) {
            (false, Side::Start) => Some((draft, Side::Start)),
            (true, Side::End) => Some((draft, Side::End)),
            _ => None,
        }
    };
    let reverse = o.strand.is_reverse();
    // In the target's orientation: target then query, or query then target.
    if o.target_start > 0 && o.target_end == t_len && q_start == 0 && q_end < q_len {
        // The target's last base runs into the (oriented) query's first.
        let from = end_of(o.target, Side::End)?;
        let to = end_of(o.query, if reverse { Side::End } else { Side::Start })?;
        Some(Join {
            from,
            to,
            to_overlap: q_end,
            between: Vec::new(),
        })
    } else if q_start > 0 && q_end == q_len && o.target_start == 0 && o.target_end < t_len {
        let from = end_of(o.query, if reverse { Side::Start } else { Side::End })?;
        let to = end_of(o.target, Side::Start)?;
        Some(Join {
            from,
            to,
            to_overlap: o.target_end,
            between: Vec::new(),
        })
    } else {
        None
    }
}

/// Merges `drafts` along `joins` (each end in at most one).
fn merge(drafts: Vec<Draft>, joins: &[Join]) -> Vec<Draft> {
    // The end a join leads into from `end`, how many of its bases the overlap covers, and the
    // bases before it.
    let partner = |end: End| -> Option<(End, u32, Vec<u8>)> {
        joins.iter().find_map(|j| {
            if j.from == end {
                Some((j.to, j.to_overlap, j.between.clone()))
            } else if j.to == end {
                // The same join read the other way: `from`'s side of it. The overlap's length
                // on `from`'s draft is taken to be the same.
                Some((j.from, j.to_overlap, reverse_complement(&j.between)))
            } else {
                None
            }
        })
    };
    let mut used = vec![false; drafts.len()];
    let mut out = Vec::new();
    for first in 0..drafts.len() {
        if used[first] || drafts[first].circular {
            continue;
        }
        // Walk back to the start of the chain (or round a circle back to `first`).
        let mut start = (first, Side::Start);
        let mut circle = false;
        while let Some(((prev, prev_side), _, _)) = partner(start) {
            if prev == first {
                circle = true;
                break;
            }
            // Entering `prev` at `prev_side` means it is read towards that end.
            start = (prev, opposite(prev_side));
        }
        if circle {
            start = (first, Side::Start);
        }
        // Walk forward, reading each draft from `start`'s side to the other.
        let mut seq: Vec<u8> = Vec::new();
        let mut at = start;
        let mut skip = 0_usize;
        // The chain's depth: its drafts', weighted by their lengths.
        let (mut weighted, mut bases) = (0.0, 0.0);
        loop {
            used[at.0] = true;
            #[expect(
                clippy::cast_precision_loss,
                reason = "a mean depth, not an exact count"
            )]
            let len = drafts[at.0].seq.len() as f64;
            weighted += drafts[at.0].depth * len;
            bases += len;
            let oriented = match at.1 {
                Side::Start => drafts[at.0].seq.clone(),
                Side::End => reverse_complement(&drafts[at.0].seq),
            };
            seq.extend_from_slice(&oriented[skip.min(oriented.len())..]);
            let leaving = (at.0, opposite(at.1));
            match partner(leaving) {
                Some((next, overlap, between)) if !used[next.0] => {
                    seq.extend_from_slice(&between);
                    skip = usize::try_from(overlap).expect("fits");
                    at = next;
                }
                Some((next, overlap, between)) if circle && next.0 == start.0 => {
                    // Close the circle: the start's first bases are already at the end.
                    seq.extend_from_slice(&between);
                    let overlap = usize::try_from(overlap).expect("fits");
                    seq.drain(..overlap.min(seq.len()));
                    break;
                }
                _ => break,
            }
        }
        let mut merged = Draft::new(seq, circle);
        merged.depth = if bases > 0.0 { weighted / bases } else { 0.0 };
        if !circle {
            // The chain's outer ends are its first draft's leading end and its last's trailing
            // one, each read from the side it was entered at.
            let side = |(d, entered): End, leading: bool| {
                let first_end = (entered == Side::Start) == leading;
                drafts[d].shallow_ends[usize::from(!first_end)]
            };
            merged.shallow_ends = [side(start, true), side(at, false)];
        }
        out.push(merged);
    }
    out.extend(drafts.into_iter().filter(|d| d.circular));
    out
}

/// A linear draft shorter than this may be one replicon read round more than once.
const MAX_TANDEM_LEN: usize = 100_000;

/// The period of a tandem draft is voted on by the k-mers of this many bases...
const TANDEM_K: usize = 16;

/// ...that start in this many bases in its middle (or a quarter of it, if shorter)...
const TANDEM_PROBE_SPAN: usize = 2_000;

/// ...and at least this share of them (in percent) must find each other at the period.
const TANDEM_SUPPORT_PERCENT: usize = 30;

/// The period's exact length is taken from the k-mers in the first this many bases of the
/// probes.
const TANDEM_EXACT_SPAN: usize = 500;

/// Collapses the short linear `drafts` that hold a small replicon two or more times over (see
/// [`collapse_tandem`]), keeping their depth; returns the drafts and how many were collapsed. Run
/// again once the drafts are polished: on the first round's drafts the copies of a 9.2 kb
/// *A. tenax* plasmid differed too much for their k-mers to agree on its length, and it came out
/// as a 20 kb line of 2.2 copies.
pub(crate) fn collapse_tandems(drafts: Vec<Draft>) -> (Vec<Draft>, usize) {
    let mut collapsed = 0;
    let drafts = drafts
        .into_iter()
        .map(|draft| {
            let (depth, was_circular) = (draft.depth, draft.circular);
            let mut done = collapse_tandem(draft);
            if done.circular && !was_circular {
                collapsed += 1;
                done.depth = depth;
            }
            done
        })
        .collect();
    (drafts, collapsed)
}

/// Returns `draft` as one circle of its unit if it is a short linear draft that holds a small
/// replicon two or more times over (a small plasmid's reads, each round the circle more than
/// once, laid out end to end), else as it is. The unit's length is the commonest distance from
/// a k-mer in the draft's middle to the same k-mer elsewhere, if enough of them agree; the
/// copies need not be identical (the draft is polished once by now, each copy from a share of
/// the reads, and worst at its ends).
fn collapse_tandem(draft: Draft) -> Draft {
    // Its raw read bases past the polished ones are no copy of the unit.
    let seq = draft.deep_part();
    let len = seq.len();
    if draft.circular || !(4 * TANDEM_K..=MAX_TANDEM_LEN).contains(&len) {
        return draft;
    }
    // Probes from the middle of the draft: its ends are where it is polished worst.
    let span = TANDEM_PROBE_SPAN.min(len / 4);
    let mid = len / 2 - span / 2;
    let mut at: std::collections::HashMap<&[u8], Vec<usize>> = std::collections::HashMap::new();
    for (i, k) in seq.windows(TANDEM_K).enumerate() {
        at.entry(k).or_default().push(i);
    }
    let distances = |i: usize| {
        at[&seq[i..i + TANDEM_K]]
            .iter()
            .map(move |&j| j.abs_diff(i))
            .filter(|&d| d > 2 * TANDEM_K)
    };
    // Distances binned to 8 bases, so a copy with an indel still votes with the rest.
    let mut votes: BTreeMap<usize, usize> = BTreeMap::new();
    for i in mid..mid + span {
        let bins: BTreeSet<usize> = distances(i).map(|d| d / 8).collect();
        for bin in bins {
            *votes.entry(bin).or_default() += 1;
        }
    }
    // The shortest well-supported distance (multiples of the unit are supported too).
    let needed = (span * TANDEM_SUPPORT_PERCENT).div_ceil(100);
    let Some(bin) = votes
        .iter()
        .find(|&(&bin, &n)| n + votes.get(&(bin + 1)).copied().unwrap_or(0) >= needed)
        .map(|(&bin, _)| bin)
    else {
        return draft;
    };
    // The exact period: the commonest distance within the bin and the next, from the k-mers
    // at the start of the probes (where the unit is cut; a copy further on may have lost or
    // gained a base).
    let mut exact: BTreeMap<usize, usize> = BTreeMap::new();
    for i in mid..mid + span.min(TANDEM_EXACT_SPAN) {
        for d in distances(i).filter(|d| (bin * 8..(bin + 2) * 8).contains(d)) {
            *exact.entry(d).or_default() += 1;
        }
    }
    let Some((&period, _)) = exact
        .iter()
        .max_by_key(|&(&d, &n)| (n, std::cmp::Reverse(d)))
    else {
        return draft;
    };
    if 2 * period > len + len / 10 || !repeats_throughout(seq, period, &at) {
        return draft;
    }
    // The unit from the middle, forward or back as the draft allows.
    let start = if mid + period <= len {
        mid
    } else {
        mid.saturating_sub(period)
    };
    log::debug!("bacteria: a linear draft of {len} bases is a circle of {period} read round");
    Draft::new(seq[start..start + period].to_vec(), true)
}

/// Of the k-mers sampled along a tandem draft, at least this share (in percent) are found
/// again one period on or back: a chromosome fragment whose middle lies in a repeat (an IS
/// element with another copy nearby) agrees with a period only there.
const TANDEM_THROUGHOUT_PERCENT: usize = 50;

/// Returns whether k-mers all along `seq` (one every 64 bases) recur `period` bases on or back
/// (give or take 8, an indel a copy keeps), as they do in a draft that is one unit read round
/// and round.
fn repeats_throughout(
    seq: &[u8],
    period: usize,
    at: &std::collections::HashMap<&[u8], Vec<usize>>,
) -> bool {
    let samples: Vec<usize> = (0..seq.len() - TANDEM_K).step_by(64).collect();
    let recurring = samples
        .iter()
        .filter(|&&i| {
            at[&seq[i..i + TANDEM_K]]
                .iter()
                .any(|&j| j.abs_diff(i).abs_diff(period) <= 8)
        })
        .count();
    recurring * 100 >= samples.len() * TANDEM_THROUGHOUT_PERCENT
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

fn piece_index(i: usize) -> u32 {
    u32::try_from(i).expect("fewer than 4 billion draft ends")
}

fn opposite(side: Side) -> Side {
    match side {
        Side::Start => Side::End,
        Side::End => Side::Start,
    }
}

#[cfg(test)]
mod tests {
    use crate::overlap::Edit;

    use super::*;

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

    #[test]
    fn reads_across_gaps_join_the_drafts_either_side() {
        // A circle of 60 kb laid out as two drafts with an 800-base gap after each, which only
        // three reads (each way round) cross.
        let genome = genome(60_000, 11);
        let drafts = vec![
            Draft::new(genome[1_000..29_200].to_vec(), false),
            Draft::new(genome[30_000..59_200].to_vec(), false),
        ];
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for at in [24_000, 25_000, 26_000, 54_000, 55_000, 56_000] {
            let read: Vec<u8> = genome
                .iter()
                .cycle()
                .skip(at)
                .take(9_000)
                .copied()
                .collect();
            seqs.push(if at % 2_000 == 0 {
                read
            } else {
                reverse_complement(&read)
            });
        }
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let backend = PoaBackend::new(crate::consensus::PoaParams::default()).expect("defaults");
        let (drafts, made) = join_by_reads(drafts, &reads, &map, &backend).expect("valid");
        assert_eq!(made, 2);
        assert_eq!(drafts.len(), 1);
        assert!(drafts[0].circular);
        assert!(is_rotation(&drafts[0].seq, &genome));
    }

    /// The reads of `genome` (a circle) from each `(start, len)`, every other one reversed.
    fn reads_at(genome: &[u8], at: &[(usize, usize)]) -> Vec<Vec<u8>> {
        at.iter()
            .enumerate()
            .map(|(i, &(start, len))| {
                let read: Vec<u8> = genome
                    .iter()
                    .cycle()
                    .skip(start)
                    .take(len)
                    .copied()
                    .collect();
                if i % 2 == 0 {
                    read
                } else {
                    reverse_complement(&read)
                }
            })
            .collect()
    }

    #[test]
    fn one_read_across_a_gap_joins_where_the_reads_running_into_it_agree() {
        // A circle of 60 kb laid out as two drafts with 800-base gaps between: three reads
        // cross the second, one alone the first, and one read from either side runs 700
        // bases into that gap and stops.
        let genome = genome(60_000, 17);
        let drafts = vec![
            Draft::new(genome[1_000..29_200].to_vec(), false),
            Draft::new(genome[30_000..59_200].to_vec(), false),
        ];
        let seqs = reads_at(
            &genome,
            &[
                (25_000, 9_000),
                (24_000, 5_900),
                (29_300, 5_000),
                (54_000, 9_000),
                (55_000, 9_000),
                (56_000, 9_000),
            ],
        );
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let backend = PoaBackend::new(crate::consensus::PoaParams::default()).expect("defaults");
        let (drafts, made) = join_by_reads(drafts, &reads, &map, &backend).expect("valid");
        assert_eq!((made, drafts.len()), (2, 1));
        assert!(drafts[0].circular);
        assert!(is_rotation(&drafts[0].seq, &genome));
    }

    #[test]
    fn one_read_lacking_a_gap_the_other_reads_run_into_joins_nothing() {
        // As above, but the first gap is 4.8 kb, and the one read across lacks 3 kb of it (a
        // deletion); a read running 600 bases into the gap agrees with the read's first
        // bases there, but nothing vouches for the rest, so only the second gap joins.
        let genome = genome(60_000, 19);
        let drafts = vec![
            Draft::new(genome[1_000..29_200].to_vec(), false),
            Draft::new(genome[34_000..59_200].to_vec(), false),
        ];
        let mut seqs = reads_at(
            &genome,
            &[
                (26_000, 3_800),
                (54_000, 9_000),
                (55_000, 9_000),
                (56_000, 9_000),
            ],
        );
        let mut lacking = genome[25_000..30_000].to_vec();
        lacking.extend_from_slice(&genome[33_000..39_000]);
        seqs.push(lacking);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let backend = PoaBackend::new(crate::consensus::PoaParams::default()).expect("defaults");
        let (drafts, made) = join_by_reads(drafts, &reads, &map, &backend).expect("valid");
        assert_eq!((made, drafts.len()), (1, 1));
        assert!(!drafts[0].circular);
        assert_eq!(drafts[0].seq.len(), 60_000 - 4_800);
    }

    #[test]
    fn reads_the_length_of_a_short_draft_join_both_its_ends() {
        // A circle of 60 kb laid out as a long draft and a short one, with gaps between, which
        // only two reads cross; each runs the length of the short draft, across both gaps.
        let genome = genome(60_000, 13);
        let drafts = vec![
            Draft::new(genome[2_000..50_000].to_vec(), false),
            Draft::new(genome[51_000..58_000].to_vec(), false),
        ];
        let seqs: Vec<Vec<u8>> = [47_000, 47_500]
            .iter()
            .map(|&at| {
                let read: Vec<u8> = genome
                    .iter()
                    .cycle()
                    .skip(at)
                    .take(17_000)
                    .copied()
                    .collect();
                if at % 1_000 == 0 {
                    read
                } else {
                    reverse_complement(&read)
                }
            })
            .collect();
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let backend = PoaBackend::new(crate::consensus::PoaParams::default()).expect("defaults");
        let (drafts, made) = join_by_reads(drafts, &reads, &map, &backend).expect("valid");
        assert_eq!((made, drafts.len()), (2, 1));
        assert!(drafts[0].circular);
        assert!(is_rotation(&drafts[0].seq, &genome));
    }

    #[test]
    fn reads_across_a_gap_close_a_draft_or_join_a_reversed_one() {
        let genome = genome(60_000, 12);
        let reads_at = |starts: &[usize]| -> Vec<Vec<u8>> {
            starts
                .iter()
                .map(|&at| {
                    genome
                        .iter()
                        .cycle()
                        .skip(at)
                        .take(9_000)
                        .copied()
                        .collect()
                })
                .collect()
        };
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let backend = PoaBackend::new(crate::consensus::PoaParams::default()).expect("defaults");
        // One draft lacking 800 bases where its end would run into its start.
        let seqs = reads_at(&[54_000, 55_000, 56_000]);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let drafts = vec![Draft::new(genome[..59_200].to_vec(), false)];
        let (drafts, made) = join_by_reads(drafts, &reads, &map, &backend).expect("valid");
        assert_eq!((made, drafts.len()), (1, 1));
        assert!(drafts[0].circular);
        assert!(is_rotation(&drafts[0].seq, &genome));
        // Two drafts, the second turned round, so the gaps join end to end and start to start.
        let seqs = reads_at(&[24_000, 25_000, 26_000, 54_000, 55_000, 56_000]);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let drafts = vec![
            Draft::new(genome[1_000..29_200].to_vec(), false),
            Draft::new(reverse_complement(&genome[30_000..59_200]), false),
        ];
        let (drafts, made) = join_by_reads(drafts, &reads, &map, &backend).expect("valid");
        assert_eq!((made, drafts.len()), (2, 1));
        assert!(drafts[0].circular);
        let seq = &drafts[0].seq;
        assert!(is_rotation(seq, &genome) || is_rotation(&reverse_complement(seq), &genome));
    }

    /// A circle of three 40 kb stretches with an 8 kb repeat after each, laid out as three
    /// drafts that each end in copies of the repeat (every end overlaps the other two drafts'
    /// along all 8 kb), and reads that span the repeat at the junctions `spanned` (by their
    /// first base in the genome) plus reads of 6 kb inside every copy.
    fn three_pieces_ending_in_one_repeat(spanned: &[usize]) -> (Vec<u8>, Vec<Draft>, Vec<Vec<u8>>) {
        let repeat = genome(8_000, 41);
        let unique: Vec<Vec<u8>> = (0..3).map(|i| genome(40_000, 42 + i)).collect();
        let genome: Vec<u8> = unique
            .iter()
            .flat_map(|u| [u.as_slice(), &repeat].concat())
            .collect();
        let at = |from: usize, to: usize| -> Vec<u8> {
            genome
                .iter()
                .cycle()
                .skip(from)
                .take(to - from)
                .copied()
                .collect()
        };
        // Repeat copies at 40-48 kb, 88-96 kb and 136-144 kb.
        let drafts = vec![
            Draft::new(at(136_000, 192_000), false),
            Draft::new(at(40_000, 96_000), false),
            Draft::new(at(88_000, 144_000), false),
        ];
        let mut placed: Vec<(usize, usize)> = spanned
            .iter()
            .flat_map(|&start| [(start - 1_500, 11_000), (start - 1_000, 10_500)])
            .collect();
        placed.extend([(41_000, 6_000), (89_000, 6_000), (137_000, 6_000)]);
        let reads = reads_at(&genome, &placed);
        (genome, drafts, reads)
    }

    #[test]
    fn ends_in_copies_of_one_repeat_join_where_reads_span_it() {
        // Reads span the first two copies only; the third junction is then the one left.
        let (genome, drafts, seqs) = three_pieces_ending_in_one_repeat(&[40_000, 88_000]);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let params = OverlapParams {
            extension: crate::overlap::Extension::Full,
            ..map
        };
        // Without the reads the ends are tied.
        let (drafts, made) = join_ends(drafts, &params).expect("valid");
        assert_eq!((made, drafts.len()), (0, 3));
        let (drafts, made) = join_across_repeats(drafts, &reads, &params, &map).expect("valid");
        assert_eq!((made, drafts.len()), (2, 1));
        let (drafts, _) = join_more(drafts, &params).expect("valid");
        assert!(drafts[0].circular);
        let seq = &drafts[0].seq;
        assert!(is_rotation(seq, &genome) || is_rotation(&reverse_complement(seq), &genome));
    }

    #[test]
    fn ends_in_copies_of_one_repeat_no_read_spans_stay_apart() {
        let (_, drafts, seqs) = three_pieces_ending_in_one_repeat(&[]);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let params = OverlapParams {
            extension: crate::overlap::Extension::Full,
            ..map
        };
        let (drafts, made) = join_across_repeats(drafts, &reads, &params, &map).expect("valid");
        assert_eq!((made, drafts.len()), (0, 3));
    }

    #[test]
    fn ends_that_overlap_past_a_raw_end_join_without_it() {
        // A circle of 60 kb laid out round past its start by 13 kb, then 1.3 kb of a read's
        // unpolished bases that match nothing: too long an overhang, until left out.
        let circle = genome(60_000, 15);
        let junk = genome(1_300, 16);
        let seq = [circle.as_slice(), &circle[..13_000], &junk].concat();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        for (shallow, closed) in [(0, false), (1_300, true)] {
            let mut draft = Draft::new(seq.clone(), false);
            draft.shallow_ends = [0, shallow];
            let (drafts, _) = join_more(vec![draft], &map).expect("valid");
            assert_eq!(drafts[0].circular, closed, "{shallow} shallow bases");
            if closed {
                assert!(is_rotation(&drafts[0].seq, &circle));
            }
        }
    }

    #[test]
    fn ends_that_overlap_with_an_element_in_one_close_through_anchors() {
        // A circle of 200 kb laid out round past its start by 40 kb, the end's copy of those
        // 40 kb with an 8 kb element inside that the start's lacks: no dovetail, but the
        // anchors close it, with the element.
        let circle = genome(200_000, 17);
        let element = genome(8_000, 18);
        let seq = [
            &circle[..],
            &circle[..20_000],
            &element,
            &circle[20_000..40_000],
        ]
        .concat();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let (drafts, made) = close_by_anchors(vec![Draft::new(seq, false)], &map).expect("valid");
        assert_eq!(made, 1);
        assert!(drafts[0].circular);
        let with = [&circle[..20_000], &element, &circle[20_000..]].concat();
        assert!(is_rotation(&drafts[0].seq, &with));
        // A draft whose ends do not meet is left linear.
        let (drafts, made) =
            close_by_anchors(vec![Draft::new(circle[..150_000].to_vec(), false)], &map)
                .expect("valid");
        assert_eq!((made, drafts[0].circular), (0, false));
    }

    #[test]
    fn bases_a_chimeric_read_put_on_the_ends_are_left_out_when_closing() {
        // A 200 kb circle in two pieces whose ends overlap by 15 kb, with 8 kb of another
        // replicon on the first piece's start and 3 kb on the second's end: no dovetail,
        // but the anchors join the pieces and close the circle without the other bases.
        let circle = genome(200_000, 23);
        let (junk, more) = (genome(8_000, 24), genome(3_000, 25));
        let first = [&junk[..], &circle[..115_000]].concat();
        let second = [&circle[100_000..], &circle[..15_000], &more[..]].concat();
        let map = crate::BactParams::new(crate::ReadKind::Ont).map;
        let mut drafts = vec![Draft::new(first, false), Draft::new(second, false)];
        for d in &mut drafts {
            d.depth = 50.0;
        }
        let (drafts, joined) = join_by_anchors(drafts, &map).expect("valid");
        assert_eq!((joined, drafts.len()), (1, 1));
        let (drafts, closed) = close_by_anchors(drafts, &map).expect("valid");
        assert_eq!(closed, 1);
        assert!(is_rotation(&drafts[0].seq, &circle));
    }

    #[test]
    fn an_end_with_rivals_is_joined_by_a_clearly_longer_overlap() {
        let join = |from: End, to: End, to_overlap: u32| Join {
            from,
            to,
            to_overlap,
            between: Vec::new(),
        };
        // Draft 2's end overlaps draft 0's end by 18 kb and draft 3's start by 12.5 kb (a
        // repeat); draft 3's start also overlaps draft 1's start by 21 kb.
        let joins = [
            join((0, Side::End), (2, Side::End), 18_000),
            join((2, Side::End), (3, Side::Start), 12_500),
            join((1, Side::Start), (3, Side::Start), 21_000),
            join((4, Side::End), (5, Side::Start), 9_000),
            join((4, Side::End), (6, Side::Start), 8_000),
        ];
        // The two long ones win; draft 4's end, with two near-equal rivals, is left alone.
        assert_eq!(dominant(&joins), [joins[0].clone(), joins[2].clone()]);
    }

    #[test]
    fn a_draft_of_a_plasmid_read_round_three_times_becomes_one_circle() {
        let unit = genome(3_514, 3);
        let mut seq = unit.repeat(2);
        seq.extend_from_slice(&unit[..3_000]);
        let d = collapse_tandem(Draft::new(seq, false));
        assert!(d.circular);
        assert!(is_rotation(&d.seq, &unit));
    }

    #[test]
    fn polished_tandem_lines_are_collapsed_and_keep_their_depth() {
        // A 9.2 kb plasmid read round 2.2 times, a line of other bases, and a circle.
        let unit = genome(9_244, 7);
        let mut tandem = unit.repeat(2);
        tandem.extend_from_slice(&unit[..1_800]);
        let mut drafts = vec![
            Draft::new(tandem, false),
            Draft::new(genome(30_000, 8), false),
            Draft::new(genome(40_000, 9), true),
        ];
        for d in &mut drafts {
            d.depth = 120.0;
        }
        let (drafts, collapsed) = collapse_tandems(drafts);
        assert_eq!(collapsed, 1);
        assert!(drafts[0].circular && is_rotation(&drafts[0].seq, &unit));
        assert!((drafts[0].depth - 120.0).abs() < f64::EPSILON);
        assert!(!drafts[1].circular && drafts[1].seq.len() == 30_000);
        assert!(drafts[2].circular && drafts[2].seq.len() == 40_000);
    }

    /// Whether `a` is `b` turned round.
    fn is_rotation(a: &[u8], b: &[u8]) -> bool {
        a.len() == b.len() && [b, b].concat().windows(a.len()).any(|w| w == a)
    }

    #[test]
    fn copies_that_differ_a_little_still_collapse() {
        let unit = genome(3_514, 5);
        let mut second = unit.clone();
        for p in (100..3_400).step_by(300) {
            second[p] = if second[p] == b'A' { b'C' } else { b'A' };
        }
        second.remove(1_000);
        let mut seq = unit.clone();
        seq.extend_from_slice(&second);
        seq.extend_from_slice(&unit[..2_500]);
        let d = collapse_tandem(Draft::new(seq, false));
        assert!(d.circular);
        assert!(
            d.seq.len().abs_diff(unit.len()) <= 1,
            "a unit of {} bases",
            d.seq.len()
        );
    }

    #[test]
    fn a_fragment_with_a_repeat_in_its_middle_is_left_alone() {
        // 90 kb whose middle 2.5 kb (an IS element) occurs again 40 kb on.
        let is = genome(2_500, 7);
        let mut seq = genome(43_750, 8);
        seq.extend_from_slice(&is);
        seq.extend(genome(37_500, 9));
        seq.extend_from_slice(&is);
        seq.extend(genome(3_750, 10));
        let d = collapse_tandem(Draft::new(seq.clone(), false));
        assert!(!d.circular);
        assert_eq!(d.seq, seq);
    }

    #[test]
    fn a_draft_without_a_period_is_left_alone() {
        let seq = genome(20_000, 4);
        let d = collapse_tandem(Draft::new(seq.clone(), false));
        assert!(!d.circular);
        assert_eq!(d.seq, seq);
    }

    #[test]
    fn only_polished_columns_are_held_to_the_divergence_cap_unless_too_few() {
        // Two 4 kb pieces alike for 1 kb, then differing at one column in ten.
        let seq = genome(4_000, 9);
        let pieces: Reads = [seq.as_slice(), seq.as_slice()].into_iter().collect();
        let mut edits = vec![Edit::new(EditOp::Match, 1_000)];
        for _ in 0..300 {
            edits.push(Edit::new(EditOp::Mismatch, 1));
            edits.push(Edit::new(EditOp::Match, 9));
        }
        let o = Overlap {
            query: ReadIndex(0),
            target: ReadIndex(1),
            strand: Strand::Forward,
            query_start: 0,
            query_end: 4_000,
            target_start: 0,
            target_end: 4_000,
            edits,
        };
        // The differing 3 kb a raw read's: the polished kilobase agrees.
        let raw_tail = [(0, 1_000), (0, 4_000)];
        assert_eq!(
            divergence(&o, &pieces, &raw_tail),
            Divergence::Polished(0.0)
        );
        // All polished: 300 substitutions in 4,000 columns.
        let polished = [(0, 4_000), (0, 4_000)];
        assert_eq!(
            divergence(&o, &pieces, &polished),
            Divergence::Polished(0.075)
        );
        assert_eq!(
            within_cap(std::slice::from_ref(&o), &pieces, &polished),
            [false]
        );
        // Too few polished columns to tell: all of them count, and an overlap over the cap is
        // kept only if no closer one competes at its pieces.
        let few = [(0, 500), (0, 4_000)];
        assert_eq!(divergence(&o, &pieces, &few), Divergence::Raw(0.075));
        assert_eq!(within_cap(std::slice::from_ref(&o), &pieces, &few), [true]);
        let mut closer = o.clone();
        closer.edits = vec![Edit::new(EditOp::Match, 4_000)];
        assert_eq!(within_cap(&[o, closer], &pieces, &few), [false, true]);
    }
}
