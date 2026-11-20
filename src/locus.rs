//! Loci where most reads have another version than the contig.
//!
//! A bacterial population can carry two versions of a locus: a segment between inverted repeats
//! that recombination turns round (phase variation), or a stretch of a few kilobases that differs
//! outright. The string graph may keep either, even one only a read or two have. Polishing cannot
//! mend that: reads of the other version align up to the locus's ends and stop. Here the reads'
//! mappings onto the polished contigs are searched for such places: pairs of positions at most
//! [`MAX_SEGMENT`] apart where many reads stop aligning. For each, the reads near it are mapped
//! onto the flanks either side, and the bases each read has between them are taken
//! ([`other_between`]). If their consensus is another version than the contig's, and most reads
//! have it, it replaces the contig's. The locus's exact edges need not be known: the reads' own
//! bases between the flanks hold them, and polishing follows. A read across the locus does not
//! align there whole, so the reads that map nowhere, or clipped, are mapped again, each up to the
//! place it stops aligning ([`with_parts`]).
//!
//! On a circular contig the flanks run on round its start, and on a small circle they are
//! shortened to fit ([`MIN_FLANK`]): a 13.3 kb *Pediococcus siamensis* LMG 24279 plasmid held a
//! 2.2 kb IS copy 9 reads have and 47 lack, 2 kb from its end, and 5 kb flanks did not fit.

use std::collections::HashSet;

use crate::consensus::PoaBackend;
use crate::overlap::{
    InvalidOverlapParamsError, Overlap, OverlapParams, ReadIndex, Reads, map_reads,
};

use crate::polish::Draft;

/// A read that maps nowhere, or with bases left over, is mapped again with its chains split
/// this many times at most ([`OverlapParams::max_chain_splits`]).
const MAX_CHAIN_SPLITS: u32 = 2;

/// Returns `mapped` (the mappings of `reads` onto `targets`) with each read that maps nowhere,
/// or with at least [`MIN_CLIP`] bases left over at an end, mapped again with its chains split
/// ([`MAX_CHAIN_SPLITS`]): a read across a locus where it has another version than the contig
/// then maps up to it, not to a copy of the rest elsewhere. Only for finding loci: polishing
/// keeps to whole chains, as the pieces of reads across many small differences mislead it.
///
/// The reads are mapped too onto each circular target (`circular`) turned half round, and a
/// read keeps that mapping if it aligns more of its bases there: a read across one's start then
/// maps whole, not to one side only, and its mapping can run past the target's length, which
/// [`repair`] folds back. On a 13.3 kb *P. siamensis* plasmid most reads cross the start, and
/// an IS copy 2 kb from its end showed 4 reads stopping at its far edge, against 18 at its near
/// one. (Mapping onto each circle twice round, end to end, did that too, but held every seed
/// twice: a *S. flexneri* chromosome lost a locus 24 of 24 reads across showed, for 2,342 errors.)
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn with_parts(
    targets: &Reads,
    circular: &[bool],
    reads: &Reads,
    mapped: &[Option<Overlap>],
    map: &OverlapParams,
) -> Result<Vec<Option<Overlap>>, InvalidOverlapParamsError> {
    let unmapped: Vec<ReadIndex> = reads
        .indices()
        .zip(mapped)
        .filter_map(|(r, m)| {
            let clipped = |m: &Overlap| {
                m.query_start >= MIN_CLIP || reads.read_len(r) - m.query_end >= MIN_CLIP
            };
            m.as_ref().is_none_or(clipped).then_some(r)
        })
        .collect();
    let again: Reads = unmapped.iter().map(|&r| reads.get(r)).collect();
    let split = OverlapParams {
        max_chain_splits: MAX_CHAIN_SPLITS,
        ..*map
    };
    // Each circular target turned half round: its base `i` is the target's `(i + shift) % len`.
    let shifts: Vec<u32> = targets
        .indices()
        .zip(circular)
        .map(|(t, &c)| if c { targets.read_len(t) / 2 } else { 0 })
        .collect();
    let turned: Vec<Vec<u8>> = targets
        .indices()
        .zip(&shifts)
        .map(|(t, &shift)| {
            let (head, tail) = targets.get(t).split_at(shift as usize);
            [tail, head].concat()
        })
        .collect();
    let turned: Reads = turned.iter().map(Vec::as_slice).collect();
    let as_is = map_reads(targets, &again, &split)?;
    let round = if shifts.iter().any(|&s| s > 0) {
        map_reads(&turned, &again, &split)?
    } else {
        vec![None; again.len()]
    };
    let aligned = |m: &Overlap| m.query_end - m.query_start;
    let mut out = mapped.to_vec();
    for ((&r, m), t) in unmapped.iter().zip(as_is).zip(round) {
        let t = t.map(|t| {
            let len = targets.read_len(t.target);
            let start = (t.target_start + shifts[t.target.index()]) % len;
            Overlap {
                target_start: start,
                target_end: start + (t.target_end - t.target_start),
                ..t
            }
        });
        let best = match (m, t) {
            (Some(m), Some(t)) if aligned(&t) > aligned(&m) => Some(t),
            (Some(m), _) => Some(m),
            (None, t) => t,
        };
        if let Some(m) = best {
            out[r.index()] = Some(Overlap { query: r, ..m });
        }
    }
    Ok(out)
}

/// A locus is at most this long.
pub(crate) const MAX_SEGMENT: u32 = 20_000;

/// A read stops aligning at a place when this many of its bases on that side are left out.
const MIN_CLIP: u32 = 300;

/// Stops are counted in bins of this many bases...
const BIN: u32 = 200;

/// ...and a bin is a candidate end when at least this many reads stop in it and its two
/// neighbours, more than in either of them...
const MIN_STOPS: usize = 5;

/// ...and at least this share (in percent) of the reads over it. A locus most reads have
/// another version of makes most of them stop at its ends; a fixed count alone is reached by
/// chance everywhere in deep data (at 92-116x, reads 7% off or of 650 bases on *Aerococcus*,
/// almost every bin), and each candidate is weighed by mapping the reads near it.
const MIN_STOP_PERCENT: usize = 20;

/// A read is weighed on a locus only if it is at least this many bases longer than the locus:
/// it must reach into the flanks either side to tell the versions apart.
const MIN_READ_PAST: u32 = 1_000;

/// A locus is weighed on at most this many reads near it, the longest: a majority of 100 is as
/// clear as one of thousands, and mapping thousands of reads per candidate took 79% of an
/// assembly at 288x.
const MAX_LOCUS_READS: usize = 100;

/// The window round a candidate holds this many bases either side of it...
const FLANK: u32 = 5_000;

/// ...or, on a circular contig too short for that, as many as fit (running on round its start),
/// but at least this many.
const MIN_FLANK: u32 = 2_000;

/// A read reaches a flank's inner end when its alignment stops at most this far short of it.
const FLANK_SLACK: u32 = 50;

/// The other version is kept only if this many more reads have it than the contig's (and at
/// least half as many again).
const MIN_GAIN: usize = 3;

/// The reads' bases between the flanks are taken from this far outside the places they stop.
const MARGIN: u32 = 600;

/// K-mers this long tell one version from the other.
const K: usize = 15;

/// The reads' version is the consensus of at most this many reads' bases.
const MAX_CONSENSUS_READS: usize = 20;

/// The versions differ in at least this many k-mers each (about 85 bases of other sequence),
/// or they are taken as the same: the consensus of a few reads has a few errors of its own.
/// One version may instead be the other without a stretch of at least [`MIN_CLIP`] bases (an
/// IS copy some cells carry): then only the longer needs this many, the shorter's own being
/// the few across its join.
const MIN_UNIQUE_KMERS: usize = 100;

fn kmers(seq: &[u8]) -> HashSet<&[u8]> {
    seq.windows(K).collect()
}

/// Of `reads`, those with the contig's version (`.0`) and those with the other version
/// (`.1`), given the k-mers only the contig's version has and those only the other has.
fn tally(
    reads: &[Vec<u8>],
    only_contig: &HashSet<&[u8]>,
    only_other: &HashSet<&[u8]>,
) -> (usize, usize) {
    let (mut contig, mut other) = (0, 0);
    for read in reads {
        let found = |set: &HashSet<&[u8]>| read.windows(K).filter(|k| set.contains(k)).count();
        // Shares of each version's own k-mers, compared without dividing: a * |B| vs b * |A|.
        let (a, b) = (
            found(only_contig) * only_other.len(),
            found(only_other) * only_contig.len(),
        );
        if a > b {
            contig += 1;
        } else if b > a {
            other += 1;
        }
    }
    (contig, other)
}

/// Returns what most reads have between `x` and `y` of `seq` if that is another version than
/// the contig's; `None` if most have the contig's.
///
/// The reads near are mapped onto the flanks either side (`flank` bases), and of each read
/// that maps onto both on one strand, up to the flanks' inner ends, the bases between are
/// taken. Their consensus is the version most of them have. If it differs from the contig's
/// bases there (in at least [`MIN_UNIQUE_KMERS`] k-mers each way, or the longer one does and
/// the other is at least [`MIN_CLIP`] bases shorter), every read is counted for
/// the version whose own k-mers it has more of, and the reads' version is returned if at least
/// [`MIN_GAIN`] more reads, and half as many again, have it.
fn other_between(
    seq: &[u8],
    (x, y): (u32, u32),
    flank: u32,
    reads: &Reads,
    map: &OverlapParams,
    backend: &PoaBackend,
) -> Option<Vec<u8>> {
    let between = bases_between(seq, (x, y), flank, reads, map)?;
    if between.len() < MIN_GAIN {
        return None;
    }
    let here = &seq[x as usize..y as usize];
    let consensus = backend.consensus_segments(&between[..between.len().min(MAX_CONSENSUS_READS)]);
    let (contig_kmers, other_kmers) = (kmers(here), kmers(&consensus));
    let only_contig: HashSet<&[u8]> = contig_kmers.difference(&other_kmers).copied().collect();
    let only_other: HashSet<&[u8]> = other_kmers.difference(&contig_kmers).copied().collect();
    let (longer, shorter) = if here.len() >= consensus.len() {
        (&only_contig, &only_other)
    } else {
        (&only_other, &only_contig)
    };
    let without_stretch = here.len().abs_diff(consensus.len()) >= MIN_CLIP as usize;
    let differ = if without_stretch {
        longer.len() >= MIN_UNIQUE_KMERS && !shorter.is_empty()
    } else {
        only_contig.len() >= MIN_UNIQUE_KMERS && only_other.len() >= MIN_UNIQUE_KMERS
    };
    if !differ {
        return None;
    }
    let (same, other) = tally(&between, &only_contig, &only_other);
    log::debug!("bacteria: {x}..{y}: {same} reads have the contig's version, {other} the other");
    (other >= same + MIN_GAIN && 2 * other >= 3 * same).then_some(consensus)
}

/// The bases each of `reads` has between `x` and `y` of `seq`: the reads are mapped onto the
/// flanks either side (`flank` bases), and of each that maps onto both on one strand, up to
/// the flanks' inner ends, the bases between are taken (`None` if `map` cannot be used).
fn bases_between(
    seq: &[u8],
    (x, y): (u32, u32),
    flank: u32,
    reads: &Reads,
    map: &OverlapParams,
) -> Option<Vec<Vec<u8>>> {
    let (xs, ys, f) = (x as usize, y as usize, flank as usize);
    let left: Reads = [&seq[xs - f..xs]].into_iter().collect();
    let right: Reads = [&seq[ys..ys + f]].into_iter().collect();
    let on_left = map_reads(&left, reads, map).ok()?;
    let on_right = map_reads(&right, reads, map).ok()?;

    let mut between = Vec::new();
    for ((read, l), r) in reads.indices().zip(&on_left).zip(&on_right) {
        let (Some(l), Some(r)) = (l, r) else { continue };
        if l.strand != r.strand
            || l.target_end + FLANK_SLACK < flank
            || r.target_start > FLANK_SLACK
        {
            continue;
        }
        let bases = reads.get(read);
        let len = reads.read_len(read);
        let (oriented, l_end, r_start) = if l.strand.is_reverse() {
            (
                crate::seq::reverse_complement(bases),
                len - l.query_start,
                len - r.query_end,
            )
        } else {
            (bases.to_vec(), l.query_end, r.query_start)
        };
        let from = (l_end + (flank - l.target_end)) as usize;
        let to = (r_start as usize).saturating_sub(r.target_start as usize);
        if from < to {
            between.push(oriented[from..to].to_vec());
        }
    }
    log::debug!(
        "bacteria: {x}..{y}: {} reads near, {} across both flanks",
        reads.len(),
        between.len()
    );
    Some(between)
}

/// `seq`'s bases from `flank` before `x` to `flank` past `y`, running on round its start (a
/// circular contig at least `y - x + 2 * flank` long).
fn round_window(seq: &[u8], (x, y): (u32, u32), flank: u32) -> Vec<u8> {
    let (len, x, y, f) = (seq.len(), x as usize, y as usize, flank as usize);
    (x + len - f..y + len + f).map(|i| seq[i % len]).collect()
}

/// Replaces loci of `drafts` with the version most of `reads` have (see the module docs), given
/// the reads' `mappings` onto them; returns the loci replaced.
pub(crate) fn repair(
    drafts: &mut [Draft],
    reads: &Reads,
    mappings: &[Option<Overlap>],
    map: &OverlapParams,
    backend: &PoaBackend,
) -> usize {
    let mut repaired = 0;
    for (contig, draft) in drafts.iter_mut().enumerate() {
        let len = u32::try_from(draft.seq.len()).expect("contigs are shorter than 4 Gb");
        let min_flank = if draft.circular { MIN_FLANK } else { FLANK };
        if len < 2 * (min_flank + MARGIN + BIN) {
            continue;
        }
        // Where reads mapped onto this contig stop aligning with bases of theirs left over, and
        // how many reads lie over each bin.
        let mut stops = vec![0_usize; (len / BIN + 1) as usize];
        let mut over = vec![0_i64; (len / BIN + 2) as usize];
        let mut near: Vec<(u32, u32, ReadIndex)> = Vec::new();
        for (read, m) in reads.indices().zip(mappings) {
            let Some(m) = m else { continue };
            if m.target.index() != contig {
                continue;
            }
            let read_len = reads.read_len(read);
            let (left, right) = if m.strand.is_reverse() {
                (read_len - m.query_end, m.query_start)
            } else {
                (m.query_start, read_len - m.query_end)
            };
            // A circular contig's mapping may run on past its length ([`with_parts`]): its
            // start folded back onto the contig, its end at most a length further on.
            let start = m.target_start % len;
            let end = start + (m.target_end - m.target_start);
            let folded_end = if end > len { end - len } else { end };
            if left >= MIN_CLIP && start >= MIN_CLIP {
                stops[(start / BIN) as usize] += 1;
            }
            if right >= MIN_CLIP && len - folded_end >= MIN_CLIP {
                stops[(folded_end / BIN) as usize] += 1;
            }
            over[(start / BIN) as usize] += 1;
            if end > len {
                over[(len / BIN) as usize + 1] -= 1;
                over[0] += 1;
            }
            over[(folded_end / BIN) as usize + 1] -= 1;
            near.push((start, end, read));
        }
        let depth: Vec<usize> = over
            .iter()
            .scan(0_i64, |sum, &d| {
                *sum += d;
                Some(usize::try_from(*sum).unwrap_or(0))
            })
            .collect();
        // Stops at a locus's end spread over neighbouring bins (inside an inverted repeat, or
        // where the versions start to differ): each bin counts its neighbours' too, and the
        // local peaks are the candidates.
        let smooth: Vec<usize> = (0..stops.len())
            .map(|b| {
                stops[b.saturating_sub(1)..(b + 2).min(stops.len())]
                    .iter()
                    .sum()
            })
            .collect();
        let ends: Vec<u32> = (0..smooth.len())
            .filter(|&b| {
                smooth[b] >= MIN_STOPS
                    && smooth[b] * 100 >= MIN_STOP_PERCENT * depth[b]
                    && (b == 0 || smooth[b] >= smooth[b - 1])
                    && (b + 1 == smooth.len() || smooth[b] > smooth[b + 1])
            })
            .map(|b| u32::try_from(b).expect("fits") * BIN + BIN / 2)
            .collect();
        log::debug!(
            "bacteria: contig {contig}: {} places where reads stop: {:?}",
            ends.len(),
            ends.iter().take(40).collect::<Vec<_>>()
        );
        // Candidate loci: two such places, the nearest first, whose window (and its left
        // flank) starts past the last one replaced. Replacements are made last to first, so
        // earlier places do not move.
        let mut replacements: Vec<(usize, usize, Vec<u8>)> = Vec::new();
        let mut done_until = 0;
        for (a, &cx) in ends.iter().enumerate() {
            if cx < done_until + MARGIN {
                continue;
            }
            for &cy in ends[a + 1..]
                .iter()
                .filter(|&&y| y - cx >= 2 * BIN && y - cx <= MAX_SEGMENT)
            {
                // Edges safely outside the locus: the reads' own bases between the flanks
                // settle where it starts and ends.
                let (Some(x), y) = (cx.checked_sub(MARGIN), cy + MARGIN) else {
                    continue;
                };
                // The flanks: inside the contig, or round a circle's start.
                let inside = x >= FLANK && y + FLANK <= len;
                let flank = if inside {
                    FLANK
                } else if draft.circular {
                    FLANK.min((len - (y - x)) / 2)
                } else {
                    0
                };
                if flank < MIN_FLANK {
                    continue;
                }
                let window;
                let (seq, at) = if inside {
                    (draft.seq.as_slice(), (x, y))
                } else {
                    window = round_window(&draft.seq, (x, y), flank);
                    (window.as_slice(), (flank, flank + y - x))
                };
                // The reads near long enough to reach past the locus into both flanks, the
                // longest first.
                let (lo, hi, l) = (
                    i64::from(x) - i64::from(flank),
                    i64::from(y) + i64::from(flank),
                    i64::from(len),
                );
                let mut long: Vec<ReadIndex> = near
                    .iter()
                    .filter(|&&(s, e, r)| {
                        let (s, e) = (i64::from(s), i64::from(e));
                        [-l, 0, l]
                            .iter()
                            .any(|&turn| s < hi + turn && e > lo + turn)
                            && reads.read_len(r) >= y - x + MIN_READ_PAST
                    })
                    .map(|&(_, _, r)| r)
                    .collect();
                if long.len() < MIN_GAIN {
                    continue;
                }
                long.sort_by_key(|&r| (std::cmp::Reverse(reads.read_len(r)), r));
                long.truncate(MAX_LOCUS_READS);
                let reads_here: Reads = long.iter().map(|&r| reads.get(r)).collect();
                if let Some(seq) = other_between(seq, at, flank, &reads_here, map, backend) {
                    log::info!(
                        "bacteria: contig {contig}: the {} bases at {x} replaced by the {} most \
                         reads have",
                        y - x,
                        seq.len()
                    );
                    replacements.push((x as usize, y as usize, seq));
                    repaired += 1;
                    done_until = y + FLANK;
                    break;
                }
            }
        }
        for (x, y, seq) in replacements.into_iter().rev() {
            draft.seq.splice(x..y, seq);
        }
    }
    repaired
}

/// A tandem copy is at least this long...
const MIN_UNIT: usize = 500;

/// ...and at most this long.
const MAX_UNIT: usize = 10_000;

/// Tandem copies are found from k-mers this long that a contig holds exactly twice.
const UNIT_K: usize = 21;

/// Two k-mers' pairs belong to one pair of copies if their spacings differ by at most this
/// much (the copies' indels)...
const UNIT_SLACK: usize = 20;

/// ...and they lie at most this far apart.
const UNIT_GAP: usize = 200;

/// Places in `seq` where a stretch is followed by a copy of itself: `(start, unit)`, the first
/// copy at `start..start + unit` and the second right after, from runs of k-mers held exactly
/// twice, `unit` apart, along most of a unit.
fn tandem_copies(seq: &[u8]) -> Vec<(usize, usize)> {
    let mut seen: std::collections::HashMap<&[u8], (usize, usize, u8)> =
        std::collections::HashMap::new();
    for (i, k) in seq.windows(UNIT_K).enumerate() {
        let e = seen.entry(k).or_insert((i, i, 0));
        e.1 = i;
        e.2 = e.2.saturating_add(1);
    }
    let mut pairs: Vec<(usize, usize)> = seen
        .values()
        .filter(|&&(a, b, n)| n == 2 && (MIN_UNIT..=MAX_UNIT).contains(&(b - a)))
        .map(|&(a, b, _)| (b - a, a))
        .collect();
    pairs.sort_unstable();
    let mut found = Vec::new();
    let mut i = 0;
    while i < pairs.len() {
        // A run of pairs at about one spacing, each near the last.
        let (unit, start) = pairs[i];
        let (mut end, mut j) = (start, i + 1);
        while j < pairs.len()
            && pairs[j].0 <= unit + UNIT_SLACK
            && pairs[j].1 <= end + UNIT_GAP
            && pairs[j].1 >= start
        {
            end = end.max(pairs[j].1);
            j += 1;
        }
        // Along most of a unit, and the two copies side by side.
        if end - start >= unit * 4 / 5 && start + 2 * unit <= seq.len() {
            found.push((start, unit));
        }
        i = j;
    }
    found.sort_unstable();
    found.dedup_by(|b, a| b.0 < a.0 + a.1);
    found
}

/// Leaves one copy where a contig holds a stretch twice in a row ([`tandem_copies`]) but most
/// reads across it have it once (part of the cells carry a duplication: a 3.2 kb one in a
/// *Borrelia* chromosome, which a quarter of the reads have). The reads near are taken between
/// flanks outside the two copies ([`bases_between`]); each counts for the version whose length
/// it is closer to, and one copy goes if at least [`MIN_GAIN`] more reads, and half as many
/// again, have one. The k-mer vote of [`other_between`] cannot tell them apart: the two
/// versions share all but a few k-mers. Returns the copies left out.
pub(crate) fn drop_minor_copies(
    drafts: &mut [Draft],
    reads: &Reads,
    mappings: &[Option<Overlap>],
    map: &OverlapParams,
) -> usize {
    let mut dropped = 0;
    for (contig, draft) in drafts.iter_mut().enumerate() {
        let len = draft.seq.len();
        let mut cuts: Vec<(usize, usize)> = Vec::new();
        for (start, unit) in tandem_copies(&draft.seq) {
            let (margin, flank) = (MARGIN as usize, FLANK as usize);
            let (Some(x), y) = (start.checked_sub(margin), start + 2 * unit + margin) else {
                continue;
            };
            if x < flank || y + flank > len {
                continue;
            }
            let near: Reads = reads
                .indices()
                .zip(mappings)
                .filter_map(|(r, m)| {
                    let m = m.as_ref()?;
                    let (s, e) = (m.target_start as usize, m.target_end as usize);
                    (m.target.index() == contig && s < y + flank && e + flank > x)
                        .then(|| reads.get(r))
                })
                .collect();
            let (xu, yu) = (
                u32::try_from(x).expect("contigs are shorter than 4 Gb"),
                u32::try_from(y).expect("contigs are shorter than 4 Gb"),
            );
            let Some(between) = bases_between(&draft.seq, (xu, yu), FLANK, &near, map) else {
                continue;
            };
            let (once, twice) = copies_in(&between, y - x - unit, unit);
            log::debug!(
                "bacteria: contig {contig}: a {unit}-base stretch twice at {start}: {twice} reads \
                 have it twice, {once} once"
            );
            if once >= twice + MIN_GAIN && 2 * once >= 3 * twice {
                log::info!(
                    "bacteria: contig {contig}: one copy of the {unit}-base stretch at {start} left \
                     out ({once} reads have it once, {twice} twice)"
                );
                cuts.push((start, unit));
            }
        }
        for &(start, unit) in cuts.iter().rev() {
            draft.seq.drain(start..start + unit);
            dropped += 1;
        }
    }
    dropped
}

/// Of the reads' bases `between` two flanks, those that hold a stretch of `unit` bases once
/// (`.0`: about `one` bases long) and twice (`.1`): each counts for the length it is closer to,
/// if within a quarter of a unit.
fn copies_in(between: &[Vec<u8>], one: usize, unit: usize) -> (usize, usize) {
    let two = one + unit;
    let (mut once, mut twice) = (0, 0);
    for b in between {
        let l = b.len();
        if l.abs_diff(one) < l.abs_diff(two) && l.abs_diff(one) <= unit / 4 {
            once += 1;
        } else if l.abs_diff(two) <= l.abs_diff(one) && l.abs_diff(two) <= unit / 4 {
            twice += 1;
        }
    }
    (once, twice)
}

/// A circle that is a copy of a stretch a longer draft holds once, to be put back beside it.
struct Lost {
    circle: usize,
    host: usize,
    /// Where the host's copy of the circle's bases starts, and they themselves from there on
    /// round the circle, in the host's orientation.
    at: usize,
    copy: Vec<u8>,
}

/// Puts back a tandem copy laid out as a circle of its own: a circular draft of
/// [`MIN_UNIT`]..=[`MAX_UNIT`] bases that a longer draft holds once, where most of the cells
/// have the stretch twice in a row (two 6 kb rRNA operons side by side in an *S. pyogenes*
/// chromosome: the string graph walks the loop once, and the second copy's reads lay out a
/// circle, which was left out at the end as an element excised from the chromosome).
///
/// Reads across both copies and their flanks are rare (12 kb and more), so the reads near are
/// mapped onto the host's stretch twice in a row, flanks either side, and those across the join
/// of the two copies are weighed against those across the join of a copy and a flank (each
/// past it by [`MARGIN`]): a cell with the stretch once or twice has the latter, only one with
/// it twice the former. The copy goes back if at least [`MIN_GAIN`] reads are across the copies'
/// join, and at least half as many as across a flank's. An element that excises leaves the
/// flanks joined to each other in the cells it left: the copy stays out if more than a tenth as
/// many reads have the host's stretch not at all ([`bases_between`]) as are across a flank's
/// join. Returns which drafts are kept and how many copies were put back.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn restore_copies(
    drafts: &mut [Draft],
    reads: &Reads,
    mappings: &[Option<Overlap>],
    map: &OverlapParams,
) -> Result<(Vec<bool>, usize), InvalidOverlapParamsError> {
    let mut kept = vec![true; drafts.len()];
    let (margin, flank) = (MARGIN as usize, FLANK as usize);
    let mut lost: Vec<Lost> = Vec::new();
    for circle in 0..drafts.len() {
        let unit = drafts[circle].seq.len();
        if !drafts[circle].circular || !(MIN_UNIT..=MAX_UNIT).contains(&unit) {
            continue;
        }
        let hosts: Vec<usize> = (0..drafts.len())
            .filter(|&h| drafts[h].seq.len() > unit + 2 * (flank + margin))
            .collect();
        if hosts.is_empty() {
            continue;
        }
        // Twice round, so the circle maps whole wherever the host's copy starts.
        let targets: Reads = hosts.iter().map(|&h| drafts[h].seq.as_slice()).collect();
        let round = [drafts[circle].seq.as_slice(), &drafts[circle].seq].concat();
        let query: Reads = [round.as_slice()].into_iter().collect();
        let Some(m) = map_reads(&targets, &query, map)?.pop().flatten() else {
            continue;
        };
        let host = hosts[m.target.index()];
        let (from, to) = (pos(m.target_start), pos(m.target_end));
        let seq = &drafts[host].seq;
        if (to - from).abs_diff(unit) > unit / 20
            || from < flank + margin
            || to + margin + flank > seq.len()
        {
            continue;
        }
        // The circle's bases from the one at `from` on, once round: a copy that, put in at
        // `from`, makes the stretch twice in a row exactly, whatever the mapping's ends miss.
        let (oriented, start) = if m.strand.is_reverse() {
            (
                crate::seq::reverse_complement(&round),
                2 * unit - pos(m.query_end),
            )
        } else {
            (round, pos(m.query_start))
        };
        let copy = oriented[start % unit..start % unit + unit].to_vec();
        let (x, y) = (from - margin, to + margin);
        // With the reads that map nowhere: a read without the stretch, beside a host that has
        // it, may differ from it too much along the gap to map.
        let near: Reads = reads
            .indices()
            .zip(mappings)
            .filter_map(|(r, m)| {
                let Some(m) = m else {
                    return Some(reads.get(r));
                };
                let (s, e, d) = (pos(m.target_start), pos(m.target_end), m.target.index());
                (d == circle || (d == host && s < y + flank && e + flank > x)).then(|| reads.get(r))
            })
            .collect();
        // The stretch twice, flanks either side: the copies' join at `join`, the flanks' about
        // a unit either side.
        let window = [&seq[x - flank..from], &copy, &seq[from..y + flank]].concat();
        let join = from + unit - (x - flank);
        let on_window: Reads = [window.as_slice()].into_iter().collect();
        let across = |at: usize, placed: &[Option<Overlap>]| {
            placed
                .iter()
                .flatten()
                .filter(|m| pos(m.target_start) + margin <= at && pos(m.target_end) >= at + margin)
                .count()
        };
        let placed = map_reads(&on_window, &near, map)?;
        let copies = across(join, &placed);
        let flanks = across(join - unit, &placed).max(across(join + (to - from), &placed));
        let (xu, yu) = (
            u32::try_from(x).expect("contigs are shorter than 4 Gb"),
            u32::try_from(y).expect("contigs are shorter than 4 Gb"),
        );
        let Some(between) = bases_between(seq, (xu, yu), FLANK, &near, map) else {
            continue;
        };
        let none = between
            .iter()
            .filter(|b| b.len().abs_diff(y - x - unit) <= unit / 4)
            .count();
        log::debug!(
            "bacteria: circle {circle} ({unit} bases) is {from}..{to} of draft {host}: {copies} \
             reads across the copies' join, {flanks} across a flank's, {none} without it"
        );
        if copies >= MIN_GAIN && 2 * copies >= flanks && 10 * none <= flanks {
            lost.push(Lost {
                circle,
                host,
                at: from,
                copy,
            });
        }
    }
    // From the end of each host back, so the places of the stretches still to copy hold.
    lost.sort_unstable_by_key(|c| (c.host, std::cmp::Reverse(c.at)));
    let mut restored = 0;
    for c in lost {
        if !kept[c.circle] || !kept[c.host] {
            continue;
        }
        let unit = c.copy.len();
        drafts[c.host].seq.splice(c.at..c.at, c.copy);
        kept[c.circle] = false;
        restored += 1;
        log::info!(
            "bacteria: draft {}: the {unit}-base stretch at {} put back twice; circle {} left out",
            c.host,
            c.at,
            c.circle
        );
    }
    Ok((kept, restored))
}

/// The stretch between the two places a version's ends lie in a longer draft is at most this
/// long.
const MAX_VERSION_SPAN: usize = 50_000;

/// A linear draft that is a second version of a stretch of a longer draft.
struct Version {
    /// The linear draft.
    line: usize,
    /// The longer draft, and where in it the stretch lies.
    host: usize,
    stretch: (usize, usize),
    /// The linear draft's own bases in its place, in the longer draft's orientation.
    bases: Vec<u8>,
    /// Reads across the linear draft's version, and across the longer draft's.
    with_line: usize,
    with_host: usize,
}

/// Settles loci laid out in two versions: a linear draft whose first and last [`FLANK`] bases
/// lie whole in a longer draft, on one strand, in order and at most [`MAX_VERSION_SPAN`] apart,
/// with other bases between (part of the cells carry another version of an 11.5 kb stretch of a
/// *B. arboris* chromosome). The reads of each version map onto the draft that has it whole, so
/// neither shows reads stopping and [`repair`] cannot see the locus. The reads across each
/// version (past its flanks by [`MARGIN`]) are counted; if at least [`MIN_GAIN`] more, and half
/// as many again, are across the linear draft's, its bases replace the longer draft's, and if
/// either version wins so, the linear draft is left out. Returns which drafts are kept and how
/// many versions were settled.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn settle_versions(
    drafts: &mut [Draft],
    mappings: &[Option<Overlap>],
    map: &OverlapParams,
) -> Result<(Vec<bool>, usize), InvalidOverlapParamsError> {
    let flank = FLANK as usize;
    let lines: Vec<usize> = (0..drafts.len())
        .filter(|&i| {
            let d = &drafts[i];
            !d.circular && d.hairpins == [false; 2] && d.seq.len() >= 3 * flank
        })
        .collect();
    let mut kept = vec![true; drafts.len()];
    if lines.is_empty() {
        return Ok((kept, 0));
    }
    // Each line's own ends masked, so its anchors map elsewhere or nowhere.
    let mut masked: Vec<Vec<u8>> = drafts.iter().map(|d| d.seq.clone()).collect();
    for &l in &lines {
        let seq = &mut masked[l];
        let len = seq.len();
        seq[..flank].fill(b'N');
        seq[len - flank..].fill(b'N');
    }
    let targets: Reads = masked.iter().map(Vec::as_slice).collect();
    drop(masked);
    let mut anchors = Reads::new();
    for &l in &lines {
        let seq = &drafts[l].seq;
        anchors.push(&seq[..flank]);
        anchors.push(&seq[seq.len() - flank..]);
    }
    let placed = map_reads(&targets, &anchors, map)?;
    let whole = |m: &Overlap| pos(m.query_end - m.query_start) + 2 * MARGIN as usize >= flank;
    // Reads across `from..to` of draft `d`, past it by a margin either side.
    let across = |d: usize, from: usize, to: usize| {
        mappings
            .iter()
            .flatten()
            .filter(|m| {
                m.target.index() == d
                    && pos(m.target_start) + MARGIN as usize <= from
                    && pos(m.target_end) >= to + MARGIN as usize
            })
            .count()
    };
    let mut versions: Vec<Version> = Vec::new();
    for (k, &l) in lines.iter().enumerate() {
        let (Some(first), Some(last)) = (&placed[2 * k], &placed[2 * k + 1]) else {
            continue;
        };
        let host = first.target.index();
        let len = drafts[l].seq.len();
        if host == l
            || last.target.index() != host
            || first.strand != last.strand
            || !whole(first)
            || !whole(last)
            || drafts[host].seq.len() <= len
        {
            continue;
        }
        // The stretch between the two anchors in the host, and the line's bases there.
        let middle = &drafts[l].seq[flank..len - flank];
        let (from, to, bases) = if first.strand.is_reverse() {
            (
                pos(last.target_end),
                pos(first.target_start),
                crate::seq::reverse_complement(middle),
            )
        } else {
            (
                pos(first.target_end),
                pos(last.target_start),
                middle.to_vec(),
            )
        };
        if from > to || to - from > MAX_VERSION_SPAN || drafts[host].seq[from..to] == bases[..] {
            continue;
        }
        // The reads across the part that differs: long reads across the whole stretch are few.
        let (same_start, same_end) = alike_ends(&drafts[host].seq[from..to], &bases);
        let (line_from, line_to) = if first.strand.is_reverse() {
            (flank + same_end, len - flank - same_start)
        } else {
            (flank + same_start, len - flank - same_end)
        };
        versions.push(Version {
            line: l,
            host,
            stretch: (from, to),
            with_line: across(l, line_from, line_to.max(line_from)),
            with_host: across(
                host,
                from + same_start,
                (to - same_end).max(from + same_start),
            ),
            bases,
        });
    }
    // From the end of each host back, so the places of the stretches still to replace hold.
    versions.sort_unstable_by_key(|v| (v.host, std::cmp::Reverse(v.stretch.0)));
    let mut settled = 0;
    let mut last_change: Option<(usize, usize)> = None;
    for v in versions {
        let wins = |a: usize, b: usize| a >= b + MIN_GAIN && 2 * a >= 3 * b;
        log::debug!(
            "bacteria: draft {} is another version of {}..{} of draft {}: {} reads across it, {} \
             across the draft's",
            v.line,
            v.stretch.0,
            v.stretch.1,
            v.host,
            v.with_line,
            v.with_host
        );
        if !kept[v.line] {
            continue;
        }
        if wins(v.with_line, v.with_host) {
            // One replacement per stretch of a host: another line's places would have moved.
            if last_change.is_some_and(|(h, from)| h == v.host && v.stretch.1 > from) {
                continue;
            }
            let (from, to) = v.stretch;
            drafts[v.host].seq.splice(from..to, v.bases);
            last_change = Some((v.host, from));
        } else if !wins(v.with_host, v.with_line) {
            continue;
        }
        kept[v.line] = false;
        settled += 1;
    }
    Ok((kept, settled))
}

/// Two versions' bases are alike up to where [`ALIKE_RUN`] bases in a row past a difference
/// match on no shift up to [`ALIKE_SHIFT`] (a polishing error is not where they differ).
const ALIKE_RUN: usize = 50;
const ALIKE_SHIFT: usize = 3;

/// How many bases at the start, and at the end, `a` and `b` have alike (small differences
/// aside), the start's and the end's not crossing.
fn alike_ends(a: &[u8], b: &[u8]) -> (usize, usize) {
    let alike = |a: &[u8], b: &[u8]| -> usize {
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            if a[i] == b[j] {
                (i, j) = (i + 1, j + 1);
                continue;
            }
            // A substitution, or an indel of up to ALIKE_SHIFT bases either way, then a run.
            let run = |x: usize, y: usize| {
                x + ALIKE_RUN <= a.len()
                    && y + ALIKE_RUN <= b.len()
                    && a[x..x + ALIKE_RUN] == b[y..y + ALIKE_RUN]
            };
            let next = (0..=ALIKE_SHIFT)
                .flat_map(|d| {
                    [
                        (i + 1, j + 1 + d),
                        (i + 1 + d, j + 1),
                        (i + d, j),
                        (i, j + d),
                    ]
                })
                .find(|&(x, y)| (x, y) != (i, j) && run(x, y));
            match next {
                Some((x, y)) => (i, j) = (x, y),
                None => break,
            }
        }
        i.min(j)
    };
    let start = alike(a, b);
    let (ra, rb): (Vec<u8>, Vec<u8>) = (
        a.iter().rev().copied().collect(),
        b.iter().rev().copied().collect(),
    );
    let end = alike(&ra, &rb).min(a.len().min(b.len()) - start);
    (start, end)
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use crate::consensus::PoaParams;

    use super::*;
    use crate::{BactParams, ReadKind};

    fn random(len: usize, seed: u64) -> Vec<u8> {
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

    /// A genome of 40 kb with, at 18 kb, `locus`.
    fn with_locus(locus: &[u8]) -> Vec<u8> {
        let flanks = random(40_000, 7);
        [&flanks[..18_000], locus, &flanks[18_000..]].concat()
    }

    /// Reads of 9 kb along `genome`, starting every `step` bases, every other one turned round.
    fn reads_of(genome: &[u8], step: usize) -> Vec<Vec<u8>> {
        (0..genome.len() - 9_000)
            .step_by(step)
            .enumerate()
            .map(|(i, at)| {
                let read = &genome[at..at + 9_000];
                if i % 2 == 0 {
                    read.to_vec()
                } else {
                    crate::seq::reverse_complement(read)
                }
            })
            .collect()
    }

    #[test]
    fn a_locus_most_reads_have_otherwise_gets_their_version() {
        let (rare, common) = (random(1_800, 1), random(1_700, 2));
        let mut all = reads_of(&with_locus(&common), 300);
        all.extend(reads_of(&with_locus(&rare), 3_000));
        let reads: Reads = all.iter().map(Vec::as_slice).collect();
        let params = BactParams::new(ReadKind::Ont);
        let backend = PoaBackend::new(PoaParams::default()).expect("default settings");
        let mut drafts = [Draft::new(with_locus(&rare), false)];
        let targets: Reads = [drafts[0].seq.as_slice()].into_iter().collect();
        let mapped = map_reads(&targets, &reads, &params.map).expect("valid settings");
        let mappings =
            with_parts(&targets, &[false], &reads, &mapped, &params.map).expect("valid settings");
        assert_eq!(
            repair(&mut drafts, &reads, &mappings, &params.map, &backend),
            1
        );
        assert_eq!(drafts[0].seq, with_locus(&common));
    }

    #[test]
    fn a_locus_most_reads_have_as_the_contig_is_left_alone() {
        let (rare, common) = (random(1_800, 1), random(1_700, 2));
        let mut all = reads_of(&with_locus(&common), 300);
        all.extend(reads_of(&with_locus(&rare), 3_000));
        let reads: Reads = all.iter().map(Vec::as_slice).collect();
        let params = BactParams::new(ReadKind::Ont);
        let backend = PoaBackend::new(PoaParams::default()).expect("default settings");
        let mut drafts = [Draft::new(with_locus(&common), false)];
        let targets: Reads = [drafts[0].seq.as_slice()].into_iter().collect();
        let mapped = map_reads(&targets, &reads, &params.map).expect("valid settings");
        let mappings =
            with_parts(&targets, &[false], &reads, &mapped, &params.map).expect("valid settings");
        assert_eq!(
            repair(&mut drafts, &reads, &mappings, &params.map, &backend),
            0
        );
        assert_eq!(drafts[0].seq, with_locus(&common));
    }

    /// Reads of 9 kb round the circle `genome`, starting every `step` bases, every other one
    /// turned round.
    fn reads_round(genome: &[u8], step: usize) -> Vec<Vec<u8>> {
        let twice = [genome, genome].concat();
        (0..genome.len())
            .step_by(step)
            .enumerate()
            .map(|(i, at)| {
                let read = &twice[at..at + 9_000];
                if i % 2 == 0 {
                    read.to_vec()
                } else {
                    crate::seq::reverse_complement(read)
                }
            })
            .collect()
    }

    #[test]
    fn a_stretch_most_reads_lack_near_a_small_circles_start_is_cut() {
        // A 13 kb plasmid with a 2.2 kb insertion 2 kb before its end, which few reads have.
        let plasmid = random(11_000, 3);
        let without = plasmid.clone();
        let with = [&plasmid[..9_000], &random(2_200, 4), &plasmid[9_000..]].concat();
        let mut all = reads_round(&without, 250);
        all.extend(reads_round(&with, 2_500));
        let reads: Reads = all.iter().map(Vec::as_slice).collect();
        let params = BactParams::new(ReadKind::Ont);
        let backend = PoaBackend::new(PoaParams::default()).expect("default settings");
        let mut drafts = [Draft::new(with.clone(), true)];
        let targets: Reads = [drafts[0].seq.as_slice()].into_iter().collect();
        let mapped = map_reads(&targets, &reads, &params.map).expect("valid settings");
        let mappings =
            with_parts(&targets, &[true], &reads, &mapped, &params.map).expect("valid settings");
        assert_eq!(
            repair(&mut drafts, &reads, &mappings, &params.map, &backend),
            1
        );
        assert_eq!(drafts[0].seq, without);

        // The same circle, its insertion kept when most reads have it.
        let mut all = reads_round(&with, 250);
        all.extend(reads_round(&without, 2_500));
        let reads: Reads = all.iter().map(Vec::as_slice).collect();
        let mut drafts = [Draft::new(with.clone(), true)];
        let mapped = map_reads(&targets, &reads, &params.map).expect("valid settings");
        let mappings =
            with_parts(&targets, &[true], &reads, &mapped, &params.map).expect("valid settings");
        assert_eq!(
            repair(&mut drafts, &reads, &mappings, &params.map, &backend),
            0
        );
        assert_eq!(drafts[0].seq, with);
    }

    #[test]
    fn a_window_round_a_circles_start_runs_on_past_its_end() {
        assert_eq!(round_window(b"ABCDEFGHIJ", (1, 3), 2), b"JABCDE");
        assert_eq!(round_window(b"ABCDEFGHIJ", (7, 9), 2), b"FGHIJA");
    }

    #[test]
    fn reads_are_counted_for_the_version_whose_own_kmers_they_have() {
        let (contig, other) = (random(2_000, 1), random(2_000, 2));
        let (c, o) = (kmers(&contig), kmers(&other));
        let only_contig: HashSet<&[u8]> = c.difference(&o).copied().collect();
        let only_other: HashSet<&[u8]> = o.difference(&c).copied().collect();
        let reads = vec![
            other.clone(),
            other[100..1_500].to_vec(),
            contig[..900].to_vec(),
            b"ACGT".to_vec(),
        ];
        assert_eq!(tally(&reads, &only_contig, &only_other), (1, 2));
    }

    #[test]
    fn a_tandem_copy_most_reads_lack_is_left_out() {
        // A stretch of 2 kb twice in a row in the contig; in the reads once (every 300 bases)
        // or twice (every 1,500), and the other way round.
        let unit = random(2_000, 3);
        let once = with_locus(&unit);
        let twice = with_locus(&[unit.as_slice(), &unit].concat());
        let params = BactParams::new(ReadKind::Ont);
        for (common, rare, dropped) in [(&once, &twice, 1), (&twice, &once, 0)] {
            let mut all = reads_of(common, 300);
            all.extend(reads_of(rare, 1_500));
            let reads: Reads = all.iter().map(Vec::as_slice).collect();
            let mut drafts = [Draft::new(twice.clone(), false)];
            let targets: Reads = [drafts[0].seq.as_slice()].into_iter().collect();
            let mapped = map_reads(&targets, &reads, &params.map).expect("valid settings");
            assert_eq!(tandem_copies(&twice), vec![(18_000, 2_000)]);
            assert_eq!(
                drop_minor_copies(&mut drafts, &reads, &mapped, &params.map),
                dropped
            );
            assert_eq!(
                drafts[0].seq,
                if dropped == 1 {
                    once.clone()
                } else {
                    twice.clone()
                }
            );
        }
    }

    #[test]
    fn a_tandem_copy_laid_out_as_a_circle_is_put_back_but_an_excised_one_is_not() {
        // A stretch of 2 kb twice in a row in the genome, laid out once in the contig and as a
        // circle of its own. Then the same stretch once, excised in part of the cells: reads
        // without it, and round its circle.
        let unit = random(2_000, 3);
        let once = with_locus(&unit);
        let twice = with_locus(&[unit.as_slice(), &unit].concat());
        let round = unit.repeat(6);
        let params = BactParams::new(ReadKind::Ont);
        let tandem = reads_of(&twice, 300);
        let mut excised = reads_of(&once, 300);
        excised.extend(reads_of(&with_locus(&[]), 600));
        excised.extend(reads_of(&round, 300));
        for (all, restored) in [(tandem, 1), (excised, 0)] {
            let reads: Reads = all.iter().map(Vec::as_slice).collect();
            let mut drafts = [
                Draft::new(once.clone(), false),
                Draft::new(unit.clone(), true),
            ];
            let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
            let mapped = map_reads(&targets, &reads, &params.map).expect("valid settings");
            let (kept, n) =
                restore_copies(&mut drafts, &reads, &mapped, &params.map).expect("valid settings");
            assert_eq!(n, restored);
            assert_eq!(kept, [true, restored == 0]);
            assert!(
                drafts[0].seq
                    == if restored == 1 {
                        twice.clone()
                    } else {
                        once.clone()
                    }
            );
        }
    }

    #[test]
    fn versions_are_alike_up_to_where_they_differ_outright() {
        let a = random(3_000, 21);
        let mut b = a.clone();
        b[500] = if b[500] == b'A' { b'C' } else { b'A' };
        b.insert(1_000, b'G');
        // 1,500..2,000 (of a) another sequence.
        b.splice(1_501..2_001, random(500, 22));
        let (start, end) = alike_ends(&a, &b);
        assert!((1_495..=1_500).contains(&start), "start {start}");
        assert!((995..=1_000).contains(&end), "end {end}");
        assert_eq!(alike_ends(&a, &a), (3_000, 0));
    }
}
