//! Linear replicons whose ends are hairpins (*Borrelia*'s chromosome and linear plasmids, the
//! prophage plasmids of N15-like phages), laid out with their ends folded back.
//!
//! At a hairpin end the two strands are one molecule, so a read can run round the end and back
//! along the other strand. The layout follows such reads: a linear plasmid with hairpins at
//! both ends becomes a circle twice its length, the plasmid and then its reverse complement;
//! and a linear contig that ends in a hairpin goes on past it with its own reverse complement.
//! Here a circle that is its own reverse complement turned round (pieces at [`PIECES`] places
//! round it, each reverse complemented, map whole back onto the rest of it, all mirrored about
//! one point)
//! is cut at its two mirror points to the one replicon, linear with hairpins at both ends; and
//! a linear draft whose start (or end) is mirrored about a point [`FOLD_PIECE`] or more in
//! (the bases either side of the point mirror each other) is cut there, with a hairpin at that
//! end. A circle of a replicon and its reverse complement is also what a hairpin-ended
//! replicon replicates through (an inverted dimer, which a telomere resolvase cuts in two), and
//! reads cannot tell the two apart: the one replicon is what is given.

use picoasm_overlap::{InvalidOverlapParamsError, Overlap, OverlapParams, Reads, map_reads};
use picoasm_seq::reverse_complement;

use crate::polish::Draft;

/// A circle is tested with this many pieces, spread evenly round it...
const PIECES: usize = 8;

/// ...of which all but this many must map whole (one across a mirror point cannot), all
/// mirrored about one point.
const MISSES: usize = 2;

/// ...each this long (or a quarter of the circle).
const PIECE: usize = 5_000;

/// The mirror points the pieces give agree within this many bases.
const MIRROR_SLACK: usize = 200;

/// A piece maps whole if at most this many of its bases are left out...
const CLIP: u32 = 300;

/// ...or, for the bases either side of a fold, this many: next to the fold lie the hairpin's
/// loop and the errors of reads where they turn (*Borrelia* lp21: 306 bases).
const FOLD_CLIP: u32 = 500;

/// A linear draft's end is tested with a piece this long: a fold shorter than it is left as
/// it is (the reads that turned there are a shallow end).
const FOLD_PIECE: usize = 2_000;

/// A linear draft's end is looked for in its first (or last) this many bases.
const FOLD_REACH: usize = 200_000;

/// Cuts the drafts laid out round hairpin ends (see the module docs); returns the drafts and
/// how many were cut.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn unfold(
    drafts: Vec<Draft>,
    map: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let mut drafts = drafts;
    let mut cut = 0;
    for draft in &mut drafts {
        if draft.circular {
            if let Some(start) = mirror_of_circle(&draft.seq, map)? {
                let half = draft.seq.len() / 2;
                log::info!(
                    "bacteria: a circle of {} bases is a linear replicon of {half} with hairpin \
                     ends, read round them",
                    draft.seq.len()
                );
                draft.rotate(start);
                draft.seq.truncate(half);
                draft.circular = false;
                draft.hairpins = [true, true];
                draft.shallow_ends = [0, 0];
                cut += 1;
            }
            continue;
        }
        // The whole draft: the folded copy is held by the few reads that turned round, and
        // so lies in its shallow end.
        if let Some(at) = fold_at_start(&draft.seq, map)? {
            log::info!(
                "bacteria: a linear draft folds back on itself {at} bases in: a hairpin end"
            );
            draft.seq.drain(..at);
            draft.shallow_ends[0] = 0;
            draft.hairpins[0] = true;
            cut += 1;
        }
        let back = reverse_complement(&draft.seq);
        if let Some(at) = fold_at_start(&back, map)? {
            log::info!(
                "bacteria: a linear draft folds back on itself {at} bases from its end: a hairpin end"
            );
            let keep = draft.seq.len().saturating_sub(at);
            draft.seq.truncate(keep);
            draft.shallow_ends[1] = 0;
            draft.hairpins[1] = true;
            cut += 1;
        }
    }
    Ok((drafts, cut))
}

/// Whether the linear draft `seq` folds back on itself at either end: [`unfold`] cuts it there.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn folds_back(
    seq: &[u8],
    map: &OverlapParams,
) -> Result<bool, InvalidOverlapParamsError> {
    Ok(fold_at_start(seq, map)?.is_some()
        || fold_at_start(&reverse_complement(seq), map)?.is_some())
}

/// Where the circle `seq` is to start for its first half to be one strand of a linear replicon
/// with hairpin ends (see the module docs), if it is one.
fn mirror_of_circle(
    seq: &[u8],
    map: &OverlapParams,
) -> Result<Option<usize>, InvalidOverlapParamsError> {
    let n = seq.len();
    let piece = PIECE.min(n / PIECES);
    if piece < FOLD_PIECE {
        return Ok(None);
    }
    let mut mirrors = Vec::with_capacity(PIECES);
    for i in 0..PIECES {
        let s = i * n / PIECES;
        // The circle from just past the piece round to its start: the piece itself is out, so
        // its reverse complement can map only where the circle mirrors it.
        let rest: Vec<u8> = (s + piece..s + n).map(|x| seq[x % n]).collect();
        let own: Vec<u8> = (s..s + piece).map(|x| seq[x % n]).collect();
        let Some(m) = map_one(&rest, &reverse_complement(&own), map)?.filter(|m| whole(m, piece))
        else {
            continue;
        };
        // Base `s + piece - 1 - j` reverse complemented lies at `rest`'s `target_start + j`,
        // that is at `s + piece + target_start + j` round the circle: each base `x` is
        // mirrored at `c - x`.
        let at = s + piece + pos(m.target_start) - pos(m.query_start);
        mirrors.push((s + piece - 1 + at) % n);
    }
    if mirrors.len() + MISSES < PIECES {
        return Ok(None);
    }
    let first = mirrors[0];
    let close = |c: usize| {
        let d = c.abs_diff(first);
        d.min(n - d) <= MIRROR_SLACK
    };
    if !mirrors.iter().all(|&c| close(c)) {
        return Ok(None);
    }
    // The mirror points are at c/2 and c/2 + n/2; the replicon runs from one to the other.
    Ok(Some(sharpen(seq, (first / 2 + 1) % n, Shape::Circle)))
}

/// How many bases at the start of `seq` are the reverse complement of those after them (the
/// draft runs back along itself from a hairpin there), if at least [`FOLD_PIECE`].
fn fold_at_start(
    seq: &[u8],
    map: &OverlapParams,
) -> Result<Option<usize>, InvalidOverlapParamsError> {
    let reach = seq.len().min(FOLD_REACH);
    if reach < 3 * FOLD_PIECE {
        return Ok(None);
    }
    let head = reverse_complement(&seq[..FOLD_PIECE]);
    let after = &seq[FOLD_PIECE..reach];
    let Some(m) = map_one(after, &head, map)?.filter(|m| whole(m, FOLD_PIECE)) else {
        return Ok(None);
    };
    // The first base is mirrored at `FOLD_PIECE + target_end - 1`: the fold lies half way.
    let mirror = FOLD_PIECE
        + pos(m.target_end)
        + pos(u32::try_from(FOLD_PIECE).expect("short") - m.query_end);
    let fold = sharpen(seq, mirror.div_ceil(2), Shape::Line);
    if fold < FOLD_PIECE || fold + FOLD_PIECE > reach {
        return Ok(None);
    }
    // The bases either side of the fold mirror each other (an inverted repeat apart does not).
    let before = reverse_complement(&seq[fold - FOLD_PIECE..fold]);
    let next = &seq[fold..(fold + 2 * FOLD_PIECE).min(reach)];
    let piece = u32::try_from(FOLD_PIECE).expect("pieces are short");
    match map_one(next, &before, map)? {
        Some(m)
            if m.query_end - m.query_start + FOLD_CLIP >= piece
                && pos(m.target_start) <= pos(FOLD_CLIP) =>
        {
            Ok(Some(fold))
        }
        _ => Ok(None),
    }
}

/// The fold points are moved to where, within this many bases, the bases either side mirror
/// each other best: a fold from the reads' alignments is off by the indels near it (*Borrelia*
/// lp21, cut 152 bases past its telomere).
const SHARPEN: usize = 500;

/// The bases either side of a fold are compared over this many.
const MIRROR_SPAN: usize = 100;

/// Whether a sequence runs on round its end (to sharpen a fold on it).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Shape {
    Line,
    Circle,
}

/// The point within [`SHARPEN`] of `at` in `seq` where the [`MIRROR_SPAN`] bases before it are
/// most like the reverse complement of those after it (the nearest such, of equals).
fn sharpen(seq: &[u8], at: usize, shape: Shape) -> usize {
    let n = seq.len();
    let base = |i: isize| -> Option<u8> {
        match shape {
            Shape::Circle => {
                Some(seq[usize::try_from(i.rem_euclid(isize::try_from(n).ok()?)).ok()?])
            }
            Shape::Line => usize::try_from(i).ok().and_then(|i| seq.get(i).copied()),
        }
    };
    let complement = |b: u8| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    };
    let score = |p: isize| -> usize {
        (0..isize::try_from(MIRROR_SPAN).expect("short"))
            .filter(|&j| match (base(p - 1 - j), base(p + j)) {
                (Some(a), Some(b)) => a == complement(b),
                _ => false,
            })
            .count()
    };
    let at = isize::try_from(at).expect("positions fit in isize");
    let reach = isize::try_from(SHARPEN).expect("short");
    let best = (at - reach..=at + reach)
        .filter(|&p| shape == Shape::Circle || (p >= 0 && p <= isize::try_from(n).expect("fits")))
        .max_by_key(|&p| (score(p), std::cmp::Reverse((p - at).abs())))
        .unwrap_or(at);
    match shape {
        Shape::Circle => usize::try_from(best.rem_euclid(isize::try_from(n).expect("fits")))
            .expect("non-negative"),
        Shape::Line => usize::try_from(best).expect("non-negative"),
    }
}

/// `query` mapped onto `target` alone, forward.
fn map_one(
    target: &[u8],
    query: &[u8],
    map: &OverlapParams,
) -> Result<Option<Overlap>, InvalidOverlapParamsError> {
    let targets: Reads = [target].into_iter().collect();
    let queries: Reads = [query].into_iter().collect();
    Ok(map_reads(&targets, &queries, map)?
        .pop()
        .flatten()
        .filter(|m| !m.strand.is_reverse()))
}

/// Whether `m` maps a piece of `len` bases whole, forward.
fn whole(m: &Overlap, len: usize) -> bool {
    let len = u32::try_from(len).expect("pieces are short");
    !m.strand.is_reverse() && m.query_end - m.query_start + CLIP >= len
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
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

    #[test]
    fn a_circle_of_a_replicon_and_its_reverse_complement_is_cut_to_it() {
        // A linear plasmid of 30 kb read round its hairpins, as a circle of 60 kb turned to
        // start 7 kb in; and a plain circle of 60 kb, left alone.
        let plasmid = random(30_000, 31);
        let doubled = [plasmid.as_slice(), &reverse_complement(&plasmid)].concat();
        let turned = [&doubled[7_000..], &doubled[..7_000]].concat();
        let map = BactParams::new(ReadKind::Ont).map;
        let drafts = vec![
            Draft::new(turned, true),
            Draft::new(random(60_000, 32), true),
        ];
        let (drafts, cut) = unfold(drafts, &map).expect("valid");
        assert_eq!(cut, 1);
        assert!(!drafts[0].circular && drafts[0].hairpins == [true, true]);
        let len = drafts[0].seq.len();
        assert!(len.abs_diff(30_000) <= 2, "{len}");
        // The plasmid, on one strand or the other, but for a base or two at its ends.
        let found = &drafts[0].seq;
        let middle = &plasmid[100..29_900];
        let holds = |s: &[u8]| s.windows(middle.len()).any(|w| w == middle);
        assert!(holds(found) || holds(&reverse_complement(found)));
        assert!(drafts[1].circular);
    }

    #[test]
    fn a_linear_draft_folded_back_at_its_ends_is_cut_at_the_folds() {
        // A chromosome of 200 kb with 20 kb folded back at its start and 15 kb at its end;
        // and one with an inverted repeat 30 kb apart near its start, left alone.
        let chromosome = random(200_000, 33);
        let start_fold = reverse_complement(&chromosome[..20_000]);
        let end_fold = reverse_complement(&chromosome[200_000 - 15_000..]);
        let seq = [start_fold.as_slice(), &chromosome, &end_fold].concat();
        let repeat = random(3_000, 34);
        let with_repeat = [
            &repeat[..],
            &random(30_000, 35),
            &reverse_complement(&repeat),
            &random(170_000, 36),
        ]
        .concat();
        let map = BactParams::new(ReadKind::Ont).map;
        let drafts = vec![
            Draft::new(seq, false),
            Draft::new(with_repeat.clone(), false),
        ];
        let (drafts, cut) = unfold(drafts, &map).expect("valid");
        assert_eq!(cut, 2);
        assert_eq!(drafts[0].hairpins, [true, true]);
        let len = drafts[0].seq.len();
        assert!(len.abs_diff(200_000) <= 2, "{len}");
        assert_eq!(drafts[1].seq, with_repeat);
        assert_eq!(drafts[1].hairpins, [false, false]);
    }
}
