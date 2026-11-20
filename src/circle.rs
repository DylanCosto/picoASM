//! Linear drafts whose end runs on into their own start.
//!
//! The joins close a draft whose end overlaps its start ([`crate::join`]), but only where the
//! overlap aligns from end to end, or holds a 10 kb stretch both ends have alike. Where a segment
//! inside the overlap is inverted in one copy and not the other (an invertible element, in both
//! orientations in the population), neither holds, and the draft stays linear with the overlap
//! twice: an *A. mictus* chromosome came out 11.3 kb too long, its two ends each holding the
//! stretch, a 2.8 kb segment inside it turned round in one. (Reads ran from the end's unique bases
//! 9 kb into the stretch, up to that segment, but none across all of it: the copy the end held was
//! not needed, but whether the chromosome is a circle the reads did not say, and it was left
//! linear.) Here each linear draft's first [`PIECE`] bases are looked for near its end and its last
//! near its start, at places that agree on the overlap's length, and most of the overlap must map
//! onto the start ([`MIN_MATCHED_SHARE`]), and neither end may lie in a copy another draft holds
//! too ([`crate::bridge::end_unique`]). The end's copy is then left out. If reads run from the
//! end's unique bases across the overlap into the start's ([`MIN_SPAN_FLANK`] either side), the
//! draft is a circle. So it is if [`MIN_JUNCTION_READS`] run from the end's unique bases into
//! the start's copy ([`MIN_SPAN_FLANK`] either side of the junction): a linear replicon with the
//! stretch at both ends would have it there the same way round, but bacteria's linear replicons
//! end in inverted repeats (and hairpins, left to [`crate::hairpin`]), which the forward overlap
//! here does not take. The *A. mictus* chromosome had 0 reads across the whole overlap and 133
//! across the junction (an independent read-based check, 24 with 2 kb either side); it is a
//! circle. With neither, the draft stays linear, its duplicate copy still left out.

use crate::overlap::{InvalidOverlapParamsError, Overlap, OverlapParams, Reads, map_reads};

use crate::join::{MIN_SPAN_FLANK, MIN_SPANNING_READS};
use crate::polish::Draft;

/// The pieces of a draft's ends looked for in each other are this long.
const PIECE: usize = 2_000;

/// A draft's end is looked for in its first this many bases (or a third of it), and its start
/// in its last as many.
const SEARCH: usize = 50_000;

/// An overlap is at least this long: shorter ones the joins settle.
const MIN_OVERLAP: usize = 3_000;

/// The two pieces' places may imply overlaps this far apart (or a tenth of the overlap): the
/// copies can differ by indels.
const MIN_LENGTH_SLACK: usize = 1_000;

/// At least this share (in percent) of the end's copy, in [`PIECE`]s, maps onto the start's at
/// the overlap's offset: two copies of a repeat at the ends would agree on the pieces at the
/// very ends only.
const MIN_MATCHED_SHARE: usize = 50;

/// The unique bases either side of the overlap that reads are weighed on.
const FLANK: usize = 2_000;

/// The draft's last piece is looked for this far in from its end, the first that maps whole
/// (bases the two copies differ in can keep the very last from mapping).
const TAIL_SKIPS: [usize; 5] = [0, 1_000, 2_000, 3_000, 4_000];

/// A piece maps whole if at most this many of its bases are left out.
const PIECE_CLIP: u32 = 200;

/// A draft is a circle if at least this many reads run across its junction alone (see the
/// module docs).
const MIN_JUNCTION_READS: usize = 3;

/// What became of the linear drafts whose ends overlap.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Closed {
    /// Closed into circles: reads run across.
    pub(crate) circles: usize,
    /// Left linear, the end's copy of the overlap left out: no read runs across.
    pub(crate) trimmed: usize,
}

/// Closes, or trims, the linear drafts whose end runs on into their own start (see the module
/// docs).
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn close_self_overlaps(
    drafts: Vec<Draft>,
    reads: &Reads,
    map: &OverlapParams,
) -> Result<(Vec<Draft>, Closed), InvalidOverlapParamsError> {
    let mut closed = Closed::default();
    // Each draft's replacement, if its ends overlap.
    let mut replaced: Vec<Option<Draft>> = vec![None; drafts.len()];
    for (i, draft) in drafts.iter().enumerate() {
        if draft.circular || draft.hairpins != [false; 2] {
            continue;
        }
        let deep = draft.deep_part();
        let Some(overlap) = overlap_of(deep, map)? else {
            continue;
        };
        // A linear replicon read round a hairpin end runs back along itself, and its ends can
        // overlap: *Borrelia* lp21, laid out 62.8 kb long folded at both ends, its ends 8.2 kb
        // alike, was closed into a 54.9 kb circle instead of cut to its 27.5 kb.
        if crate::hairpin::folds_back(&draft.seq, map)? {
            continue;
        }
        // Ends in copies of a repeat another draft holds too run on into whichever copy.
        if !(crate::bridge::end_unique(&drafts, i, false, map)?
            && crate::bridge::end_unique(&drafts, i, true, map)?)
        {
            continue;
        }
        let kept = deep[..deep.len() - overlap].to_vec();
        let Spans { across, junction } = spanning(&kept, overlap, reads, map)?;
        let circular = across >= MIN_SPANNING_READS || junction >= MIN_JUNCTION_READS;
        log::info!(
            "bacteria: a linear draft of {} bases whose ends overlap by {overlap} bases: {} \
             ({across} reads across, {junction} across the junction)",
            draft.seq.len(),
            if circular {
                "closed"
            } else {
                "trimmed, left linear"
            }
        );
        if circular {
            closed.circles += 1;
        } else {
            closed.trimmed += 1;
        }
        let mut done = Draft::new(kept, circular);
        done.depth = draft.depth;
        replaced[i] = Some(done);
    }
    let out = drafts
        .into_iter()
        .zip(replaced)
        .map(|(draft, done)| done.unwrap_or(draft))
        .collect();
    Ok((out, closed))
}

/// How many bases `seq`'s end repeats of its start, if its first [`PIECE`] lie near its end and
/// its last near its start at places that agree, and most of the overlap matches.
fn overlap_of(seq: &[u8], map: &OverlapParams) -> Result<Option<usize>, InvalidOverlapParamsError> {
    let len = seq.len();
    let window = SEARCH.min(len / 3);
    if window < MIN_OVERLAP + PIECE {
        return Ok(None);
    }
    let (start, end) = (&seq[..window], &seq[len - window..]);
    // Where the draft's first bases lie in its end: the overlap runs from there to the end.
    let Some(head) = whole_forward(&seq[..PIECE], end, map)? else {
        return Ok(None);
    };
    let overlap = window - pos(head.target_start);
    if overlap < MIN_OVERLAP {
        return Ok(None);
    }
    // Its last bases lie in the start's copy of the overlap: on the other strand, if the copies
    // have a segment there each its own way round (the *A. mictus* chromosome's last 4 kb).
    let slack = MIN_LENGTH_SLACK.max(overlap / 10);
    let mut agrees = false;
    for skip in TAIL_SKIPS {
        let to = len - skip;
        if let Some(tail) = whole(&seq[to - PIECE..to], start, map)? {
            agrees = pos(tail.target_end) <= overlap + slack;
            break;
        }
    }
    if !agrees {
        return Ok(None);
    }
    // Most of the end's copy, piece by piece, at about its place in the start's.
    let copy = &seq[len - overlap..];
    let pieces: Vec<usize> = (0..copy.len() / PIECE).map(|k| k * PIECE).collect();
    let queries: Reads = pieces.iter().map(|&at| &copy[at..at + PIECE]).collect();
    let target: Reads = [&seq[..(overlap + slack).min(len)]].into_iter().collect();
    let matched = map_reads(&target, &queries, map)?
        .iter()
        .zip(&pieces)
        .filter(|(m, at)| {
            m.as_ref().is_some_and(|m| {
                !m.strand.is_reverse()
                    && pos(m.target_start).abs_diff(**at + pos(m.query_start)) <= slack
            })
        })
        .count();
    Ok((matched * 100 >= MIN_MATCHED_SHARE * pieces.len().max(1)).then_some(overlap))
}

/// `piece`'s mapping onto `region`, if it maps whole and forward.
fn whole_forward(
    piece: &[u8],
    region: &[u8],
    map: &OverlapParams,
) -> Result<Option<Overlap>, InvalidOverlapParamsError> {
    Ok(whole(piece, region, map)?.filter(|m| !m.strand.is_reverse()))
}

/// `piece`'s mapping onto `region`, if it maps whole (on either strand).
fn whole(
    piece: &[u8],
    region: &[u8],
    map: &OverlapParams,
) -> Result<Option<Overlap>, InvalidOverlapParamsError> {
    let targets: Reads = [region].into_iter().collect();
    let queries: Reads = [piece].into_iter().collect();
    let len = u32::try_from(piece.len()).expect("pieces are short");
    Ok(map_reads(&targets, &queries, map)?
        .pop()
        .flatten()
        .filter(|m| m.query_end - m.query_start + PIECE_CLIP >= len))
}

/// How many reads run across a circle's junction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Spans {
    /// From the unique bases before the stretch held twice to those past it.
    across: usize,
    /// From the unique bases before the stretch into it.
    junction: usize,
}

/// How many `reads` run across the junction of `circle` (its end into its start, the first
/// `overlap` bases being the stretch the draft held twice), from [`MIN_SPAN_FLANK`] of the
/// unique bases before it to as many past the stretch, or as many into it.
fn spanning(
    circle: &[u8],
    overlap: usize,
    reads: &Reads,
    map: &OverlapParams,
) -> Result<Spans, InvalidOverlapParamsError> {
    let before = FLANK.min(circle.len().saturating_sub(overlap));
    let after = (overlap + FLANK).min(circle.len());
    let window = [&circle[circle.len() - before..], &circle[..after]].concat();
    let targets: Reads = [window.as_slice()].into_iter().collect();
    let (from, to) = (before, before + overlap);
    let mut spans = Spans {
        across: 0,
        junction: 0,
    };
    for m in map_reads(&targets, reads, map)?.iter().flatten() {
        if pos(m.target_start) + MIN_SPAN_FLANK > from {
            continue;
        }
        if pos(m.target_end) >= from + MIN_SPAN_FLANK {
            spans.junction += 1;
        }
        if pos(m.target_end) >= to + MIN_SPAN_FLANK && to + MIN_SPAN_FLANK <= window.len() {
            spans.across += 1;
        }
    }
    Ok(spans)
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use crate::seq::reverse_complement;

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

    fn is_rotation(a: &[u8], b: &[u8]) -> bool {
        a.len() == b.len() && [b, b].concat().windows(a.len()).any(|w| w == a)
    }

    /// A 200 kb circle laid out as a line holding its first 10 kb again at its end, with 2 kb
    /// inside that copy turned round; and reads of 15 kb, and two of 18-19 kb across the
    /// junction if `across`.
    fn overlapping(across: bool) -> (Vec<u8>, Draft, Vec<Vec<u8>>) {
        let genome = random(200_000, 31);
        let mut copy = genome[..10_000].to_vec();
        let turned = reverse_complement(&copy[4_000..6_000]);
        copy[4_000..6_000].copy_from_slice(&turned);
        let line = [genome.as_slice(), &copy].concat();
        let mut reads: Vec<Vec<u8>> = (0..40)
            .map(|i| genome[20_000 + i * 4_000..35_000 + i * 4_000].to_vec())
            .collect();
        if across {
            // From 4-5 kb before the end, across the junction, to 14 kb past it.
            for start in [195_000, 196_000] {
                reads.push([&genome[start..], &genome[..14_000]].concat());
            }
        }
        (genome, Draft::new(line, false), reads)
    }

    #[test]
    fn a_draft_whose_end_runs_into_its_start_is_closed_where_reads_cross() {
        let (genome, draft, seqs) = overlapping(true);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed) = close_self_overlaps(vec![draft], &reads, &map).expect("valid");
        assert_eq!(
            closed,
            Closed {
                circles: 1,
                trimmed: 0
            }
        );
        assert!(drafts[0].circular);
        assert!(is_rotation(&drafts[0].seq, &genome));
    }

    #[test]
    fn reads_across_the_junction_alone_close_the_draft() {
        // Three reads from 4-6 kb before the end to 3 kb into the stretch held twice: none
        // across the whole stretch.
        let (genome, draft, mut seqs) = overlapping(false);
        for start in [194_000, 195_000, 196_000] {
            seqs.push([&genome[start..], &genome[..3_000]].concat());
        }
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed) = close_self_overlaps(vec![draft], &reads, &map).expect("valid");
        assert_eq!(
            closed,
            Closed {
                circles: 1,
                trimmed: 0
            }
        );
        assert!(drafts[0].circular);
        assert!(is_rotation(&drafts[0].seq, &genome));
    }

    #[test]
    fn with_no_read_across_the_duplicate_is_left_out_and_the_draft_stays_linear() {
        let (genome, draft, seqs) = overlapping(false);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed) = close_self_overlaps(vec![draft], &reads, &map).expect("valid");
        assert_eq!(
            closed,
            Closed {
                circles: 0,
                trimmed: 1
            }
        );
        assert!(!drafts[0].circular);
        assert_eq!(drafts[0].seq, genome);
    }

    #[test]
    fn ends_in_a_repeat_another_draft_holds_are_left_as_they_are() {
        // The same 12 kb at both ends of a draft, and in another draft with other bases on:
        // which copy the ends run on into is a guess.
        let repeat = random(12_000, 32);
        let line = [&repeat[..], &random(150_000, 33), &repeat[..]].concat();
        let other = [&random(40_000, 34), &repeat[..], &random(40_000, 35)].concat();
        let drafts = vec![Draft::new(line.clone(), false), Draft::new(other, false)];
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed) = close_self_overlaps(drafts, &Reads::new(), &map).expect("valid");
        assert_eq!(closed, Closed::default());
        assert_eq!(drafts[0].seq, line);
    }

    /// `seq` with about 1 in 32 bases changed to another.
    fn relative(seq: &[u8], seed: u64) -> Vec<u8> {
        let mut state = seed;
        seq.iter()
            .map(|&b| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                match (state >> 59, b) {
                    (1.., _) => b,
                    (0, b'A') => b'C',
                    (0, b'C') => b'G',
                    (0, b'G') => b'T',
                    (0, _) => b'A',
                }
            })
            .collect()
    }

    #[test]
    fn ends_relatives_hold_some_percent_apart_still_close() {
        // Two of a plasmid family's relatives, 3% apart in substitutions, each holding the
        // draft's first 20 kb: not copies the ends could run on into.
        let (genome, draft, seqs) = overlapping(true);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = BactParams::new(ReadKind::Ont).map;
        let drafts = vec![
            draft,
            Draft::new(relative(&genome[..20_000], 36), false),
            Draft::new(relative(&genome[..20_000], 38), false),
        ];
        let (drafts, closed) = close_self_overlaps(drafts, &reads, &map).expect("valid");
        assert_eq!(closed.circles, 1);
        assert!(is_rotation(&drafts[0].seq, &genome));
    }

    #[test]
    fn a_line_folded_back_at_its_hairpin_ends_is_left_to_be_cut_there() {
        // A 30 kb linear replicon laid out from 8 kb into its folded start, across it, and back
        // along all but its first 500 bases: the ends hold 7.5 kb alike, on the same strand.
        let replicon = random(30_000, 37);
        let line = [
            reverse_complement(&replicon[..8_000]),
            replicon.clone(),
            reverse_complement(&replicon[500..]),
        ]
        .concat();
        let map = BactParams::new(ReadKind::Ont).map;
        assert!(overlap_of(&line, &map).expect("valid").is_some());
        let drafts = vec![Draft::new(line.clone(), false)];
        let (drafts, closed) = close_self_overlaps(drafts, &Reads::new(), &map).expect("valid");
        assert_eq!(closed, Closed::default());
        assert_eq!(drafts[0].seq, line);
    }
}
