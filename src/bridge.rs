//! Closing a linear draft whose two ends lie inside another draft, and leaving out linear
//! drafts whose bases the other drafts hold.
//!
//! Where reads are scarce, or a stretch has a second version in part of the population (an element
//! that also excises as a circle), the string graph can end the chromosome's contig at both sides
//! of the stretch and lay the stretch out as a contig of its own, longer than the gap, with other
//! bases on. Its ends then overlap neither of the chromosome's, and the end joins miss it (the more
//! so for the few raw read bases past its dead ends, which the drafts' [`Draft::shallow_ends`] mark
//! and the end pieces here leave out). Here a linear draft at least [`CONTAINED_PERCENT`] of whose
//! bases map onto the other drafts adds nothing and is left out, the shortest first, and so is a
//! circular draft another circle of about its length holds (a high-copy plasmid laid out twice;
//! a circle's bases run on across its cut, so that a copy cut elsewhere maps whole). Then each
//! linear draft's first and last [`PIECE`] bases are mapped onto each other linear draft; where
//! both lie inside it, on one strand, in the order that closes the draft, the bases between them
//! are put in and the draft is a circle ([`bridge`]); and the draft that held them is left out in
//! turn.

use std::collections::HashSet;

use crate::overlap::{EditOp, InvalidOverlapParamsError, Overlap, OverlapParams, Reads, map_reads};
use crate::seq::{Anchor, SketchParams, Sketcher, reverse_complement};
use rayon::prelude::*;

use crate::polish::Draft;

/// A draft's end pieces are this long (or a quarter of it): long enough to pass repeats (an
/// rRNA operon is 5 kb), short enough that the stretch next to the gap decides where a piece
/// maps, not more of the draft that the other draft also holds elsewhere.
const PIECE: usize = 10_000;

/// An end piece maps up to the draft's end when at most this many of its bases are left over
/// there, and maps whole when at most this many are left over in all (not along a repeat
/// inside it).
const MAX_END_CLIP: u32 = 500;

/// At most this many bases are put in between a draft's ends.
const MAX_BRIDGE: u32 = 200_000;

/// A linear draft is left out when this share (in percent) of its bases maps onto the others.
const CONTAINED_PERCENT: u64 = 90;

/// A draft is mapped onto the others in pieces this long...
const CHUNK: usize = 5_000;

/// ...and, once the drafts are closed, in pieces this long: a line stitched from stretches the
/// circles hold (a variant arrangement in part of the cells) has a junction every few
/// kilobases, and a piece across one maps on one side only (*Burkholderia* lines of 10-23 kb,
/// 70-90% held in 5 kb pieces).
pub(crate) const FINE_CHUNK: usize = 2_000;

/// A draft's piece is held where it maps with substitutions at no more than this share of its
/// columns: a plasmid of a family (*Borrelia*'s cp32s, 4-13% apart, 2.3-12% in substitutions)
/// maps onto its relatives along most of its length and was left out as held by them. Indels
/// are not counted: a copy polished from few or short reads (a *V. parahaemolyticus* with
/// reads of 1.2 kb) differs from the deeper draft at 1-2% in indels and about 1% in
/// substitutions, and kept as not held, its ends were rivals to the right joins.
const MAX_HELD_DIVERGENCE: f64 = 0.02;

/// A circular draft is a copy of a circle whose length is within this share (in percent) of
/// its own, if that circle holds it.
const MAX_COPY_LEN_PERCENT: usize = 10;

/// Closes linear drafts whose ends lie in another (see the module docs), then leaves out those
/// the others hold; returns the drafts, how many were closed and how many left out.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn bridge(
    drafts: Vec<Draft>,
    map: &OverlapParams,
) -> Result<(Vec<Draft>, usize, usize), InvalidOverlapParamsError> {
    let (mut drafts, before) = without_contained(drafts, map, CHUNK)?;
    let held_seeds = all_seeds(drafts.iter().map(|d| d.seq.as_slice()), map)?;
    let mut order: Vec<usize> = (0..drafts.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(drafts[i].seq.len()));
    let mut closed = 0;
    for &i in &order {
        if drafts[i].circular {
            continue;
        }
        let [head, tail] = end_pieces(drafts[i].deep_part()).map(|p| seeds(p, map));
        let (head, tail) = (head?, tail?);
        // A draft that shares no seed with an end piece holds none of it.
        let holds_both =
            |j: usize| !held_seeds[j].is_disjoint(&head) && !held_seeds[j].is_disjoint(&tail);
        let mut closing = None;
        for j in (0..drafts.len()).filter(|&j| j != i && !drafts[j].circular && holds_both(j)) {
            closing = closed_by(drafts[i].deep_part(), &drafts[j].seq, map)?;
            if closing.is_some() {
                break;
            }
        }
        if closing.is_some() && !ends_unique(&drafts, i, map)? {
            log::info!(
                "bacteria: a draft of {} bases not closed by another's bases: an end lies in a \
                 repeat",
                drafts[i].seq.len()
            );
            closing = None;
        }
        if let Some(seq) = closing {
            log::info!(
                "bacteria: a draft of {} bases closed into a circle of {} by another's bases",
                drafts[i].seq.len(),
                seq.len()
            );
            drafts[i] = Draft::new(seq, true);
            closed += 1;
        }
    }
    let (drafts, after) = without_contained(drafts, map, CHUNK)?;
    Ok((drafts, closed, before + after))
}

/// Returns `seq` closed into a circle by the bases between its two ends in `other`, if both
/// its end pieces map inside `other`, on one strand, up to its ends, in the order that closes
/// it.
fn closed_by(
    seq: &[u8],
    other: &[u8],
    map: &OverlapParams,
) -> Result<Option<Vec<u8>>, InvalidOverlapParamsError> {
    let pieces = end_pieces(seq);
    let n = pieces[0].len();
    if n == 0 {
        return Ok(None);
    }
    let len = u32::try_from(seq.len()).expect("drafts are shorter than 4 Gb");
    let piece = u32::try_from(n).expect("pieces are short");
    let pieces: Reads = pieces.into_iter().collect();
    let targets: Reads = [other].into_iter().collect();
    let mapped = map_reads(&targets, &pieces, map)?;
    let (Some(head), Some(tail)) = (&mapped[0], &mapped[1]) else {
        return Ok(None);
    };
    let whole = |m: &Overlap| m.query_end - m.query_start + MAX_END_CLIP >= piece;
    if head.strand != tail.strand
        || head.query_start > MAX_END_CLIP
        || tail.query_end + MAX_END_CLIP < piece
        || !whole(head)
        || !whole(tail)
    {
        return Ok(None);
    }
    // The draft's bases from its head's first aligned base to its tail's last.
    let kept = &seq[head.query_start as usize..(len - piece + tail.query_end) as usize];
    let between = between(head, tail, other);
    Ok(between.map(|b| [kept, &b].concat()))
}

/// Whether both end pieces of draft `i` map whole onto one place only among the drafts (see
/// [`end_unique`]).
fn ends_unique(
    drafts: &[Draft],
    i: usize,
    map: &OverlapParams,
) -> Result<bool, InvalidOverlapParamsError> {
    Ok(end_unique(drafts, i, false, map)? && end_unique(drafts, i, true, map)?)
}

/// Whether the first (or, `at_end`, the last) [`PIECE`] bases of draft `i`'s polished bases map
/// whole onto one place only among the drafts: onto one draft other than `i`, or onto the rest
/// of `i`, but not both nor two of them. An end that lies in a copy of a repeat longer than the
/// piece maps whole onto every copy, and which copy's continuation it runs on into is then a
/// guess: on a *C. difficile* chromosome with a 34.5 kb repeat in two copies, which no read
/// spans, the guess cut 1.6 Mb out of it as a circle, and on a *S. europaeiscabiei* one with a
/// 16 kb repeat, 1.8 Mb. A place counts only where the piece differs from it in substitutions
/// at no more than [`MAX_HELD_DIVERGENCE`] of its columns: a *Borrelia* cp32 plasmid's end maps
/// whole onto its relatives, 4-5% apart, and a cp32 laid out round 1.4 times stayed a 42 kb
/// line, closed neither through its end anchors nor along its own ends.
pub(crate) fn end_unique(
    drafts: &[Draft],
    i: usize,
    at_end: bool,
    map: &OverlapParams,
) -> Result<bool, InvalidOverlapParamsError> {
    let seq = drafts[i].deep_part();
    let [head, tail] = end_pieces(seq);
    let n = head.len();
    let piece = u32::try_from(n).expect("pieces are short");
    // The piece, and the rest of its draft without it.
    let (end, rest) = if at_end {
        (tail, &seq[..seq.len() - n])
    } else {
        (head, &seq[n..])
    };
    let queries: Reads = [end].into_iter().collect();
    let targets = std::iter::once(rest).chain(
        (0..drafts.len())
            .filter(|&j| j != i)
            .map(|j| drafts[j].seq.as_slice()),
    );
    let mut places = 0;
    for target in targets {
        let one: Reads = [target].into_iter().collect();
        if let Some(m) = &map_reads(&one, &queries, map)?[0]
            && m.query_end - m.query_start + MAX_END_CLIP >= piece
            && within(m, MAX_HELD_DIVERGENCE)
        {
            places += 1;
        }
    }
    log::debug!(
        "bacteria: the {} piece of a draft of {} bases maps whole onto {places} places",
        if at_end { "last" } else { "first" },
        seq.len()
    );
    Ok(places <= 1)
}

/// A draft's first and last [`PIECE`] bases (or a quarter of it).
fn end_pieces(seq: &[u8]) -> [&[u8]; 2] {
    let n = PIECE.min(seq.len() / 4);
    [&seq[..n], &seq[seq.len() - n..]]
}

/// The seeds `map` looks up for `seq`. A piece and a target that share none cannot map onto
/// each other, so a target sharing none with any piece can be left out of an index unseen: it
/// adds no occurrence of the pieces' seeds, and the pieces map just as they would with it.
fn seeds(seq: &[u8], map: &OverlapParams) -> Result<HashSet<u64>, InvalidOverlapParamsError> {
    let mut sketcher = Sketcher::new(SketchParams::new(map.kmer_len, map.density, 1)?);
    let mut anchors: Vec<Anchor> = Vec::new();
    sketcher
        .sketch(seq, &mut anchors)
        .expect("drafts are short enough to sketch");
    Ok(anchors.iter().map(|a| a.hash).collect())
}

/// [`seeds`] of each of `seqs`.
fn all_seeds<'seq>(
    seqs: impl Iterator<Item = &'seq [u8]>,
    map: &OverlapParams,
) -> Result<Vec<HashSet<u64>>, InvalidOverlapParamsError> {
    let seqs: Vec<&[u8]> = seqs.collect();
    seqs.par_iter().map(|s| seeds(s, map)).collect()
}

/// The bases of `other` that run on from `tail`'s alignment to `head`'s, in the draft's
/// orientation; `None` if they are in the wrong order or too many.
fn between(head: &Overlap, tail: &Overlap, other: &[u8]) -> Option<Vec<u8>> {
    let (from, to) = if head.strand.is_reverse() {
        // The draft runs backwards along `other`: past its tail's alignment lie the bases
        // before it, down to the end of its head's.
        (head.target_end, tail.target_start)
    } else {
        (tail.target_end, head.target_start)
    };
    if from > to || to - from > MAX_BRIDGE {
        return None;
    }
    let bases = &other[from as usize..to as usize];
    Some(if head.strand.is_reverse() {
        reverse_complement(bases)
    } else {
        bases.to_vec()
    })
}

/// Leaves out linear drafts (not the longest draft) of which at least [`CONTAINED_PERCENT`]
/// maps onto the drafts kept, and circular ones that a kept circle of about their length
/// ([`MAX_COPY_LEN_PERCENT`]) holds so (a high-copy plasmid laid out more than once), the
/// shortest (then shallowest) first, each mapped in pieces of `chunk` bases (see
/// [`held_in_pieces`]); returns the drafts kept and how many were left out.
pub(crate) fn without_contained(
    drafts: Vec<Draft>,
    map: &OverlapParams,
    chunk: usize,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let kept = not_contained(&drafts, map, chunk)?;
    let left_out = kept.iter().filter(|&&k| !k).count();
    let drafts = drafts
        .into_iter()
        .zip(kept)
        .filter_map(|(d, k)| k.then_some(d))
        .collect();
    Ok((drafts, left_out))
}

/// Which of `drafts` [`without_contained`] keeps.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn not_contained(
    drafts: &[Draft],
    map: &OverlapParams,
    chunk: usize,
) -> Result<Vec<bool>, InvalidOverlapParamsError> {
    let longest = (0..drafts.len()).max_by_key(|&i| (drafts[i].seq.len(), std::cmp::Reverse(i)));
    let held_seeds = all_seeds(drafts.iter().map(Draft::deep_part), map)?;
    // A circle's bases run on across its cut, so that a piece of a copy cut elsewhere maps.
    let held_bases: Vec<Vec<u8>> = drafts
        .iter()
        .map(|d| {
            let deep = d.deep_part();
            if d.circular {
                [deep, &deep[..CHUNK.min(deep.len())]].concat()
            } else {
                deep.to_vec()
            }
        })
        .collect();
    let mut order: Vec<usize> = (0..drafts.len()).collect();
    order.sort_by(|&a, &b| {
        (drafts[a].seq.len(), drafts[a].depth)
            .partial_cmp(&(drafts[b].seq.len(), drafts[b].depth))
            .expect("depths are not NaN")
    });
    let mut kept = vec![true; drafts.len()];
    for &i in &order {
        if Some(i) == longest {
            continue;
        }
        let circular = drafts[i].circular;
        let copy_of = |j: usize| {
            let (a, b) = (drafts[i].seq.len(), drafts[j].seq.len());
            drafts[j].circular && a.abs_diff(b) * 100 <= a.max(b) * MAX_COPY_LEN_PERCENT
        };
        // The others' polished bases only: a draft's raw end can hold a stretch loosely that
        // only this one holds well. Those sharing no seed with it hold none of it.
        let own = seeds(&drafts[i].seq, map)?;
        let targets: Reads = (0..drafts.len())
            .filter(|&j| j != i && kept[j] && !held_seeds[j].is_disjoint(&own))
            .filter(|&j| !circular || copy_of(j))
            .map(|j| held_bases[j].as_slice())
            .collect();
        let shape = if circular {
            Shape::Circular
        } else {
            Shape::Linear
        };
        if held_in_pieces(
            &drafts[i].seq,
            shape,
            &targets,
            map,
            chunk,
            Some(MAX_HELD_DIVERGENCE),
        )? {
            log::info!(
                "bacteria: a {} draft of {} bases left out: the others hold it",
                if circular { "circular" } else { "linear" },
                drafts[i].seq.len()
            );
            kept[i] = false;
        }
    }
    Ok(kept)
}

/// Whether a sequence is a line or a circle (a circle's pieces run on across its cut).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Shape {
    Linear,
    Circular,
}

/// Whether at least [`CONTAINED_PERCENT`] of `seq` maps onto `targets`, in pieces of
/// [`CHUNK`] bases.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn held(
    seq: &[u8],
    shape: Shape,
    targets: &Reads,
    map: &OverlapParams,
) -> Result<bool, InvalidOverlapParamsError> {
    held_in_pieces(seq, shape, targets, map, CHUNK, Some(MAX_HELD_DIVERGENCE))
}

/// Whether at least [`CONTAINED_PERCENT`] of `seq`'s bases map onto `targets`, in pieces of
/// `chunk` bases that start every half piece (short ones find a variant of a target with a
/// stretch left out or put in; overlapping ones leave no end too short to map). A circle's
/// pieces run on across its cut, and the bases they cover are counted round it. A piece counts
/// only where it differs at no more than `max_divergence` of its columns, if given (polished
/// drafts against each other: [`MAX_HELD_DIVERGENCE`]).
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn held_in_pieces(
    seq: &[u8],
    shape: Shape,
    targets: &Reads,
    map: &OverlapParams,
    chunk: usize,
    max_divergence: Option<f64>,
) -> Result<bool, InvalidOverlapParamsError> {
    if targets.is_empty() || seq.is_empty() {
        return Ok(false);
    }
    let len = seq.len();
    let chunk = chunk.min(len);
    let step = (chunk / 2).max(1);
    let run_on: Vec<u8>;
    let (bases, starts): (&[u8], Vec<usize>) = match shape {
        // Every half piece to the end; a last piece under half a piece long is left out (the
        // one before covers it).
        Shape::Linear => (
            seq,
            (0..len)
                .step_by(step)
                .filter(|&s| s == 0 || len - s >= step)
                .collect(),
        ),
        Shape::Circular => {
            run_on = [seq, &seq[..chunk]].concat();
            (&run_on, (0..len).step_by(step).collect())
        }
    };
    let pieces: Reads = starts
        .iter()
        .map(|&s| &bases[s..(s + chunk).min(bases.len())])
        .collect();
    // The bases of `seq` some piece maps over.
    let mut covered = vec![false; len];
    for (s, m) in starts.iter().zip(map_reads(targets, &pieces, map)?) {
        let close = |m: &Overlap| max_divergence.is_none_or(|d| within(m, d));
        if let Some(m) = m.filter(close) {
            for q in m.query_start..m.query_end {
                covered[(s + q as usize) % len] = true;
            }
        }
    }
    let held = covered.iter().filter(|&&c| c).count();
    Ok(held as u64 * 100 >= len as u64 * CONTAINED_PERCENT)
}

/// The columns of `o` where both sides have a base and the bases differ.
/// Whether `o` differs in substitutions at no more than `max_divergence` of its columns.
fn within(o: &Overlap, max_divergence: f64) -> bool {
    #[expect(clippy::cast_precision_loss, reason = "a share, not an exact count")]
    let close = substitutions(o) as f64 <= max_divergence * o.columns() as f64;
    close
}

fn substitutions(o: &Overlap) -> u64 {
    o.edits
        .iter()
        .filter(|e| e.op() == EditOp::Mismatch)
        .map(|e| u64::from(e.len()))
        .sum()
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

    fn is_rotation(a: &[u8], b: &[u8]) -> bool {
        a.len() == b.len() && [b, b].concat().windows(a.len()).any(|w| w == a)
    }

    #[test]
    fn a_draft_whose_ends_lie_in_another_is_closed_and_the_other_left_out() {
        // A circle of 200 kb; the chromosome's draft lacks 30 kb of it, which a draft of
        // 80 kb (the gap and 25 kb either side) holds, on either strand.
        let genome = random(200_000, 5);
        let chromosome: Vec<u8> = genome[100_000..]
            .iter()
            .chain(&genome[..70_000])
            .copied()
            .collect();
        for bridging in [
            genome[45_000..125_000].to_vec(),
            reverse_complement(&genome[45_000..125_000]),
        ] {
            let drafts = vec![
                Draft::new(chromosome.clone(), false),
                Draft::new(bridging, false),
            ];
            let map = BactParams::new(ReadKind::Ont).map;
            let (drafts, closed, left_out) = bridge(drafts, &map).expect("valid settings");
            assert_eq!((closed, left_out), (1, 1));
            assert_eq!(drafts.len(), 1);
            assert!(drafts[0].circular);
            assert!(is_rotation(&drafts[0].seq, &genome));
        }
    }

    #[test]
    fn a_draft_whose_end_lies_in_a_repeat_is_not_closed() {
        // As above, on the forward strand, but 30 kb about the chromosome draft's start (where
        // it maps into the bridging draft, and where the bridging draft's end maps into it) lie
        // in a third draft too, with other bases around: a copy of a repeat, so which copy those
        // ends run on from is a guess, and neither draft is closed through the other.
        let genome = random(200_000, 5);
        let chromosome: Vec<u8> = genome[100_000..]
            .iter()
            .chain(&genome[..70_000])
            .copied()
            .collect();
        let copy = [
            random(30_000, 11),
            genome[98_000..128_000].to_vec(),
            random(30_000, 12),
        ]
        .concat();
        let drafts = vec![
            Draft::new(chromosome.clone(), false),
            Draft::new(genome[45_000..125_000].to_vec(), false),
            Draft::new(copy, false),
        ];
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed, _) = bridge(drafts, &map).expect("valid settings");
        assert_eq!(closed, 0);
        assert!(drafts.iter().all(|d| !d.circular));
        assert!(drafts.iter().any(|d| d.seq == chromosome));
    }

    #[test]
    fn drafts_with_ends_elsewhere_are_left_as_they_are() {
        let (a, b) = (random(100_000, 6), random(60_000, 7));
        let drafts = vec![Draft::new(a.clone(), false), Draft::new(b.clone(), false)];
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed, left_out) = bridge(drafts, &map).expect("valid settings");
        assert_eq!((closed, left_out), (0, 0));
        assert_eq!(drafts[0].seq, a);
        assert_eq!(drafts[1].seq, b);
    }

    #[test]
    fn copies_of_a_small_circle_cut_elsewhere_are_left_out() {
        // A chromosome, a plasmid of 3.8 kb laid out three times: as a circle, as the circle
        // cut elsewhere on the other strand, and as a line a little longer than it, cut
        // elsewhere again; and a second plasmid of its length, another sequence.
        let chromosome = random(100_000, 8);
        let plasmid = random(3_800, 9);
        let other = random(3_800, 10);
        let cut = |at: usize| -> Vec<u8> { [&plasmid[at..], &plasmid[..at]].concat() };
        let mut drafts = vec![
            Draft::new(chromosome, true),
            Draft::new(plasmid.clone(), true),
            Draft::new(reverse_complement(&cut(1_500)), true),
            Draft::new([cut(2_500), plasmid[2_500..2_620].to_vec()].concat(), false),
            Draft::new(other.clone(), true),
        ];
        drafts[1].depth = 900.0;
        drafts[2].depth = 600.0;
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, closed, left_out) = bridge(drafts, &map).expect("valid settings");
        assert_eq!((closed, left_out), (0, 2));
        let seqs: Vec<&[u8]> = drafts.iter().map(|d| d.seq.as_slice()).collect();
        assert_eq!(seqs.len(), 3);
        assert!(seqs.contains(&plasmid.as_slice()));
        assert!(seqs.contains(&other.as_slice()));
    }
}
