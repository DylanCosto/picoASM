//! Undoing a layout's guess at a repeat longer than the reads that cross it.
//!
//! Where a repeat's copies are longer than most reads, the string graph cannot tell which flank
//! of one copy runs on into which flank of the other, and a contig can run through a copy from
//! the flank before one into the flank after the other (on a *S. europaeiscabiei* chromosome, a
//! 16 kb repeat in two copies). The two flanks left over then end in the repeat as drafts of
//! their own: one arrives at the repeat's start in the contig, its earlier bases those of
//! another place, and the other leaves from the repeat's end. Here each draft is searched for
//! such a pair of arrivals ([`Arrival`]), the reads are weighed on the two ways of pairing the
//! four flanks (the draft's own, and the other), and the draft is rewired where reads span the
//! other way and not the draft's own ([`rewire_at_repeats`]). With no read across either way
//! nothing is changed.

use picoasm_overlap::{InvalidOverlapParamsError, OverlapParams, Reads, map_reads};
use picoasm_seq::reverse_complement;

use crate::join::{MIN_SPAN_FLANK, MIN_SPANNING_READS, SPANNING_DOMINANCE};
use crate::polish::Draft;

/// A draft end is looked for in other drafts over its last this many bases (or half the draft).
const PROBE: usize = 30_000;

/// An end arrives at a repeat if its bases before the stretch that maps into the other draft
/// are at least this many (they are another place's)...
const MIN_DIVERGED: u32 = 1_000;

/// ...the stretch that maps is at least this long...
const MIN_MATCHED: u32 = 3_000;

/// ...and it reaches the end's last base but for at most this many.
const END_SLACK: u32 = 200;

/// A repeat is at most this long: longer, no read would span it.
const MAX_REPEAT: usize = 60_000;

/// The flanks either side of a repeat that reads are weighed on are this long.
const FLANK: usize = 2_000;

/// Which boundary of a repeat in a draft an end arrives at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Boundary {
    /// The repeat's start (the end runs along the draft forward).
    Start,
    /// The repeat's end (the end runs along the draft backward).
    End,
}

/// A linear draft's end whose last bases lie inside another draft and whose earlier ones do
/// not: it arrives at a boundary of a repeat there from another place.
#[derive(Clone, Debug)]
struct Arrival {
    /// The draft the end is of, and whether it is the draft's end (else its start).
    from: (usize, bool),
    /// The draft it arrives in.
    into: usize,
    boundary: Boundary,
    /// Where in `into` the repeat's boundary lies (its first base, or one past its last).
    at: usize,
    /// The end's draft read outward (its end last), without its shallow bases there.
    outward: Vec<u8>,
    /// Where in `outward` its bases stop matching `into` (its bases from here on lie in the
    /// repeat).
    diverged: usize,
    /// Where in `into` the end's last base lies (`at` lies on the other side of it).
    last_in: usize,
}

/// Rewires linear drafts that run through a repeat from one copy's flank into another's where
/// reads span the other pairing (see the module docs); returns the drafts and how many were
/// rewired.
///
/// # Errors
///
/// Returns an error if `map` cannot be used.
pub(crate) fn rewire_at_repeats(
    drafts: Vec<Draft>,
    reads: &Reads,
    map: &OverlapParams,
) -> Result<(Vec<Draft>, usize), InvalidOverlapParamsError> {
    let arrivals = arrivals(&drafts, map)?;
    // Each draft with one arrival at a repeat's start and one at its end, from two other
    // drafts, the start before the end and the repeat no longer than reads could span.
    let mut triples: Vec<(usize, &Arrival, &Arrival)> = Vec::new();
    for d in 0..drafts.len() {
        let here: Vec<&Arrival> = arrivals.iter().filter(|a| a.into == d).collect();
        let starts: Vec<&&Arrival> = here
            .iter()
            .filter(|a| a.boundary == Boundary::Start)
            .collect();
        let ends: Vec<&&Arrival> = here
            .iter()
            .filter(|a| a.boundary == Boundary::End)
            .collect();
        if let ([s], [e]) = (starts.as_slice(), ends.as_slice())
            && s.from.0 != e.from.0
            && s.at < e.at
            && e.at - s.at <= MAX_REPEAT
        {
            triples.push((d, s, e));
        }
    }
    // Each draft takes part in one rewiring at most.
    let mut used = vec![false; drafts.len()];
    triples.retain(|&(d, s, e)| {
        let all = [d, s.from.0, e.from.0];
        if all.iter().any(|&k| used[k]) {
            return false;
        }
        for k in all {
            used[k] = true;
        }
        true
    });
    if triples.is_empty() {
        return Ok((drafts, 0));
    }
    // Four ways across each repeat, spelled with their flanks: the draft's own, the two
    // arrivals' (which a join of their ends would make), and the two crossed.
    let mut windows: Vec<Vec<u8>> = Vec::new();
    let mut spans: Vec<(usize, usize)> = Vec::new();
    for &(d, s, e) in &triples {
        for (left, right) in [
            (Side::Draft, Side::Draft),
            (Side::Arrival, Side::Arrival),
            (Side::Draft, Side::Arrival),
            (Side::Arrival, Side::Draft),
        ] {
            let (seq, from, to) = spell(&drafts[d].seq, s, e, left, right);
            windows.push(seq);
            spans.push((from, to));
        }
    }
    let targets: Reads = windows.iter().map(Vec::as_slice).collect();
    let mut spanning = vec![0_usize; windows.len()];
    for m in map_reads(&targets, reads, map)?.into_iter().flatten() {
        let k = m.target.index();
        let (from, to) = spans[k];
        if pos(m.target_start) + MIN_SPAN_FLANK <= from
            && pos(m.target_end) >= to + MIN_SPAN_FLANK
            && to + MIN_SPAN_FLANK <= windows[k].len()
        {
            spanning[k] += 1;
        }
    }
    let mut drafts: Vec<Option<Draft>> = drafts.into_iter().map(Some).collect();
    let mut made: Vec<Draft> = Vec::new();
    for (t, &(d, s, e)) in triples.iter().enumerate() {
        let n = &spanning[4 * t..4 * t + 4];
        let (kept, crossed) = (n[0] + n[1], n[2] + n[3]);
        log::debug!(
            "bacteria: a draft of {} bases runs through a repeat of {} bases at {}; reads span \
             its own way {} and {} times, the other {} and {}",
            drafts[d].as_ref().map_or(0, |x| x.seq.len()),
            e.at - s.at,
            s.at,
            n[0],
            n[1],
            n[2],
            n[3]
        );
        // The other way needs reads across, and twice as many as the draft's own.
        if crossed < MIN_SPANNING_READS || crossed < SPANNING_DOMINANCE * kept {
            continue;
        }
        let draft = drafts[d].take().expect("each draft rewired once");
        // The arrivals' drafts are rebuilt into the two below.
        drafts[s.from.0] = None;
        drafts[e.from.0] = None;
        // The draft's bases up to where the second arrival's end lies, then that arrival's
        // draft read inward; and the first arrival's draft up to its end, then the draft's
        // bases on from where that end lies.
        let mut a = draft.seq[..e.last_in].to_vec();
        a.extend_from_slice(&reverse_complement(&e.outward));
        let mut b = s.outward.clone();
        b.extend_from_slice(&draft.seq[s.last_in..]);
        for seq in [a, b] {
            let mut rewired = Draft::new(seq, false);
            rewired.depth = draft.depth;
            made.push(rewired);
        }
        log::info!("bacteria: a draft run through a repeat into another copy's flank rewired");
    }
    let rewired = made.len() / 2;
    let mut out: Vec<Draft> = drafts.into_iter().flatten().collect();
    out.extend(made);
    Ok((out, rewired))
}

/// Which flank a way across a repeat takes at one side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Side {
    /// The draft's own.
    Draft,
    /// The arrival's (at the repeat's start, the bases before its stretch in the draft; at
    /// its end, those after).
    Arrival,
}

/// One way across the repeat between arrivals `s` (at its start) and `e` (at its end) in
/// `draft`: [`FLANK`] bases of the `left` flank, the repeat (the draft's bases), and as many of
/// the `right` flank; and where the repeat starts and ends in it.
fn spell(
    draft: &[u8],
    s: &Arrival,
    e: &Arrival,
    left: Side,
    right: Side,
) -> (Vec<u8>, usize, usize) {
    let mut seq = match left {
        Side::Draft => draft[s.at.saturating_sub(FLANK)..s.at].to_vec(),
        Side::Arrival => s.outward[s.diverged.saturating_sub(FLANK)..s.diverged].to_vec(),
    };
    let from = seq.len();
    seq.extend_from_slice(&draft[s.at..e.at]);
    let to = seq.len();
    match right {
        Side::Draft => seq.extend_from_slice(&draft[e.at..(e.at + FLANK).min(draft.len())]),
        Side::Arrival => {
            // Read inward, the arrival's bases after the repeat follow its diverging point.
            let inward = reverse_complement(&e.outward);
            let from_here = inward.len() - e.diverged;
            seq.extend_from_slice(&inward[from_here..(from_here + FLANK).min(inward.len())]);
        }
    }
    (seq, from, to)
}

/// The arrivals at repeats among `drafts` (see [`Arrival`]).
fn arrivals(
    drafts: &[Draft],
    map: &OverlapParams,
) -> Result<Vec<Arrival>, InvalidOverlapParamsError> {
    // Each linear draft's two ends, read outward, and their probes.
    let mut ends: Vec<((usize, bool), Vec<u8>)> = Vec::new();
    for (i, d) in drafts.iter().enumerate() {
        if d.circular || d.seq.len() < 2 * MIN_MATCHED as usize {
            continue;
        }
        for at_end in [false, true] {
            let full = if at_end {
                d.seq.clone()
            } else {
                reverse_complement(&d.seq)
            };
            let shallow = d.shallow_ends[usize::from(at_end)];
            let outward = full[..full.len() - shallow.min(full.len() / 2)].to_vec();
            ends.push(((i, at_end), outward));
        }
    }
    let probes: Vec<&[u8]> = ends
        .iter()
        .map(|(_, o)| &o[o.len() - PROBE.min(o.len() / 2)..])
        .collect();
    let queries: Reads = probes.iter().copied().collect();
    let mut found = Vec::new();
    for (into, d) in drafts.iter().enumerate() {
        if d.circular || d.seq.len() < 2 * FLANK + MIN_MATCHED as usize {
            continue;
        }
        let target: Reads = [d.seq.as_slice()].into_iter().collect();
        for (k, m) in map_reads(&target, &queries, map)?.into_iter().enumerate() {
            let Some(m) = m else { continue };
            let ((from, at_end), outward) = &ends[k];
            let probe_len = u32::try_from(probes[k].len()).expect("probes are short");
            if *from == into
                || m.query_end + END_SLACK < probe_len
                || m.query_start < MIN_DIVERGED
                || m.query_end - m.query_start < MIN_MATCHED
            {
                continue;
            }
            // Where the probe starts in the end's outward bases.
            let offset = outward.len() - probes[k].len();
            let diverged = offset + pos(m.query_start);
            let (boundary, at, last_in) = if m.strand.is_reverse() {
                (Boundary::End, pos(m.target_end), pos(m.target_start))
            } else {
                (Boundary::Start, pos(m.target_start), pos(m.target_end))
            };
            // Inside the draft, with room for the flanks.
            if at < FLANK || at + FLANK > d.seq.len() {
                continue;
            }
            found.push(Arrival {
                from: (*from, *at_end),
                into,
                boundary,
                at,
                outward: outward.clone(),
                diverged,
                last_in,
            });
        }
    }
    Ok(found)
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

    /// A circle of four 40 kb stretches with a 12 kb repeat after the first and the third
    /// (copies at 40-52 kb and 132-144 kb), laid out as a draft that runs from the third
    /// stretch through a copy into the second (the other copy's flank), and the two flanks
    /// left over ending in the repeat; and reads of 16 kb, two across the true junction of the
    /// copies in `spanned` (by where they start in the genome).
    fn swapped(spanned: &[usize]) -> (Vec<u8>, Vec<Draft>, Vec<Vec<u8>>) {
        let repeat = random(12_000, 21);
        let u: Vec<Vec<u8>> = (0..4).map(|i| random(40_000, 22 + i)).collect();
        let genome = [&u[0][..], &repeat, &u[1], &u[2], &repeat, &u[3]].concat();
        let drafts = vec![
            // The third stretch's last 30 kb, the repeat, the second's first 30 kb.
            Draft::new(
                [&u[2][10_000..], &repeat[..], &u[1][..30_000]].concat(),
                false,
            ),
            // The first stretch, into 7 kb of the repeat.
            Draft::new([&u[0][..], &repeat[..7_000]].concat(), false),
            // From 5 kb into the repeat, then the fourth stretch.
            Draft::new([&repeat[5_000..], &u[3][..]].concat(), false),
        ];
        let reads = spanned
            .iter()
            .flat_map(|&start| [start - 1_500, start - 1_000])
            .map(|s| {
                genome
                    .iter()
                    .cycle()
                    .skip(s)
                    .take(16_000)
                    .copied()
                    .collect()
            })
            .collect();
        (genome, drafts, reads)
    }

    #[test]
    fn a_draft_run_into_another_copys_flank_is_rewired_where_reads_span_the_copies() {
        // Reads across the second copy, from the third stretch into the fourth.
        let (genome, drafts, seqs) = swapped(&[132_000]);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = BactParams::new(ReadKind::Ont).map;
        let (drafts, made) = rewire_at_repeats(drafts, &reads, &map).expect("valid settings");
        assert_eq!((made, drafts.len()), (1, 2));
        let circle = [genome.as_slice(), &genome].concat();
        for d in &drafts {
            let found = circle.windows(d.seq.len()).any(|w| w == d.seq.as_slice())
                || circle
                    .windows(d.seq.len())
                    .any(|w| w == reverse_complement(&d.seq).as_slice());
            assert!(
                found,
                "a draft of {} bases is not in the genome",
                d.seq.len()
            );
        }
    }

    #[test]
    fn a_draft_through_a_repeat_no_read_spans_is_left_as_it_is() {
        let (_, drafts, seqs) = swapped(&[]);
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let map = BactParams::new(ReadKind::Ont).map;
        let before = drafts.clone();
        let (drafts, made) = rewire_at_repeats(drafts, &reads, &map).expect("valid settings");
        assert_eq!(made, 0);
        assert_eq!(drafts, before);
    }
}
