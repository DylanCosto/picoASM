//! Replicons shorter than their reads: small plasmids, from reads that hold them whole.
//!
//! A plasmid of a few kilobases has reads as long as itself: with a rapid (transposase) kit
//! each read is the circle cut once, the whole of it from a random place, and a longer read
//! runs on past its own start. Such reads overlap each other at both ends, and the longer ones
//! overlap themselves, which the string graph does not lay out. Here the reads that the drafts
//! do not hold are taken as candidate units of circles: a read that repeats its start further
//! on (k-mers of its first [`PROBE_SPAN`] bases found again at one distance) is cut to one
//! turn, and any other read is taken whole. Units are grouped by their overlaps; a group of at
//! least [`MIN_UNITS`] units of about one length, whose units start all round the circle when
//! mapped onto one of them written twice ([`MIN_START_QUARTERS`]: a linear molecule's reads,
//! such as a control DNA's, all start at its end), becomes a circular draft, which polishing
//! corrects, if enough reads hold it whole ([`MIN_DEPTH_SHARE`] of the chromosome's depth)
//! and the drafts do not hold it in pieces ([`VARIANT_PIECE`]): a handful of reads is a
//! variant of a replicon in a few of its copies, or poor reads of one that missed it. A
//! circular draft of the same circle ([`same_circle`]) that the reads hardly support (a garbled
//! draft of it) is replaced.

use std::collections::HashMap;

use crate::overlap::{
    InvalidOverlapParamsError, OverlapParams, ReadIndex, Reads, find_overlaps, map_reads,
};
use rayon::prelude::*;

use crate::polish::Draft;

/// K-mers this long find a read's start again.
const K: usize = 16;

/// The probes are the k-mers starting in this many bases at the read's start...
const PROBE_SPAN: usize = 500;

/// ...and at least this share of them (in percent) must recur at one distance (give or take
/// [`PERIOD_BIN`] bases).
const MIN_SUPPORT_PERCENT: usize = 30;

/// Distances are counted in bins of this many bases.
const PERIOD_BIN: usize = 16;

/// A circle is at least this long...
const MIN_UNIT: usize = 500;

/// ...and at most this long (larger plasmids have reads shorter than themselves and lay out).
const MAX_UNIT: usize = 30_000;

/// A circle is made from at least this many reads' units...
const MIN_UNITS: usize = 5;

/// ...whose starts, mapped onto one of them written twice, fall in at least this many of its
/// four quarters.
const MIN_START_QUARTERS: usize = 3;

/// A unit's length is within this share (one over it) of its group's median.
const UNIT_LEN_DIVISOR: usize = 10;

/// A read is held by the drafts when at least this share of it (in percent) maps onto them.
const HELD_PERCENT: usize = 90;

/// A circle is checked against the drafts in pieces this long.
const VARIANT_PIECE: usize = 1_000;

/// A circle needs at least this share of the chromosome's depth (the longest draft's) in reads
/// that hold it whole.
const MIN_DEPTH_SHARE: f64 = 0.25;

/// At most this many units are grouped (a high-copy plasmid gives thousands).
const MAX_CANDIDATES: usize = 3_000;

/// Circles within this share (one over it) of each other's length...
const SAME_LEN_DIVISOR: usize = 4;

/// ...that share at least this share (in percent) of the shorter one's k-mers, on either
/// strand, are one circle.
const SAME_KMER_PERCENT: usize = 30;

/// Whether circle `a` (a read, with its errors) is circle `b`: of about its length, sharing
/// [`SAME_KMER_PERCENT`] of their k-mers whatever their starts and strands. (A raw read of a
/// short circle with adapter bases at its ends maps onto it too partially to be held.)
fn same_circle(a: &[u8], b: &[u8]) -> bool {
    if a.len().abs_diff(b.len()) * SAME_LEN_DIVISOR > a.len().max(b.len()) {
        return false;
    }
    let canonical = |s: &[u8]| -> std::collections::HashSet<Vec<u8>> {
        s.windows(K)
            .map(|k| {
                let rc = crate::seq::reverse_complement(k);
                if rc.as_slice() < k { rc } else { k.to_vec() }
            })
            .collect()
    };
    let (ka, kb) = (canonical(a), canonical(b));
    let shared = ka.intersection(&kb).count();
    shared * 100 >= ka.len().min(kb.len()) * SAME_KMER_PERCENT
}

/// Returns the length of one turn of `read` round a circle, if it repeats its start that far
/// on (see the module docs).
fn period(read: &[u8]) -> Option<usize> {
    if read.len() < MIN_UNIT + PROBE_SPAN {
        return None;
    }
    let probes: HashMap<&[u8], usize> = read[..PROBE_SPAN + K - 1]
        .windows(K)
        .enumerate()
        .map(|(i, k)| (k, i))
        .collect();
    let mut votes: HashMap<usize, usize> = HashMap::new();
    for (j, k) in read.windows(K).enumerate().skip(MIN_UNIT) {
        if let Some(&i) = probes.get(k) {
            let d = j - i;
            if (MIN_UNIT..=MAX_UNIT).contains(&d) {
                *votes.entry(d / PERIOD_BIN).or_default() += 1;
            }
        }
    }
    // The best-supported bin (with its neighbour, as a period near a bin's edge splits).
    let (&bin, _) = votes.iter().max_by_key(|&(&b, &n)| {
        (
            n + votes.get(&(b + 1)).copied().unwrap_or(0),
            std::cmp::Reverse(b),
        )
    })?;
    let support = votes[&bin] + votes.get(&(bin + 1)).copied().unwrap_or(0);
    if support * 100 < PROBE_SPAN * MIN_SUPPORT_PERCENT {
        return None;
    }
    // The exact period: the commonest distance in the two bins.
    let mut exact: HashMap<usize, usize> = HashMap::new();
    for (j, k) in read.windows(K).enumerate().skip(MIN_UNIT) {
        if let Some(&i) = probes.get(k) {
            let d = j - i;
            if (bin * PERIOD_BIN..(bin + 2) * PERIOD_BIN).contains(&d) {
                *exact.entry(d).or_default() += 1;
            }
        }
    }
    exact
        .into_iter()
        .max_by_key(|&(d, n)| (n, std::cmp::Reverse(d)))
        .map(|(d, _)| d)
}

/// A small circle found in reads, and the draft it replaces (a garbled one of it), if any.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Circle {
    pub(crate) replaces: Option<usize>,
    pub(crate) draft: Draft,
}

/// Returns circular drafts of the small replicons whose reads `drafts` do not hold (see the
/// module docs).
///
/// # Errors
///
/// Returns an error if `map` or `overlap` cannot be used.
pub(crate) fn circles_in_reads(
    reads: &Reads,
    drafts: &[Draft],
    map: &OverlapParams,
    overlap: &OverlapParams,
) -> Result<Vec<Circle>, InvalidOverlapParamsError> {
    // The reads of a small circle's length, and whether the drafts hold them: each circular
    // draft with its start again past its end, so a read across its cut maps whole.
    let candidates: Vec<ReadIndex> = reads
        .indices()
        .filter(|&r| (MIN_UNIT..=2 * MAX_UNIT).contains(&(reads.read_len(r) as usize)))
        .collect();
    let held_by: Vec<Vec<u8>> = drafts
        .iter()
        .map(|d| {
            let mut seq = d.seq.clone();
            if d.circular {
                seq.extend_from_slice(&d.seq[..d.seq.len().min(2 * MAX_UNIT)]);
            }
            seq
        })
        .collect();
    let targets: Reads = held_by.iter().map(Vec::as_slice).collect();
    let chromosome_depth = drafts
        .iter()
        .max_by_key(|d| d.seq.len())
        .map_or(0.0, |d| d.depth);
    let queries: Reads = candidates.iter().map(|&r| reads.get(r)).collect();
    let free: Vec<ReadIndex> = candidates
        .iter()
        .zip(map_reads(&targets, &queries, map)?)
        .filter(|(r, m)| {
            m.as_ref().is_none_or(|m| {
                (m.query_end - m.query_start) as usize * 100
                    < reads.read_len(**r) as usize * HELD_PERCENT
            })
        })
        .map(|(&r, _)| r)
        .collect();
    // One turn of each read that runs on past its start; the others whole.
    let mut units: Vec<Vec<u8>> = free
        .par_iter()
        .map(|&r| {
            let read = reads.get(r);
            match period(read) {
                Some(p) => read[..p].to_vec(),
                None => read.to_vec(),
            }
        })
        .filter(|u| (MIN_UNIT..=MAX_UNIT).contains(&u.len()))
        .collect();
    log::debug!(
        "bacteria: {} reads the drafts do not hold, {} units of small circles",
        free.len(),
        units.len()
    );
    if units.len() < MIN_UNITS {
        return Ok(Vec::new());
    }
    // Every so many, if there are too many: they are in read order, so evenly from the run.
    if units.len() > MAX_CANDIDATES {
        let step = units.len().div_ceil(MAX_CANDIDATES);
        units = units.into_iter().step_by(step).collect();
    }
    // Groups: units joined by an overlap, the largest group first.
    let unit_reads: Reads = units.iter().map(Vec::as_slice).collect();
    let mut parent: Vec<usize> = (0..units.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    // How many units each overlaps: a read with fewer errors overlaps more.
    let mut degree = vec![0_usize; units.len()];
    for o in find_overlaps(&unit_reads, overlap)? {
        degree[o.query.index()] += 1;
        degree[o.target.index()] += 1;
        let (a, b) = (
            root(&mut parent, o.query.index()),
            root(&mut parent, o.target.index()),
        );
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..units.len() {
        let r = root(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut groups: Vec<Vec<usize>> = groups
        .into_values()
        .filter(|g| g.len() >= MIN_UNITS)
        .collect();
    groups.sort_by_key(|g| (std::cmp::Reverse(g.len()), g[0]));
    let mut circles = Vec::new();
    for group in groups {
        // The units of about the median length; the one that overlaps most is the circle's
        // draft (the fewest errors, so it maps onto a draft that holds it already).
        let mut lens: Vec<usize> = group.iter().map(|&i| units[i].len()).collect();
        lens.sort_unstable();
        let median = lens[lens.len() / 2];
        let members: Vec<usize> = group
            .into_iter()
            .filter(|&i| units[i].len().abs_diff(median) * UNIT_LEN_DIVISOR <= median)
            .collect();
        if members.len() < MIN_UNITS {
            continue;
        }
        let seed = *members
            .iter()
            .max_by_key(|&&i| (degree[i], std::cmp::Reverse(i)))
            .expect("members are at least MIN_UNITS");
        let twice = [units[seed].as_slice(), units[seed].as_slice()].concat();
        let target: Reads = [twice.as_slice()].into_iter().collect();
        let others: Reads = members.iter().map(|&i| units[i].as_slice()).collect();
        let len = u32::try_from(units[seed].len()).expect("units are at most MAX_UNIT long");
        let mut quarters = [false; 4];
        let mut whole = 0;
        for m in map_reads(&target, &others, map)?.into_iter().flatten() {
            let mapped = (m.query_end - m.query_start) as usize;
            if mapped * 100 >= others.read_len(m.query) as usize * HELD_PERCENT {
                whole += 1;
                // Where on the circle the unit starts.
                let start = if m.strand.is_reverse() {
                    m.target_end
                } else {
                    m.target_start
                } % len;
                quarters[(4 * start / len.max(1)) as usize % 4] = true;
            }
        }
        let spread = quarters.iter().filter(|&&q| q).count();
        log::debug!(
            "bacteria: a group of {} units of {median} bases: {whole} map whole, starting in \
             {spread} quarters",
            members.len()
        );
        // As many reads hold it whole as a replicon the layout missed has (at a share of the
        // chromosome's depth); a handful are a variant of a replicon the drafts hold (a
        // plasmid with a stretch left out, in a few of its copies), or poor reads of one that
        // missed it, and a circle the drafts hold in pieces is one too.
        #[expect(clippy::cast_precision_loss, reason = "a count compared with a depth")]
        let deep_enough = whole as f64 >= MIN_DEPTH_SHARE * chromosome_depth;
        // A circular draft it is the same circle as: one the reads support stays and this is
        // its copy; one they hardly do (a garbled draft of it) is replaced.
        let same = (0..drafts.len())
            .find(|&d| drafts[d].circular && same_circle(&units[seed], &drafts[d].seq));
        let copy = same.is_some_and(|d| drafts[d].depth >= MIN_DEPTH_SHARE * chromosome_depth);
        if whole >= MIN_UNITS
            && spread >= MIN_START_QUARTERS
            && deep_enough
            && !copy
            && (same.is_some()
                || !crate::bridge::held_in_pieces(
                    &units[seed],
                    crate::bridge::Shape::Circular,
                    &targets,
                    map,
                    VARIANT_PIECE,
                    // A raw read against polished drafts: any divergence the mapping takes.
                    None,
                )?)
        {
            log::info!(
                "bacteria: a circle of {} bases from {whole} reads that hold it whole",
                len
            );
            circles.push(Circle {
                replaces: same,
                draft: Draft::new(units[seed].clone(), true),
            });
        }
    }
    Ok(circles)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn a_read_that_runs_round_a_circle_gives_its_length() {
        let circle = random(3_011, 1);
        let read: Vec<u8> = circle
            .iter()
            .cycle()
            .skip(700)
            .take(4_200)
            .copied()
            .collect();
        assert_eq!(period(&read), Some(3_011));
        assert_eq!(period(&random(5_000, 2)), None);
    }

    #[test]
    fn a_circle_is_itself_from_any_start_and_strand() {
        let circle = random(1_240, 7);
        let turned: Vec<u8> = circle
            .iter()
            .cycle()
            .skip(500)
            .take(1_240)
            .copied()
            .collect();
        assert!(same_circle(
            &crate::seq::reverse_complement(&turned),
            &circle
        ));
        assert!(!same_circle(&random(1_240, 8), &circle));
        assert!(!same_circle(
            &circle,
            &[circle.as_slice(), circle.as_slice()].concat()
        ));
    }

    #[test]
    fn a_small_plasmid_only_in_reads_longer_than_it_becomes_a_circle() {
        let plasmid = random(3_011, 3);
        let chromosome = random(50_000, 4);
        let mut seqs: Vec<Vec<u8>> = (0..12)
            .map(|i| {
                plasmid
                    .iter()
                    .cycle()
                    .skip(i * 251)
                    .take(3_600 + i * 300)
                    .copied()
                    .collect()
            })
            .collect();
        // Chromosome reads, which repeat nothing.
        seqs.extend((0..20).map(|i| chromosome[i * 2_000..i * 2_000 + 8_000].to_vec()));
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let drafts = vec![Draft::new(chromosome.clone(), true)];
        let params = crate::BactParams::new(crate::ReadKind::Ont);
        let (map, overlap) = (params.map, params.overlap);
        let circles = circles_in_reads(&reads, &drafts, &map, &overlap).expect("valid settings");
        assert_eq!(circles.len(), 1);
        assert!(circles[0].draft.circular);
        assert!(is_rotation(&circles[0].draft.seq, &plasmid));
        // Held by a draft already: nothing new.
        let drafts = vec![Draft::new(chromosome, true), Draft::new(plasmid, true)];
        assert!(
            circles_in_reads(&reads, &drafts, &map, &overlap)
                .expect("valid")
                .is_empty()
        );
    }

    #[test]
    fn reads_each_one_turn_of_a_circle_become_it_and_a_linear_molecules_do_not() {
        // A 3 kb plasmid cut once per read at random places, and a 3.5 kb linear molecule
        // whose reads all start at its end.
        let plasmid = random(3_011, 5);
        let linear = random(3_500, 6);
        let mut seqs: Vec<Vec<u8>> = (0..10)
            .map(|i| {
                plasmid
                    .iter()
                    .cycle()
                    .skip(i * 293)
                    .take(3_011)
                    .copied()
                    .collect()
            })
            .collect();
        seqs.extend((0..10).map(|i| linear[..3_400 + i * 10].to_vec()));
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let params = crate::BactParams::new(crate::ReadKind::Ont);
        let circles =
            circles_in_reads(&reads, &[], &params.map, &params.overlap).expect("valid settings");
        assert_eq!(circles.len(), 1);
        assert!(is_rotation(&circles[0].draft.seq, &plasmid));
    }
}
