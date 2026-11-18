//! Polishing finished contigs with short reads (Illumina).
//!
//! Long-read assemblies keep a few errors that every long-read polisher shares: a homopolymer a
//! base short, a base the nanopore misreads the same way on both strands. Short reads do not
//! share them, so where an isolate also has short reads they settle each base. Each read is
//! aligned to every place it fits best (a read from identical copies of a repeat, to each), and
//! each aligned column votes for what the read has there (a base, a deletion, or a base and the
//! bases inserted after it; indels moved to their leftmost place, so reads agree on where one
//! in a homopolymer is). A position takes another version when at least
//! [`ShortParams::valid_share`] of its reads show that version and at most
//! [`ShortParams::invalid_share`] show what the contig has, with at least
//! [`ShortParams::min_depth`] reads in all: Polypolish's rule and defaults (Wick & Holt 2022,
//! PLoS Comput Biol 18:e1009802).
//!
//! Unlike Polypolish, a read does not also vote where it fits a little worse: the copies of a
//! repeat that truly differ at a base (an rRNA operon among seven) would then be outvoted at
//! that base by the reads of the others. So an error in a copy of a repeat longer than the
//! reads, whose reads fit the other copies better, is left as it is.

use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU32, Ordering};

use picoasm_overlap::Reads;
use picoasm_seq::reverse_complement_into;
use rayon::prelude::*;

use crate::assemble::Contig;

/// Seeds are exact matches of this many bases.
const K: usize = 21;

/// A read is seeded every this many bases (and at its end).
const SEED_STEP: usize = 12;

/// Seeds that occur more often than this in the contigs are skipped (a repeat this many times
/// over is still found through its less common seeds, or not at all).
const MAX_OCCURRENCES: usize = 200;

/// An alignment may leave its seeds' diagonal by this many bases (indels).
const BAND: usize = 12;

/// Alignment scores: as BWA-MEM's defaults, with a linear gap.
const MATCH: i32 = 1;
const MISMATCH: i32 = -4;
const GAP: i32 = -6;
/// Leaving a read's end unaligned costs this much (a low-quality tail, an adapter).
const CLIP: i32 = -5;

/// An alignment counts when it scores at least the read's length over this.
const MIN_SCORE_DIVISOR: usize = 2;

/// The columns this near either end of an alignment do not vote (an indel there is not placed
/// reliably), nor those of a homopolymer or short tandem repeat reaching that near.
const EDGE: usize = 5;

/// The columns this near a clipped end do not vote: the clip may stand for an indel the read
/// ends too soon after to show.
const CLIPPED_EDGE: usize = 20;

/// Repeat units up to this long count as short tandem repeats at an alignment's ends.
const MAX_UNIT: usize = 3;

/// Circular contigs are indexed this far past their end, so that reads across the end align.
const WRAP: usize = 1_000;

/// When a position takes another version.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShortParams {
    /// At least this many reads must cover the position.
    pub min_depth: u32,
    /// At least this share of them must show the new version.
    pub valid_share: f64,
    /// At most this share of them may show the contig's version.
    pub invalid_share: f64,
}

impl ShortParams {
    /// Polypolish's defaults: 5 reads, 50% for the new version, 20% at most for the old.
    #[must_use]
    pub fn new() -> Self {
        Self {
            min_depth: 5,
            valid_share: 0.5,
            invalid_share: 0.2,
        }
    }
}

impl Default for ShortParams {
    fn default() -> Self {
        Self::new()
    }
}

/// What [`polish_short`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ShortStats {
    /// Reads aligned at least once.
    pub aligned: usize,
    /// Contig positions with at least [`ShortParams::min_depth`] reads.
    pub covered: usize,
    /// All contig positions.
    pub positions: usize,
    /// Positions given another version.
    pub changed: usize,
}

/// What one read shows at one contig position: its base there (or none, deleted), and the bases
/// it has between that position and the next.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Token {
    base: Option<u8>,
    inserted: Vec<u8>,
}

/// The votes for each position of one contig: plain bases and deletions counted in place, the
/// rare versions with inserted bases in a map.
struct Votes {
    /// Per position, the reads with A, C, G, T or a deletion there and nothing inserted after.
    plain: Vec<[AtomicU32; 5]>,
    inserted: Mutex<HashMap<(usize, Token), u32>>,
}

impl Votes {
    fn new(len: usize) -> Self {
        Self {
            plain: (0..len).map(|_| Default::default()).collect(),
            inserted: Mutex::new(HashMap::new()),
        }
    }

    fn add(&self, pos: usize, token: Token) {
        if token.inserted.is_empty() {
            let slot = token.base.map_or(4, code);
            if slot < 5 {
                self.plain[pos][slot].fetch_add(1, Ordering::Relaxed);
            }
        } else {
            *self
                .inserted
                .lock()
                .expect("lock poisoned")
                .entry((pos, token))
                .or_default() += 1;
        }
    }
}

/// The 2-bit code of an upper-case base; 4 for anything else.
fn code(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

const BASES: [u8; 4] = *b"ACGT";

/// The k-mer at the start of `seq` packed 2 bits a base, or `None` if it holds an `N`.
fn pack(seq: &[u8]) -> Option<u64> {
    seq.iter().try_fold(0u64, |acc, &b| {
        let c = code(b);
        (c < 4).then_some((acc << 2) | c as u64)
    })
}

/// Every k-mer of the contigs (circular ones run on [`WRAP`] bases past their end), sorted, with
/// the contig and position it starts at.
struct Index {
    entries: Vec<(u64, u32, u32)>,
}

impl Index {
    fn new(seqs: &[Vec<u8>]) -> Self {
        let mut entries: Vec<(u64, u32, u32)> = seqs
            .par_iter()
            .enumerate()
            .flat_map_iter(|(c, seq)| {
                let c = u32::try_from(c).expect("fewer than 4 billion contigs");
                seq.windows(K).enumerate().filter_map(move |(p, w)| {
                    let p = u32::try_from(p).expect("contigs are shorter than 4 Gb");
                    pack(w).map(|kmer| (kmer, c, p))
                })
            })
            .collect();
        entries.par_sort_unstable();
        Self { entries }
    }

    fn lookup(&self, kmer: u64) -> &[(u64, u32, u32)] {
        let start = self.entries.partition_point(|e| e.0 < kmer);
        let end = start + self.entries[start..].partition_point(|e| e.0 == kmer);
        &self.entries[start..end]
    }
}

/// One alignment of a read to a contig (its indexed sequence): where the read's aligned part
/// starts on the contig, and per contig position from there what the read shows.
struct Alignment {
    score: i32,
    contig: usize,
    start: usize,
    tokens: Vec<Token>,
    /// Whether the read's first and last bases are left unaligned.
    clipped: [bool; 2],
}

impl Alignment {
    /// The tokens that vote: not near an end ([`EDGE`], [`CLIPPED_EDGE`]), nor in a homopolymer
    /// or short tandem repeat the read may end inside (whose length it then does not show), nor
    /// on the base before such a repeat at its far end (where a left-aligned indel in it lies).
    /// `contig` holds the contig's bases from the alignment's start.
    fn voting(&self, contig: &[u8]) -> std::ops::Range<usize> {
        let edge = |clipped: bool| if clipped { CLIPPED_EDGE } else { EDGE };
        let n = self.tokens.len();
        let mut lo = edge(self.clipped[0]);
        let mut hi = n.saturating_sub(edge(self.clipped[1]));
        while lo < hi && repeats_across(contig, lo) {
            lo += 1;
        }
        if lo < hi && repeats_across(contig, hi) {
            while hi > lo && repeats_across(contig, hi) {
                hi -= 1;
            }
            hi = hi.saturating_sub(1);
        }
        lo..hi.max(lo)
    }
}

/// Whether a unit of up to [`MAX_UNIT`] bases repeats across the boundary before `at`.
fn repeats_across(seq: &[u8], at: usize) -> bool {
    (1..=MAX_UNIT).any(|q| at >= q && at + q <= seq.len() && seq[at - q..at] == seq[at..at + q])
}

/// Aligns `read` (both of its strands) to every place its seeds point to in `seqs` and returns
/// the best-scoring ones, each starting inside its contig.
fn align_read(read: &[u8], seqs: &[Vec<u8>], lens: &[usize], index: &Index) -> Vec<Alignment> {
    if read.len() < K {
        return Vec::new();
    }
    let mut rc = Vec::with_capacity(read.len());
    reverse_complement_into(read, &mut rc);
    let min_score = i32::try_from(read.len() / MIN_SCORE_DIVISOR).unwrap_or(i32::MAX);
    let mut found: Vec<Alignment> = Vec::new();
    for strand in [read, rc.as_slice()] {
        let mut diagonals: Vec<(u32, i64)> = Vec::new();
        let last = strand.len() - K;
        let offsets = (0..=last).step_by(SEED_STEP).chain(std::iter::once(last));
        for offset in offsets {
            let Some(kmer) = pack(&strand[offset..offset + K]) else {
                continue;
            };
            let hits = index.lookup(kmer);
            if hits.len() > MAX_OCCURRENCES {
                continue;
            }
            diagonals.extend(
                hits.iter()
                    .map(|&(_, c, p)| (c, i64::from(p) - offset as i64)),
            );
        }
        diagonals.sort_unstable();
        diagonals.dedup();
        let mut previous: Option<(u32, i64)> = None;
        for (c, d) in diagonals {
            // Diagonals a few bases apart are one place (an indel between seeds).
            if previous.is_some_and(|(pc, pd)| pc == c && d - pd <= BAND as i64) {
                continue;
            }
            previous = Some((c, d));
            let contig = c as usize;
            if let Some(a) = align_at(strand, contig, &seqs[contig], d)
                && a.score >= min_score
                && a.start < lens[contig]
            {
                found.push(a);
            }
        }
    }
    let Some(best) = found.iter().map(|a| a.score).max() else {
        return found;
    };
    found.retain(|a| a.score == best);
    found
}

/// Aligns all of `read` (bar clipped ends) to `seq`, contig `contig`'s indexed bases, near
/// diagonal `diagonal` (the contig position of the read's first base), in a band of [`BAND`]
/// bases.
fn align_at(read: &[u8], contig: usize, seq: &[u8], diagonal: i64) -> Option<Alignment> {
    let n = read.len();
    let from = usize::try_from((diagonal - BAND as i64).max(0)).ok()?;
    let to = usize::try_from((diagonal + n as i64 + BAND as i64).min(seq.len() as i64)).ok()?;
    if from >= to {
        return None;
    }
    let window = &seq[from..to];
    let m = window.len();
    // Column j of row i is window position `i + offset + j - BAND` for read base i.
    let offset = diagonal - from as i64;
    let width = 2 * BAND + 1;
    const NONE: i32 = i32::MIN / 2;
    // Scores and moves (0 start, 1 diagonal, 2 up = read base inserted, 3 left = deletion).
    let mut score = vec![NONE; (n + 1) * width];
    let mut moves = vec![0u8; (n + 1) * width];
    let at = |i: usize, j: usize| i * width + j;
    let column = |i: usize, j: usize| -> Option<usize> {
        let x = i as i64 + offset + j as i64 - BAND as i64;
        usize::try_from(x).ok().filter(|&x| x <= m)
    };
    for i in 0..=n {
        for j in 0..width {
            let Some(x) = column(i, j) else {
                continue;
            };
            // Start here: free on the contig, a clip on the read unless at its start.
            let mut best = if i == 0 { 0 } else { CLIP };
            let mut mv = 0;
            if i > 0 && x > 0 {
                // Diagonal: same j in row i - 1 is position x - 1.
                let s = score[at(i - 1, j)];
                if s > NONE {
                    let sub = if read[i - 1] == window[x - 1] && read[i - 1] != b'N' {
                        MATCH
                    } else {
                        MISMATCH
                    };
                    if s + sub > best {
                        best = s + sub;
                        mv = 1;
                    }
                }
            }
            if i > 0 && j + 1 < width {
                // Up: row i - 1 at position x is column j + 1.
                let s = score[at(i - 1, j + 1)];
                if s > NONE && s + GAP > best {
                    best = s + GAP;
                    mv = 2;
                }
            }
            if j > 0 && x > 0 {
                let s = score[at(i, j - 1)];
                if s > NONE && s + GAP > best {
                    best = s + GAP;
                    mv = 3;
                }
            }
            score[at(i, j)] = best;
            moves[at(i, j)] = mv;
        }
    }
    // The end: the read's last base, or a clip.
    let mut end = None;
    for i in 0..=n {
        for j in 0..width {
            if column(i, j).is_none() {
                continue;
            }
            let s = score[at(i, j)] + if i == n { 0 } else { CLIP };
            if end.is_none_or(|(b, _, _)| s > b) {
                end = Some((s, i, j));
            }
        }
    }
    let (total, mut i, mut j) = end?;
    let clipped_end = i < n;
    // Trace back: pairs of (read base or None, contig position or None), last first.
    let mut path: Vec<(Option<u8>, Option<usize>)> = Vec::new();
    loop {
        match moves[at(i, j)] {
            1 => {
                let x = column(i, j).expect("on the band");
                path.push((Some(read[i - 1]), Some(x - 1)));
                i -= 1;
            }
            2 => {
                path.push((Some(read[i - 1]), None));
                i -= 1;
                j += 1;
            }
            3 => {
                let x = column(i, j).expect("on the band");
                path.push((None, Some(x - 1)));
                j -= 1;
            }
            _ => break,
        }
    }
    path.reverse();
    // Inserted bases before the first contig position are not placed; those after the last
    // belong to it.
    let first = path.iter().position(|p| p.1.is_some())?;
    let start = path[first].1.expect("found above");
    let mut tokens: Vec<Token> = Vec::new();
    for &(base, pos) in &path[first..] {
        match pos {
            Some(_) => tokens.push(Token {
                base,
                inserted: Vec::new(),
            }),
            None => tokens
                .last_mut()
                .expect("starts on a contig position")
                .inserted
                .push(base.expect("an insertion has a base")),
        }
    }
    left_align(&mut tokens, &window[start..]);
    Some(Alignment {
        score: total,
        contig,
        start: from + start,
        tokens,
        clipped: [i > 0, clipped_end],
    })
}

/// Moves each insertion and deletion in `tokens` (one per base of `contig` from its start) as
/// far left as it goes over matching bases, so that every read puts an indel in a homopolymer
/// or a tandem repeat at the same place and their votes add up.
fn left_align(tokens: &mut [Token], contig: &[u8]) {
    // Insertions after position k: while base k matches and is the insertion's last base, the
    // insertion moves before it.
    for k in 1..tokens.len() {
        let mut k = k;
        while k >= 1
            && !tokens[k].inserted.is_empty()
            && tokens[k].base == Some(contig[k])
            && tokens[k].inserted.last() == Some(&contig[k])
            && tokens[k - 1].inserted.is_empty()
        {
            let mut moved = std::mem::take(&mut tokens[k].inserted);
            moved.pop();
            moved.insert(0, contig[k]);
            tokens[k - 1].inserted = moved;
            k -= 1;
        }
    }
    // Deletions of positions a..b: while base a - 1 matches and equals base b - 1, the deletion
    // moves one left.
    let mut b = 0;
    while b < tokens.len() {
        if tokens[b].base.is_some() || !tokens[b].inserted.is_empty() {
            b += 1;
            continue;
        }
        let mut a = b;
        while b < tokens.len() && tokens[b].base.is_none() && tokens[b].inserted.is_empty() {
            b += 1;
        }
        let end = b;
        while a >= 1
            && tokens[a - 1].base == Some(contig[a - 1])
            && tokens[a - 1].inserted.is_empty()
            && contig[a - 1] == contig[b - 1]
        {
            tokens[a - 1].base = None;
            tokens[b - 1].base = Some(contig[b - 1]);
            a -= 1;
            b -= 1;
        }
        b = end;
    }
}

/// Polishes `contigs` with `reads` (see the module docs) and returns what it did.
pub fn polish_short(contigs: &mut [Contig], reads: &Reads, params: &ShortParams) -> ShortStats {
    let lens: Vec<usize> = contigs.iter().map(|c| c.seq.len()).collect();
    let indexed: Vec<Vec<u8>> = contigs
        .iter()
        .map(|c| {
            if c.circular {
                let wrap = WRAP.min(c.seq.len());
                [c.seq.as_slice(), &c.seq[..wrap]].concat()
            } else {
                c.seq.clone()
            }
        })
        .collect();
    let index = Index::new(&indexed);
    let votes: Vec<Votes> = lens.iter().map(|&len| Votes::new(len)).collect();
    let indices: Vec<_> = reads.indices().collect();
    let aligned: usize = indices
        .par_iter()
        .map(|&r| {
            let alignments = align_read(reads.get(r), &indexed, &lens, &index);
            for a in &alignments {
                let len = lens[a.contig];
                let end = (a.start + a.tokens.len()).min(indexed[a.contig].len());
                for k in a.voting(&indexed[a.contig][a.start..end]) {
                    votes[a.contig].add((a.start + k) % len, a.tokens[k].clone());
                }
            }
            usize::from(!alignments.is_empty())
        })
        .sum();
    let mut stats = ShortStats {
        aligned,
        positions: lens.iter().sum(),
        ..ShortStats::default()
    };
    for (contig, votes) in contigs.iter_mut().zip(votes) {
        let inserted = votes.inserted.into_inner().expect("lock poisoned");
        let mut extra: HashMap<usize, Vec<(Token, u32)>> = HashMap::new();
        for ((pos, token), n) in inserted {
            extra.entry(pos).or_default().push((token, n));
        }
        let mut seq = Vec::with_capacity(contig.seq.len());
        // Where each place of the contig (each base, and the end) is in the new bases: kept
        // only for a contig with places flagged, to move them.
        let follow = !contig.check.flagged.is_empty();
        let mut moved: Vec<u32> = Vec::with_capacity(if follow { contig.seq.len() + 1 } else { 0 });
        let place =
            |seq: &Vec<u8>| u32::try_from(seq.len()).expect("contigs are shorter than 4 Gb");
        for (pos, &have) in contig.seq.iter().enumerate() {
            if follow {
                moved.push(place(&seq));
            }
            let plain: Vec<u32> = votes.plain[pos]
                .iter()
                .map(|a| a.load(Ordering::Relaxed))
                .collect();
            let others = extra.get(&pos).map_or(&[][..], Vec::as_slice);
            let depth: u32 = plain.iter().sum::<u32>() + others.iter().map(|o| o.1).sum::<u32>();
            if depth < params.min_depth {
                seq.push(have);
                continue;
            }
            stats.covered += 1;
            let current = plain.get(code(have)).copied().filter(|_| code(have) < 4);
            // The version most reads show: a plain base, a deletion, or one with insertions.
            let mut best: (u32, Option<u8>, &[u8]) = (0, Some(have), &[]);
            for (slot, &n) in plain.iter().enumerate() {
                if n > best.0 {
                    best = (n, BASES.get(slot).copied(), &[]);
                }
            }
            for (token, n) in others {
                if *n > best.0 {
                    best = (*n, token.base, &token.inserted);
                }
            }
            let (n, base, inserted) = best;
            let changes = (base != Some(have) || !inserted.is_empty())
                && f64::from(n) >= params.valid_share * f64::from(depth)
                && f64::from(current.unwrap_or(0)) <= params.invalid_share * f64::from(depth);
            if changes {
                stats.changed += 1;
                seq.extend(base);
                seq.extend_from_slice(inserted);
            } else {
                seq.push(have);
            }
        }
        if follow {
            moved.push(place(&seq));
            let end = place(&seq);
            contig
                .check
                .lift(|at| moved.get(at as usize).copied().unwrap_or(end));
        }
        contig.seq = seq;
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::check::Check;

    fn random(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                BASES[(x % 4) as usize]
            })
            .collect()
    }

    fn contig(seq: Vec<u8>, circular: bool) -> Contig {
        Contig {
            seq,
            circular,
            hairpins: [false; 2],
            depth: 0.0,
            check: Check::default(),
        }
    }

    /// 150-base reads from both strands of `genome` every `step` bases (round its end when
    /// circular).
    fn reads(genome: &[u8], step: usize, circular: bool) -> Reads {
        let ring = [genome, &genome[..150]].concat();
        let end = if circular {
            genome.len()
        } else {
            genome.len() - 150
        };
        let mut reads = Reads::new();
        let mut rc = Vec::new();
        for start in (0..end).step_by(step) {
            let read = &ring[start..start + 150];
            reads.push(read);
            reverse_complement_into(read, &mut rc);
            reads.push(&rc);
        }
        reads
    }

    #[test]
    fn a_missing_base_in_a_homopolymer_is_put_back() {
        let mut genome = random(20_000, 1);
        genome.splice(10_000..10_000, *b"AAAAAAAA");
        let mut draft = genome.clone();
        draft.remove(10_003);
        let mut contigs = vec![contig(draft, false)];
        let stats = polish_short(&mut contigs, &reads(&genome, 7, false), &ShortParams::new());
        assert_eq!(contigs[0].seq, genome);
        assert_eq!(stats.changed, 1);
    }

    #[test]
    fn flagged_places_move_with_the_bases_put_back() {
        let mut genome = random(20_000, 1);
        genome.splice(10_000..10_000, *b"AAAAAAAA");
        let mut draft = genome.clone();
        draft.remove(10_003);
        let len = u32::try_from(draft.len()).expect("short");
        let bin = |start, end| crate::check::FlaggedBin {
            start,
            end,
            flag: crate::check::Flag::Shallow,
            depth: 1.0,
            stops: 0,
        };
        let mut contigs = vec![Contig {
            check: Check {
                flagged: vec![bin(2_000, 3_000), bin(9_000, 10_003), bin(15_000, len)],
                ..Check::default()
            },
            ..contig(draft, false)
        }];
        polish_short(&mut contigs, &reads(&genome, 7, false), &ShortParams::new());
        assert_eq!(contigs[0].seq, genome);
        let at: Vec<(u32, u32)> = contigs[0]
            .check
            .flagged
            .iter()
            .map(|b| (b.start, b.end))
            .collect();
        // Before the base put back, unmoved; at it and after, one on; the end, the new end.
        let end = u32::try_from(genome.len()).expect("short");
        assert_eq!(at, [(2_000, 3_000), (9_000, 10_004), (15_001, end)]);
    }

    #[test]
    fn two_bases_missing_from_a_homopolymer_and_a_tandem_repeat_are_put_back() {
        let mut genome = random(20_000, 11);
        genome.splice(6_000..6_000, *b"CGGGGGGGA");
        genome.splice(12_000..12_000, *b"TAGAGAGAGAGAGC");
        let mut draft = genome.clone();
        draft.drain(12_005..12_007);
        draft.drain(6_002..6_004);
        let mut contigs = vec![contig(draft, false)];
        let stats = polish_short(&mut contigs, &reads(&genome, 7, false), &ShortParams::new());
        assert_eq!(contigs[0].seq, genome);
        assert_eq!(stats.changed, 2);
    }

    #[test]
    fn two_extra_bases_in_a_homopolymer_are_taken_out() {
        let mut genome = random(20_000, 12);
        genome.splice(6_000..6_000, *b"ATTTTTTC");
        let mut draft = genome.clone();
        draft.splice(6_003..6_003, *b"TT");
        let mut contigs = vec![contig(draft, false)];
        polish_short(&mut contigs, &reads(&genome, 7, false), &ShortParams::new());
        assert_eq!(contigs[0].seq, genome);
    }

    #[test]
    fn a_wrong_base_and_an_extra_base_are_fixed_across_a_circles_end() {
        let genome = random(30_000, 2);
        let mut draft = genome.clone();
        draft[29_990] = if draft[29_990] == b'A' { b'C' } else { b'A' };
        draft.insert(20, b'G');
        let mut contigs = vec![contig(draft, true)];
        polish_short(&mut contigs, &reads(&genome, 5, true), &ShortParams::new());
        assert_eq!(contigs[0].seq, genome);
    }

    #[test]
    fn an_error_in_one_copy_of_a_repeat_longer_than_the_reads_is_left_and_nothing_else_changes() {
        let unit = random(3_000, 3);
        let mut genome = random(10_000, 4);
        genome.extend_from_slice(&unit);
        genome.extend(random(10_000, 5));
        genome.extend_from_slice(&unit);
        genome.extend(random(10_000, 6));
        let mut draft = genome.clone();
        let at = 10_000 + 3_000 + 10_000 + 1_500;
        draft[at] = if draft[at] == b'T' { b'G' } else { b'T' };
        let mut contigs = vec![contig(draft.clone(), false)];
        let stats = polish_short(&mut contigs, &reads(&genome, 6, false), &ShortParams::new());
        assert_eq!(contigs[0].seq, draft);
        assert_eq!(stats.changed, 0);
    }

    #[test]
    fn copies_of_a_repeat_that_truly_differ_are_left_as_they_are() {
        let unit = random(3_000, 7);
        let mut other = unit.clone();
        other[1_500] = if other[1_500] == b'C' { b'A' } else { b'C' };
        let mut genome = random(10_000, 8);
        genome.extend_from_slice(&unit);
        genome.extend(random(10_000, 9));
        genome.extend_from_slice(&other);
        genome.extend(random(10_000, 10));
        let mut contigs = vec![contig(genome.clone(), false)];
        let stats = polish_short(&mut contigs, &reads(&genome, 6, false), &ShortParams::new());
        assert_eq!(contigs[0].seq, genome);
        assert_eq!(stats.changed, 0);
    }
}
