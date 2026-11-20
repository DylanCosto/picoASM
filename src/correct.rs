//! Noisy reads corrected from each other before the layout.
//!
//! Reads that differ from the genome at 8-16% (R9 reads called with a fast model) share too few
//! error-free seeds for the layout's overlaps: most reads find none, and the contigs fall
//! apart. Here the reads are first overlapped with shorter, denser seeds and a looser bound on
//! divergence ([`loose`]); then each read is cut into windows (as a contig is for polishing)
//! and each window becomes the partial-order consensus of the read's own bases and those its
//! overlapping reads align there (up to [`MAX_SEGMENTS`]), so the read keeps its own length and
//! place but loses most of its errors. The layout then overlaps the corrected reads as usual.

use crate::consensus::{ConsensusBackend, PoaBackend, Window};
use crate::overlap::{
    Edit, EditOp, InvalidOverlapParamsError, Overlap, OverlapParams, ReadIndex, Reads,
    find_overlaps,
};
use crate::seq::{Strand, reverse_complement_into};
use rayon::prelude::*;

use crate::polish::{cut_points, window_boundaries};

/// Reads are corrected in windows of this many bases.
const WINDOW_LEN: usize = 500;

/// A window is the consensus of at most this many segments (the read's own first).
const MAX_SEGMENTS: usize = 20;

/// A window is corrected only if at least this many segments (the read's own among them) cover
/// it; otherwise it keeps the read's bases.
const MIN_SEGMENTS: usize = 3;

/// Reads are corrected this many at a time, to bound the segments held at once.
const BATCH_READS: usize = 20_000;

/// The aligner's band for the correction overlaps ([`OverlapParams::max_lag`]). On
/// *A. baumannii* R9 reads (82% identity) a band of 200 cut the correction's CPU time 2.5-fold
/// and left the assembly as accurate (20,247 errors against 20,242 with no band); a band of
/// 150 saved a little more (2.9-fold) but made 22,456 errors, as indels between reads longer
/// than about half the band are aligned as mismatches.
const MAX_LAG: u32 = 200;

/// The settings the reads are overlapped with for correction: `params` with 15-base seeds, twice
/// the density, and overlaps differing at up to a quarter of their columns (two reads that each
/// differ from the genome at 12% differ from each other at about 22%).
#[must_use]
pub(crate) fn loose(params: &OverlapParams) -> OverlapParams {
    OverlapParams {
        kmer_len: 15,
        density: (params.density * 2.0).min(1.0),
        max_divergence: 0.25,
        max_lag: Some(MAX_LAG),
        ..*params
    }
}

/// Returns `reads` corrected from their overlaps under `params` (see the module docs), in the
/// same order.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
pub(crate) fn correct_reads(
    reads: &Reads,
    params: &OverlapParams,
    backend: &PoaBackend,
) -> Result<Reads, InvalidOverlapParamsError> {
    let overlaps = find_overlaps(reads, params)?;
    // Each read's overlaps, and whether it is their target (else the overlap is turned round).
    let mut of_read: Vec<Vec<(usize, bool)>> = vec![Vec::new(); reads.len()];
    for (i, o) in overlaps.iter().enumerate() {
        of_read[o.target.index()].push((i, true));
        of_read[o.query.index()].push((i, false));
    }
    log::info!(
        "bacteria: {} overlaps between the reads to correct them from; {} reads have none",
        overlaps.len(),
        of_read.iter().filter(|o| o.is_empty()).count()
    );
    let mut corrected = Reads::with_capacity(reads.len(), reads.total_bases());
    let indices: Vec<usize> = (0..reads.len()).collect();
    for batch in indices.chunks(BATCH_READS) {
        let gathered: Vec<ReadWindows> = batch
            .par_iter()
            .map(|&r| segments_of(reads, r, &of_read[r], &overlaps))
            .collect();
        let mut windows: Vec<Window> = Vec::new();
        for read in &gathered {
            for held in read.segments.iter().filter(|h| h.len() >= MIN_SEGMENTS) {
                windows.push(Window::new(held.clone()));
            }
        }
        let mut consensus = backend.consensus(&windows).into_iter();
        drop(windows);
        for (&r, windows) in batch.iter().zip(&gathered) {
            let read = reads.get(ReadIndex(index(r)));
            let mut seq = Vec::with_capacity(read.len());
            for (w, held) in windows.segments.iter().enumerate() {
                let own = &read[pos(windows.bounds[w])..pos(windows.bounds[w + 1])];
                if held.len() >= MIN_SEGMENTS {
                    match consensus.next() {
                        Some(c) if !c.is_empty() => seq.extend_from_slice(&c),
                        _ => seq.extend_from_slice(own),
                    }
                } else {
                    seq.extend_from_slice(own);
                }
            }
            corrected.push(&seq);
        }
    }
    Ok(corrected)
}

/// A read's windows and the segments over each.
struct ReadWindows {
    /// Where its windows start, and its length.
    bounds: Vec<u32>,
    /// Over each window, the read's own bases, then those other reads align there.
    segments: Vec<Vec<Vec<u8>>>,
}

/// The windows of read `r` and, over each, its own bases then those of the reads its
/// `overlaps` (indices into `all`, and whether `r` is their target) align there.
fn segments_of(
    reads: &Reads,
    r: usize,
    overlaps: &[(usize, bool)],
    all: &[Overlap],
) -> ReadWindows {
    let read = reads.get(ReadIndex(index(r)));
    let bounds = window_boundaries(read, WINDOW_LEN);
    let mut segments: Vec<Vec<Vec<u8>>> = bounds
        .windows(2)
        .map(|w| vec![read[pos(w[0])..pos(w[1])].to_vec()])
        .collect();
    let mut oriented = Vec::new();
    for &(i, as_target) in overlaps {
        let m = if as_target {
            all[i].clone()
        } else {
            turned(&all[i])
        };
        let other = reads.get(m.query);
        let oriented_start = match m.strand {
            Strand::Forward => {
                oriented.clear();
                oriented.extend_from_slice(other);
                m.query_start
            }
            Strand::Reverse => {
                reverse_complement_into(other, &mut oriented);
                reads.read_len(m.query) - m.query_end
            }
        };
        let cuts = cut_points(&m, oriented_start, &bounds);
        for pair in cuts.windows(2) {
            let ((window, from), (_, to)) = (pair[0], pair[1]);
            let held = &mut segments[window];
            if held.len() < MAX_SEGMENTS {
                held.push(oriented[pos(from)..pos(to)].to_vec());
            }
        }
    }
    ReadWindows { bounds, segments }
}

/// `o` turned round: its query as the target and its target as the query. Read backwards
/// along each other, the alignment also runs the other way.
fn turned(o: &Overlap) -> Overlap {
    let swap = |e: &Edit| match e.op() {
        EditOp::Insertion => Edit::new(EditOp::Deletion, e.len()),
        EditOp::Deletion => Edit::new(EditOp::Insertion, e.len()),
        _ => *e,
    };
    let edits = match o.strand {
        Strand::Forward => o.edits.iter().map(swap).collect(),
        Strand::Reverse => o.edits.iter().rev().map(swap).collect(),
    };
    Overlap {
        query: o.target,
        target: o.query,
        strand: o.strand,
        query_start: o.target_start,
        query_end: o.target_end,
        target_start: o.query_start,
        target_end: o.query_end,
        edits,
    }
}

fn index(r: usize) -> u32 {
    u32::try_from(r).expect("fewer than 4 billion reads")
}

fn pos(x: u32) -> usize {
    usize::try_from(x).expect("positions fit in usize")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[usize::try_from(x >> 62).expect("two bits")]
            })
            .collect()
    }

    /// `seq` with every `every`-th base changed (a stand-in for a noisy read).
    fn noisy(seq: &[u8], every: usize, offset: usize) -> Vec<u8> {
        seq.iter()
            .enumerate()
            .map(|(i, &b)| {
                if (i + offset).is_multiple_of(every) {
                    match b {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        _ => b'A',
                    }
                } else {
                    b
                }
            })
            .collect()
    }

    #[test]
    fn a_turned_overlap_maps_the_other_read() {
        let o = Overlap {
            query: ReadIndex(0),
            target: ReadIndex(1),
            strand: Strand::Forward,
            query_start: 10,
            query_end: 20,
            target_start: 30,
            target_end: 41,
            edits: vec![
                Edit::new(EditOp::Match, 5),
                Edit::new(EditOp::Deletion, 1),
                Edit::new(EditOp::Match, 5),
            ],
        };
        let t = turned(&o);
        assert_eq!((t.query.index(), t.target.index()), (1, 0));
        assert_eq!(
            (t.query_start, t.query_end, t.target_start, t.target_end),
            (30, 41, 10, 20)
        );
        assert_eq!(t.edits[1].op(), EditOp::Insertion);
        assert_eq!(t.columns(), o.columns());
    }

    #[test]
    fn reads_with_errors_in_different_places_are_corrected() {
        // Ten reads of one 6 kb stretch, each with every 25th base wrong at its own offset,
        // half of them on the other strand.
        let truth = genome(6_000, 3);
        let seqs: Vec<Vec<u8>> = (0..10)
            .map(|k| {
                let read = noisy(&truth, 25, k * 3);
                if k % 2 == 0 {
                    read
                } else {
                    crate::seq::reverse_complement(&read)
                }
            })
            .collect();
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let params = loose(&crate::BactParams::new(crate::ReadKind::Ont).overlap);
        let backend = PoaBackend::new(crate::consensus::PoaParams::default()).expect("defaults");
        let corrected = correct_reads(&reads, &params, &backend).expect("valid");
        let wrong = |a: &[u8]| a.iter().zip(&truth).filter(|(x, y)| x != y).count();
        let first = corrected.get(ReadIndex(0));
        assert_eq!(first.len(), truth.len());
        assert!(
            wrong(first) * 10 < wrong(&seqs[0]),
            "{} of {}",
            wrong(first),
            wrong(&seqs[0])
        );
    }
}
