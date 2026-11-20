//! Depth normalisation: the reads of high-copy replicons thinned before the layout.
//!
//! A small plasmid can be held in hundreds of copies per cell, so its reads come at thousands
//! of times the chromosome's depth. Their minimizers are then too frequent to seed overlaps
//! ([`crate::overlap::OverlapParams::max_occurrences`]), and if they did seed, the overlaps
//! between thousands of reads of one short circle would cost more than the rest of the
//! genome. So each read's depth is estimated from the counts, over all reads, of its
//! minimizers; the commonest count of a distinct minimizer (past the counts of errors) is the
//! genome's depth, the chromosome's; and a read deeper than [`MAX_DEPTH_RATIO`] times
//! that is kept for the layout with a chance that brings its replicon down to that depth,
//! decided by a hash of its index so that every run keeps the same reads. Where the whole genome
//! is deeper than the layout needs, only the longest reads are kept, as Flye's and Raven's
//! read selection does, and then the reads of whatever those leave uncovered: a plasmid of a
//! few kilobases has only short reads ([`MIN_COVERED`]). Polishing still uses every read.

use crate::overlap::{ReadIndex, Reads};
use crate::seq::{Anchor, InvalidSketchParamsError, SketchParams, Sketcher};
use rayon::prelude::*;

/// Reads at up to this many times the genome's depth are all kept (the default
/// [`crate::BactParams::max_depth_ratio`]).
pub(crate) const MAX_DEPTH_RATIO: f64 = 3.0;

/// A read is thinned only if its solid minimizers' (those seen [`MIN_SOLID_COUNT`] times or
/// more) count this fraction of the way up (the median) is deep too, and as much as the median of
/// all its counts says. Counting all minimizers, a quarter of an Oxford Nanopore read's hold an
/// error, so a lower quantile was about 1 for every read; a lower quantile of the solid ones
/// (one fifth) kept the reads of HiFi rRNA operons but too many of a small plasmid's for their
/// minimizers to seed overlaps. A chromosome whose rRNA-operon reads are thinned can lose its
/// layout's circle there; joining the drafts' ends (`join.rs`) closes it again.
const DEPTH_QUANTILE_DIVISOR: usize = 2;

/// The layout keeps the longest reads up to this depth, counted as the genome's depth is (the
/// default [`crate::BactParams::max_layout_depth`]): about 70x of R10.4.1 reads.
pub(crate) const MAX_LAYOUT_DEPTH: f64 = 60.0;

/// The layout keeps the longest reads up to at most this many times the genome's length in
/// bases, whatever [`MAX_LAYOUT_DEPTH`] allows. The genome's depth counts the reads that hold a
/// minimizer without an error, a share of all reads that falls with their errors (about
/// 0.98^19 of R10.4.1 reads' 19-mers, 0.93^19 of reads 7% off): on an *A. mictus* run of 288x
/// of such reads it came out at 92, the layout kept 194x, and overlapping them took 634 s. The
/// genome's length is estimated from its solid minimizers, and comes out at about 70% of it
/// (4.30 Mb for a 5.7 Mb *K. pneumoniae*, 3.30 for a 4.8 Mb *S. flexneri*, 1.49 for the 2.13 Mb
/// *A. mictus*), so this is about 105x: above what the depth already keeps of R10.4.1 reads
/// (the *S. flexneri* at 200x keeps 460 Mb either way).
pub(crate) const MAX_LAYOUT_BASE_DEPTH: f64 = 150.0;

/// After the longest reads are chosen, a read left out is taken back if its solid minimizers'
/// median count among the chosen reads is below this: its replicon (a plasmid shorter than the
/// length cut) or stretch has no long reads.
const MIN_COVERED: u32 = 5;

/// Minimizers seen fewer times than this are taken to hold an error.
const MIN_SOLID_COUNT: usize = 3;

/// What [`normalise`] chose.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Selection {
    /// The reads kept, in index order.
    pub(crate) kept: Vec<ReadIndex>,
    /// The genome's depth, as estimated.
    pub(crate) depth: f64,
}

/// Picks the reads to lay out (see the module docs), estimating depths from minimizers of
/// `kmer_len` bases at `density`; where the genome is deeper than `max_depth` (in the same
/// units), keeps only the longest reads, enough for that depth.
///
/// # Errors
///
/// Returns an error if `kmer_len` or `density` cannot be used.
pub(crate) fn normalise(
    reads: &Reads,
    kmer_len: u8,
    density: f64,
    max_ratio: f64,
    max_depth: f64,
) -> Result<Selection, InvalidSketchParamsError> {
    let sketch = SketchParams::new(kmer_len, density, 1)?;
    let indices: Vec<ReadIndex> = reads.indices().collect();
    let sketches: Vec<Vec<u64>> = indices
        .par_iter()
        .map_init(
            || (Sketcher::new(sketch), Vec::<Anchor>::new()),
            |(sketcher, anchors), &read| {
                sketcher
                    .sketch(reads.get(read), anchors)
                    .expect("read lengths were checked when added");
                anchors.iter().map(|a| a.hash).collect()
            },
        )
        .collect();
    let mut all: Vec<u64> = sketches.iter().flatten().copied().collect();
    all.par_sort_unstable();
    let count = |hash: u64| -> u32 {
        let from = all.partition_point(|&h| h < hash);
        let to = all.partition_point(|&h| h <= hash);
        u32::try_from(to - from).unwrap_or(u32::MAX)
    };
    // Each read's (median of its solid minimizers' counts, median of all its counts).
    let estimates: Vec<(u32, u32)> = sketches
        .par_iter()
        .map(|hashes| {
            let mut counts: Vec<u32> = hashes.iter().map(|&h| count(h)).collect();
            if counts.is_empty() {
                return (0, 0);
            }
            let n = counts.len();
            let median = *counts.select_nth_unstable(n / 2).1;
            // Of the solid minimizers only: those holding an error are about a quarter of an
            // Oxford Nanopore read's, and would pull every read's quantile down.
            counts.retain(|&c| c as usize >= MIN_SOLID_COUNT);
            let solid = counts.len();
            let low = if solid == 0 {
                0
            } else {
                *counts.select_nth_unstable(solid / DEPTH_QUANTILE_DIVISOR).1
            };
            (low, median)
        })
        .collect();
    // The genome's depth: the commonest count of a distinct minimizer, past the counts of
    // errors (1 and 2). The chromosome holds most distinct minimizers; a high-copy plasmid
    // holds few, however deep.
    let mut histogram: Vec<u64> = Vec::new();
    for run in all.chunk_by(|a, b| a == b) {
        let c = run.len();
        if c >= histogram.len() {
            histogram.resize(c + 1, 0);
        }
        histogram[c] += 1;
    }
    let depth = histogram_peak(&histogram);
    // The genome's length: its distinct solid minimizers over their density.
    let solid_distinct: u64 = histogram
        .iter()
        .skip(solid_trough(&histogram).unwrap_or(usize::MAX))
        .sum();
    #[expect(clippy::cast_precision_loss, reason = "an estimate")]
    let genome_len = solid_distinct as f64 / density;
    log::info!(
        "bacteria: about {:.2} Mb of genome, by its solid minimizers",
        genome_len / 1e6
    );
    log::debug!(
        "bacteria: minimizer count histogram (count: distinct minimizers): {}",
        histogram
            .iter()
            .enumerate()
            .skip(1)
            .take(80)
            .map(|(c, n)| format!("{c}:{n}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
    let depth = f64::from(depth.max(1));
    let mut kept: Vec<ReadIndex> = indices
        .into_iter()
        .zip(&estimates)
        .filter(|&(read, &(low, median))| {
            // High-copy all along: thinned to the depth its median count gives.
            let ratio = max_ratio * depth / f64::from(median.max(1));
            f64::from(low) <= max_ratio * depth || ratio >= 1.0 || unit_hash(read.0) < ratio
        })
        .map(|(read, _)| read)
        .collect();
    // Deep data: the longest reads, up to the depth the layout needs.
    let total: u64 = kept.iter().map(|&r| u64::from(reads.read_len(r))).sum();
    #[expect(clippy::cast_precision_loss, reason = "a base count")]
    let base_cap = (genome_len > 0.0)
        .then_some(MAX_LAYOUT_BASE_DEPTH * genome_len)
        .filter(|&cap| total as f64 > cap);
    if depth > max_depth || base_cap.is_some() {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss,
            reason = "a share of a base count"
        )]
        let want = (total as f64 * (max_depth / depth).min(1.0))
            .min(base_cap.unwrap_or(f64::INFINITY)) as u64;
        kept.sort_by_key(|&r| (std::cmp::Reverse(reads.read_len(r)), r));
        let mut seen = 0;
        let cut = kept
            .iter()
            .position(|&r| {
                let enough = seen >= want;
                seen += u64::from(reads.read_len(r));
                enough
            })
            .unwrap_or(kept.len());
        let left_out = kept.split_off(cut);
        // The minimizers of the longest reads, and how often each occurs among them.
        let mut chosen: Vec<u64> = kept
            .iter()
            .flat_map(|r| sketches[r.index()].iter().copied())
            .collect();
        chosen.par_sort_unstable();
        let among_chosen = |hash: u64| -> u32 {
            let from = chosen.partition_point(|&h| h < hash);
            let to = chosen.partition_point(|&h| h <= hash);
            u32::try_from(to - from).unwrap_or(u32::MAX)
        };
        let taken_back: Vec<ReadIndex> = left_out
            .into_par_iter()
            .filter(|r| {
                let mut counts: Vec<u32> = sketches[r.index()]
                    .iter()
                    .filter(|&&h| count(h) as usize >= MIN_SOLID_COUNT)
                    .map(|&h| among_chosen(h))
                    .collect();
                let n = counts.len();
                // Thinned to the layout's depth by the depth all reads give it, as above.
                let ratio = max_depth / f64::from(estimates[r.index()].1.max(1));
                n > 0
                    && *counts.select_nth_unstable(n / 2).1 < MIN_COVERED
                    && (ratio >= 1.0 || unit_hash(r.0) < ratio)
            })
            .collect();
        if !taken_back.is_empty() {
            log::info!(
                "bacteria: {} shorter reads of stretches the longest leave uncovered kept",
                taken_back.len()
            );
        }
        kept.extend(taken_back);
        kept.sort_unstable();
    }
    Ok(Selection { kept, depth })
}

/// Returns the count at the histogram's solid peak: past the first trough after the counts of
/// errors (where the number of distinct minimizers stops falling), the count held by the most;
/// 1 if the histogram never rises again.
fn histogram_peak(histogram: &[u64]) -> u32 {
    solid_trough(histogram)
        .and_then(|t| {
            histogram
                .iter()
                .enumerate()
                .skip(t)
                .max_by_key(|&(c, &n)| (n, std::cmp::Reverse(c)))
        })
        .map_or(1, |(c, _)| u32::try_from(c).unwrap_or(u32::MAX))
}

/// Returns the first count past the counts of errors where the number of distinct minimizers
/// stops falling, if it does.
fn solid_trough(histogram: &[u64]) -> Option<usize> {
    (MIN_SOLID_COUNT..histogram.len().saturating_sub(1)).find(|&c| histogram[c] < histogram[c + 1])
}

/// Returns a number in `[0, 1)` from a read index (SplitMix64's finaliser).
fn unit_hash(x: u32) -> f64 {
    let mut z = u64::from(x).wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    #[expect(clippy::cast_precision_loss, reason = "53 bits of a hash are plenty")]
    let unit = (z >> 11) as f64 / (1_u64 << 53) as f64;
    unit
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noisy_deep_data_is_cut_by_the_genomes_length() {
        // A 60 kb "chromosome" at about 270x in reads with a substitution at one base in ten:
        // few of their minimizers are free of errors, so the depth counted in minimizers is low
        // (under the layout's 200 here), and the longest reads are kept up to the base cap.
        let chromosome = genome(60_000, 3);
        let mut state = 7_u64;
        let mut noisy = |seq: &[u8]| -> Vec<u8> {
            seq.iter()
                .map(|&b| {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    if state >> 60 == 0 {
                        if b == b'A' { b'C' } else { b'A' }
                    } else {
                        b
                    }
                })
                .collect()
        };
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for i in 0..2_000 {
            let start = (i * 271) % 56_000;
            seqs.push(noisy(&chromosome[start..start + 4_000]));
        }
        for i in 0..4_000 {
            let start = (i * 97) % 58_000;
            seqs.push(noisy(&chromosome[start..start + 2_000]));
        }
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let s = normalise(&reads, 15, 0.1, MAX_DEPTH_RATIO, 200.0).expect("valid sketch");
        assert!(s.depth < 200.0, "depth {}", s.depth);
        let kept_bases: u64 = s.kept.iter().map(|&r| u64::from(reads.read_len(r))).sum();
        let total: u64 = reads.indices().map(|r| u64::from(reads.read_len(r))).sum();
        assert!(
            kept_bases * 3 < total * 2,
            "{kept_bases} of {total} bases kept"
        );
        // The longest first (short reads come back only where the long ones leave a stretch
        // thin): a far greater share of the 4 kb reads than of the 2 kb ones.
        let long = s.kept.iter().filter(|r| r.0 < 2_000).count();
        let short = s.kept.len() - long;
        assert!(
            long * 4_000 > 4 * short * 2_000,
            "{long} long and {short} short reads kept"
        );
    }

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
    fn the_peak_is_found_past_the_errors() {
        let h = [0, 900, 300, 120, 60, 40, 50, 80, 95, 70, 30, 5];
        assert_eq!(histogram_peak(&h), 8);
        assert_eq!(histogram_peak(&[0, 900, 300, 120, 60]), 1);
    }

    #[test]
    fn a_high_copy_circle_is_thinned_and_the_rest_kept() {
        // A 60 kb "chromosome" at 10x in 3 kb reads, and a 2 kb "plasmid" at 500x.
        let chromosome = genome(60_000, 1);
        let plasmid = genome(2_000, 2);
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for i in 0..200 {
            let start = (i * 2_999) % 57_000;
            seqs.push(chromosome[start..start + 3_000].to_vec());
        }
        for _ in 0..500 {
            seqs.push(plasmid.clone());
        }
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        let s = normalise(&reads, 15, 0.1, MAX_DEPTH_RATIO, 1e9).expect("valid sketch");
        let kept_chromosome = s.kept.iter().filter(|r| r.0 < 200).count();
        let kept_plasmid = s.kept.len() - kept_chromosome;
        assert_eq!(kept_chromosome, 200, "depth {}", s.depth);
        assert!(
            (10..100).contains(&kept_plasmid),
            "{kept_plasmid} plasmid reads kept"
        );
    }

    #[test]
    fn deep_data_keeps_the_longest_reads_and_a_short_plasmids() {
        // A 60 kb "chromosome" at 40x, half in 6 kb reads and half in 2 kb reads, and a 2 kb
        // "plasmid" at the same depth in 1.5 kb reads.
        let chromosome = genome(60_000, 1);
        let plasmid = genome(2_000, 2);
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for i in 0..200 {
            let start = (i * 271) % 54_000;
            seqs.push(chromosome[start..start + 6_000].to_vec());
        }
        for i in 0..600 {
            let start = (i * 97) % 58_000;
            seqs.push(chromosome[start..start + 2_000].to_vec());
        }
        for i in 0..53 {
            let start = (i * 37) % 500;
            seqs.push(plasmid[start..start + 1_500].to_vec());
        }
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        // The depth comes to 42 (in minimizer counts): room for half the bases, the long reads.
        let s = normalise(&reads, 15, 0.1, MAX_DEPTH_RATIO, 21.0).expect("valid sketch");
        let long = s.kept.iter().filter(|r| r.0 < 200).count();
        let short = s.kept.iter().filter(|r| (200..800).contains(&r.0)).count();
        let plasmid_reads = s.kept.iter().filter(|r| r.0 >= 800).count();
        assert_eq!(long, 200, "depth {}", s.depth);
        assert!(short < 100, "{short} short chromosome reads kept");
        // All the plasmid's reads are taken back, thinned to the layout's depth (about half).
        assert!(
            (15..45).contains(&plasmid_reads),
            "{plasmid_reads} plasmid reads kept"
        );
    }
}
