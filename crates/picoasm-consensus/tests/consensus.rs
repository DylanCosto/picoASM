//! Consensus on tiny simulated HiFi-like windows and on edge cases.

use picoasm_consensus::{ConsensusBackend, PoaBackend, PoaParams, Window};

/// A tiny deterministic generator (splitmix64), so the tests need no RNG dependency.
struct SplitMix(u64);

impl SplitMix {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut x = self.0;
        x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        x ^ (x >> 31)
    }

    /// Returns a value in `0..n`.
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }

    /// Returns `true` with probability `per_mille / 1000`.
    fn chance(&mut self, per_mille: u64) -> bool {
        self.below(1000) < per_mille
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[usize::try_from(self.below(4)).expect("below 4 fits in usize")]
    }

    fn other_base(&mut self, base: u8) -> u8 {
        loop {
            let other = self.base();
            if other != base {
                return other;
            }
        }
    }
}

fn random_genome(rng: &mut SplitMix, len: usize) -> Vec<u8> {
    (0..len).map(|_| rng.base()).collect()
}

/// Copies `truth` with HiFi-like errors: substitutions at about `sub_per_mille` per base, and
/// homopolymer runs one base longer or shorter at about `indel_per_mille` per run.
fn hifi_copy(
    rng: &mut SplitMix,
    truth: &[u8],
    sub_per_mille: u64,
    indel_per_mille: u64,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(truth.len() + 8);
    for run in truth.chunk_by(|a, b| a == b) {
        let mut len = run.len();
        if rng.chance(indel_per_mille) {
            if rng.chance(500) {
                len += 1;
            } else {
                len -= 1;
            }
        }
        for _ in 0..len {
            let base = if rng.chance(sub_per_mille) {
                rng.other_base(run[0])
            } else {
                run[0]
            };
            out.push(base);
        }
    }
    out
}

/// A random truth of 400-600 bases and `depth` HiFi-like copies of it.
fn simulated_window(seed: u64, depth: usize) -> (Vec<u8>, Window) {
    let mut rng = SplitMix(seed);
    let len = 400 + usize::try_from(rng.below(201)).expect("below 201 fits in usize");
    let truth = random_genome(&mut rng, len);
    let window = (0..depth)
        .map(|_| hifi_copy(&mut rng, &truth, 10, 20))
        .collect();
    (truth, window)
}

fn backend() -> PoaBackend {
    PoaBackend::new(PoaParams::default()).expect("default parameters are valid")
}

#[test]
fn recovers_the_truth_on_fixed_seeds() {
    let backend = backend();
    for seed in [1, 2, 3, 4, 5] {
        let (truth, window) = simulated_window(seed, 10);
        assert!(
            window.segments().iter().any(|segment| *segment != truth),
            "seed {seed}: the copies should carry errors"
        );
        assert_eq!(
            backend.consensus_window(&window),
            truth,
            "seed {seed}: consensus differs from the truth"
        );
    }
}

#[test]
fn recovers_the_truth_in_nearly_all_trials() {
    let backend = backend();
    let trials: Vec<(Vec<u8>, Window)> = (100..124)
        .map(|seed| simulated_window(seed, 8 + usize::try_from(seed % 5).expect("fits")))
        .collect();
    let windows: Vec<Window> = trials.iter().map(|(_, window)| window.clone()).collect();
    let consensus = backend.consensus(&windows);
    assert_eq!(consensus.len(), trials.len());
    let exact = trials
        .iter()
        .zip(&consensus)
        .filter(|((truth, _), got)| truth == *got)
        .count();
    assert!(
        exact + 1 >= trials.len(),
        "only {exact} of {} consensus sequences equal the truth",
        trials.len()
    );
}

#[test]
fn batch_output_matches_single_windows_for_any_thread_count() {
    let backend = backend();
    let windows: Vec<Window> = (200..206)
        .map(|seed| simulated_window(seed, 6).1)
        .chain([Window::default()])
        .collect();
    let single: Vec<Vec<u8>> = windows
        .iter()
        .map(|window| backend.consensus_window(window))
        .collect();
    for threads in [1, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a tiny thread pool starts");
        let batch = pool.install(|| backend.consensus(&windows));
        assert_eq!(batch, single, "{threads} threads");
    }
}

#[test]
fn empty_windows_give_empty_output() {
    let backend = backend();
    assert!(backend.consensus(&[]).is_empty());
    let windows = [Window::default(), Window::new(vec![Vec::new(), Vec::new()])];
    assert_eq!(backend.consensus(&windows), vec![Vec::<u8>::new(); 2]);
}

#[test]
fn single_segment_is_returned() {
    let backend = backend();
    let window = Window::new(vec![Vec::new(), b"ACGTTGCA".to_vec()]);
    assert_eq!(backend.consensus_window(&window), b"ACGTTGCA");
    assert_eq!(backend.consensus_segments(&[b"acgtn"]), b"ACGTN");
}

#[test]
fn identical_segments_give_that_segment() {
    let backend = backend();
    let mut rng = SplitMix(7);
    let truth = random_genome(&mut rng, 300);
    let segments = vec![truth.clone(); 5];
    assert_eq!(backend.consensus_segments(&segments), truth);
    assert_eq!(
        backend.consensus_segments(&[b"ACGT".as_slice(), b"ACGT"]),
        b"ACGT"
    );
}

#[test]
fn truncated_segments_do_not_change_the_consensus() {
    let backend = backend();
    let (truth, window) = simulated_window(11, 8);
    let mut segments = window.into_segments();
    let half = truth.len() / 2;
    segments.push(truth[..half].to_vec());
    segments.push(truth[half + 40..].to_vec());
    segments.push(truth[100..140].to_vec());
    assert_eq!(backend.consensus_segments(&segments), truth);
}

#[test]
fn an_overlong_segment_does_not_change_the_consensus() {
    let backend = backend();
    let (truth, window) = simulated_window(12, 8);
    let mut segments = window.into_segments();
    let mut rng = SplitMix(99);
    let mut long = random_genome(&mut rng, 150);
    long.extend_from_slice(&truth);
    long.extend(random_genome(&mut rng, 150));
    segments.push(long);
    assert_eq!(backend.consensus_segments(&segments), truth);
}

#[test]
fn non_acgt_bytes_are_outvoted() {
    let backend = backend();
    let (truth, window) = simulated_window(13, 8);
    let mut segments = window.into_segments();
    for (i, segment) in segments.iter_mut().enumerate().take(3) {
        let pos = 50 + 97 * i;
        segment[pos] = b'N';
        segment[pos + 1] = b'-';
    }
    let lower: Vec<u8> = segments[3].to_ascii_lowercase();
    segments[3] = lower;
    assert_eq!(backend.consensus_segments(&segments), truth);
}

#[test]
fn unknown_bases_never_match_each_other() {
    let backend = backend();
    // Every segment has N at the same place: the column has no majority base, so the
    // consensus keeps one N rather than dropping the position.
    let segments = [b"ACGTNACGT".as_slice(), b"ACGTNACGT", b"ACGTNACGT"];
    assert_eq!(backend.consensus_segments(&segments), b"ACGTNACGT");
}

#[test]
fn narrow_band_widens_instead_of_failing() {
    let params = PoaParams::default().with_band_width(1);
    let backend = PoaBackend::new(params).expect("valid parameters");
    let (truth, window) = simulated_window(14, 8);
    let mut segments = window.into_segments();
    let mut insert = truth.clone();
    insert.splice(200..200, b"GATTACAGATTACAGATTACA".iter().copied());
    segments.push(insert);
    assert_eq!(backend.consensus_segments(&segments), truth);
}

#[test]
fn canonically_sorted_segments_give_the_same_consensus_in_any_input_order() {
    let backend = backend();
    // Even depth and a two-segment pair, where the tie-break and so the input order matter.
    let (_, window) = simulated_window(15, 4);
    let segments = window.into_segments();
    let canonical = |mut segments: Vec<Vec<u8>>| {
        segments.sort();
        backend.consensus_segments(&segments)
    };
    let want = canonical(segments.clone());
    let mut rotated = segments.clone();
    for _ in 0..segments.len() {
        rotated.rotate_left(1);
        assert_eq!(canonical(rotated.clone()), want);
        let reversed: Vec<Vec<u8>> = rotated.iter().rev().cloned().collect();
        assert_eq!(canonical(reversed), want);
    }
    let pair = vec![segments[0].clone(), segments[1].clone()];
    let swapped = vec![segments[1].clone(), segments[0].clone()];
    assert_eq!(canonical(pair), canonical(swapped));
}
