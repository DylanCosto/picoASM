//! The pipeline: overlaps, the string graph, and polishing.

use std::time::Instant;

use picoasm_consensus::{InvalidPoaParamsError, PoaBackend};
use picoasm_overlap::{
    CorrectionParams, Extension, Fate, InvalidOverlapParamsError, Layout, Overlap, OverlapParams,
    ReadIndex, Reads, Usage, correct_in_rounds, find_overlaps_counted, lay_out, map_reads, usage,
};
use picoasm_seq::InvalidSketchParamsError;
use thiserror::Error;

use crate::check::{self, Check};
use crate::hang;
use crate::homopolymer::{self, Evidence};
use crate::motif::Modified;
use crate::params::{BactParams, ReadKind};
use crate::polish::{Draft, polish_round};
use crate::{
    bridge, circle, context, correct, hairpin, join, locus, repeat, select, small, trim, vote,
};

/// A contig of the assembly.
#[derive(Clone, Debug, PartialEq)]
pub struct Contig {
    /// Its bases.
    pub seq: Vec<u8>,
    /// Whether its last base runs on into its first (a whole replicon).
    pub circular: bool,
    /// Whether its start and its end are hairpins: a linear replicon's ends, which reads run
    /// round (both: a whole linear replicon).
    pub hairpins: [bool; 2],
    /// The mean depth of the reads mapped onto it in the last polishing round.
    pub depth: f64,
    /// What the reads show of it: shallow stretches, places many reads stop.
    pub check: Check,
}

impl Contig {
    /// Whether it is a whole replicon: a circle, or a line with hairpins at both ends.
    #[must_use]
    pub fn is_whole(&self) -> bool {
        self.circular || self.hairpins == [true, true]
    }
}

/// What [`assemble`] built.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Assembly {
    /// The contigs, longest first.
    pub contigs: Vec<Contig>,
    /// Contigs left out as too shallow beside the chromosome (see [`MIN_DEPTH_SHARE`] and
    /// [`MIN_LINEAR_DEPTH_SHARE`]): another sample's reads, or a chimera's; and shallow circles
    /// of a longer contig's bases ([`MAX_EXCISED_DEPTH_SHARE`]).
    pub shallow: Vec<Contig>,
    /// The overlaps the layout used, and how (for development).
    pub overlaps: Vec<(Overlap, Usage)>,
    /// The draft contigs' reads, before polishing, and what became of each read (for
    /// development). Reads are numbered as laid out: see [`Assembly::laid_reads`].
    pub layout: Layout,
    /// For each read laid out, its index among the input reads and its length once trimmed
    /// (for development: the overlaps and the layout number reads this way).
    pub laid_reads: Vec<(ReadIndex, u32)>,
    /// The reads' lengths at the polished contigs' homopolymer runs, before they are called
    /// again (for development).
    pub homopolymers: Vec<Evidence>,
}

/// [`assemble`]'s settings cannot be used.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum AssembleError {
    /// The overlap or mapping settings.
    #[error("invalid overlap settings")]
    Overlap(#[from] InvalidOverlapParamsError),
    /// The consensus settings.
    #[error("invalid consensus settings")]
    Poa(#[from] InvalidPoaParamsError),
    /// The minimizer settings.
    #[error("invalid minimizer settings")]
    Sketch(#[from] InvalidSketchParamsError),
    /// The polishing windows are shorter than [`MIN_WINDOW_LEN`].
    #[error("polishing windows of {0} bases are too short")]
    Window(usize),
    /// The polishing depths: none, or the most below the least.
    #[error("polishing depths from {0} to {1} cannot be used")]
    Depth(usize, usize),
}

/// Polishing windows are at least this long.
pub const MIN_WINDOW_LEN: usize = 32;

/// Assembles `reads` of one bacterial isolate.
///
/// Overlaps the reads ([`picoasm_overlap::find_overlaps`]; HiFi reads are first corrected in
/// [`BactParams::correction_rounds`] rounds), lays them out in a string graph
/// ([`picoasm_overlap::lay_out`]), whose unitigs are the draft contigs (a path that closes on
/// itself is a circular one, cut where its last read runs into its first), and polishes them
/// [`BactParams::polish_rounds`] times (at least twice) from all the reads. After the first
/// round, linear drafts whose ends overlap are joined, and closed into circles (`join.rs`),
/// those whose ends lie in another draft are closed by its bases (`bridge.rs`), and ends that
/// reads run between are joined through them (`join.rs`). Between
/// rounds each circular contig is turned half way round, so the place it was cut is polished
/// as any other. Then loci where most reads have
/// another version get it ([`BactParams::repair_loci`]), and the reads call homopolymer
/// lengths and vote on the other bases, motifs one strand misreads learned on the way
/// ([`BactParams::vote_margin`]).
///
/// The result is the same for any number of threads; work runs on the caller's rayon pool.
///
/// # Errors
///
/// Returns an error if `params` cannot be used.
pub fn assemble(reads: &Reads, params: &BactParams) -> Result<Assembly, AssembleError> {
    if params.window_len < MIN_WINDOW_LEN {
        return Err(AssembleError::Window(params.window_len));
    }
    if params.min_depth == 0 || params.max_depth < params.min_depth {
        return Err(AssembleError::Depth(params.min_depth, params.max_depth));
    }
    let backend = PoaBackend::new(params.poa)?;
    let started = Instant::now();
    let selection = select::normalise(
        reads,
        params.overlap.kmer_len,
        params.overlap.density,
        params.max_depth_ratio,
        params.max_layout_depth,
    )?;
    let chosen: Reads = selection.kept.iter().map(|&r| reads.get(r)).collect();
    log::info!(
        "bacteria: depth {:.0}x; {} of {} reads ({} of {} bases) laid out, high-copy replicons \
         thinned",
        selection.depth,
        chosen.len(),
        reads.len(),
        chosen.total_bases(),
        reads.total_bases()
    );
    // Whether ONT reads were found noisy and corrected (see `BactParams::noisy_divergence`).
    let mut found_noisy = false;
    // The reads laid out (trimmed or corrected): the chosen ones, moved in, are not needed after.
    let (laid, overlaps, used): (Reads, Vec<Overlap>, Vec<Usage>) = match (params.kind, chosen) {
        (ReadKind::Ont | ReadKind::OntNoisy, chosen) => {
            // Reads laid out raw: their overlaps, unless they turn out noisy.
            let mut rounds = params.correction_rounds;
            let mut raw = None;
            if rounds == 0 {
                let found = find_overlaps_counted(&chosen, &params.overlap)?;
                if params
                    .noisy_divergence
                    .is_some_and(|max| looks_noisy(&chosen, &found, max))
                {
                    rounds = BactParams::new(ReadKind::OntNoisy).correction_rounds;
                    found_noisy = true;
                    log::warn!(
                        "bacteria: the reads are noisy (R10.4.1 reads called with sup or hac \
                         overlap at 2-5%): they are corrected from each other first, as with \
                         `--reads ont-noisy` (`--reads ont-clean` lays them out raw)"
                    );
                } else {
                    raw = Some(found);
                }
            }
            // Noisy reads are corrected from each other first (see `correct.rs`).
            let mut corrected: Option<Reads> = None;
            for round in 0..rounds {
                let from = corrected.as_ref().unwrap_or(&chosen);
                let next =
                    correct::correct_reads(from, &correct::loose(&params.overlap), &backend)?;
                log::info!(
                    "bacteria: reads corrected, round {}, {:.1} s",
                    round + 1,
                    started.elapsed().as_secs_f64()
                );
                corrected = Some(next);
            }
            let chosen = corrected.as_ref().unwrap_or(&chosen);
            // Only the overlaps' places and how far their reads differ are used below (the
            // edges, the trim and the cut); their alignments are not kept.
            let found = match raw {
                Some(found) => found,
                None => find_overlaps_counted(chosen, &params.overlap)?,
            };
            let edge_divergence = if corrected.is_some() {
                let median = median_divergence(&found);
                let cut = (CORRECTED_EDGE_FACTOR * median)
                    .clamp(MIN_CORRECTED_EDGE_DIVERGENCE, params.edge_divergence);
                log::info!(
                    "bacteria: overlaps between the corrected reads differ at a median {:.2}% \
                     of their columns; edges differ at {:.2}% at most",
                    100.0 * median,
                    100.0 * cut
                );
                cut
            } else {
                params.edge_divergence
            };
            let trimmed = trim::trim(
                chosen,
                &found,
                params.overlap.min_overlap,
                trim::min_cover(selection.depth),
            );
            // Whether an overlap is an edge is decided on its alignment, before it is cut to
            // the reads' trusted stretches and carried over its overhangs.
            let (overlaps, used): (Vec<Overlap>, Vec<Usage>) = found
                .iter()
                .zip(edges(&found, edge_divergence))
                .filter_map(|(o, u)| {
                    let cut = trimmed.cut(o, params.overlap.min_overlap)?;
                    hang::snap(&trimmed.reads, &cut).map(|o| (o, u))
                })
                .unzip();
            log::info!(
                "bacteria: {} of {} overlaps between the trimmed reads ({} of {} bases kept)",
                overlaps.len(),
                found.len(),
                trimmed.reads.total_bases(),
                chosen.total_bases()
            );
            (trimmed.reads, overlaps, used)
        }
        (ReadKind::Hifi, chosen) => {
            let corrected = correct_in_rounds(
                chosen,
                &params.overlap,
                &CorrectionParams::default(),
                params.correction_rounds,
            )?;
            let used = usage(
                &corrected.overlaps,
                &corrected.consistent,
                picoasm_overlap::MAX_RESCUE_DIVERGENCE,
            );
            (corrected.reads, corrected.overlaps, used)
        }
    };
    log::info!(
        "bacteria: {} overlaps between {} reads, {:.1} s",
        overlaps.len(),
        laid.len(),
        started.elapsed().as_secs_f64()
    );
    let mut used = used;
    let mut layout = lay_out(&laid, &overlaps, &used, &params.layout);
    log_layout(&layout);
    // Contigs well below the chromosome's depth that are not whole circles are mostly a few
    // reads beside a deeper path (a chimera's, or a rare variant's); they stop the deeper
    // path from closing. Their reads are left out and the reads laid out again.
    let weak = weak_reads(&layout, &laid, &overlaps, &used, &params.map)?;
    if !weak.is_empty() {
        for (o, u) in overlaps.iter().zip(&mut used) {
            if weak.binary_search(&o.query).is_ok() || weak.binary_search(&o.target).is_ok() {
                *u = Usage::Ignore;
            }
        }
        log::info!("bacteria: {} reads of shallow contigs left out", weak.len());
        layout = lay_out(&laid, &overlaps, &used, &params.layout);
        log_layout(&layout);
    }
    let (kept_overlaps, laid_reads) = if params.keep_layout {
        let lens = laid.indices().map(|r| laid.read_len(r));
        (
            overlaps.into_iter().zip(used).collect(),
            selection.kept.iter().copied().zip(lens).collect(),
        )
    } else {
        // Not needed past the layout: memory for the polishing.
        drop((overlaps, used));
        (Vec::new(), Vec::new())
    };
    drop(laid);
    for c in layout.contigs.iter().filter(|c| c.seq.len() >= 20_000) {
        log::debug!(
            "bacteria: contig of {} bases, {} reads, first read {}, stops {:?}",
            c.seq.len(),
            c.reads.len(),
            c.reads[0].read.0,
            c.stops
        );
    }
    let mut drafts: Vec<Draft> = layout
        .contigs
        .iter()
        .map(|c| Draft::new(c.seq.clone(), c.circular))
        .collect();
    if !params.keep_layout {
        layout = Layout::default();
    }
    // Ends are joined after the first round, once the drafts no longer hold the raw reads'
    // indels, and the joins (and the circles they close) are polished by the rounds after.
    // Overlaps that reach the pieces' ends: with clipped ones, the first chain of a pair to
    // align could be a repeat's copy inside both (an rRNA operon at a circle's join), and hide
    // the overlap.
    let join_params = OverlapParams {
        extension: Extension::Full,
        chains_per_pair: JOIN_CHAINS_PER_PAIR,
        ..params.map
    };
    let rounds = params.polish_rounds.max(2);
    let mut covered = vec![0; drafts.len()];
    for round in 0..rounds {
        let started = Instant::now();
        if round == 1 {
            // A draft run through a repeat into another copy's flank, where reads say so.
            let rewired;
            (drafts, rewired) = repeat::rewire_at_repeats(drafts, reads, &params.map)?;
            let joined;
            (drafts, joined) = join::join_ends(drafts, &join_params)?;
            let (closed, left_out, crossed);
            (drafts, closed, left_out) = bridge::bridge(drafts, &params.map)?;
            let more;
            (drafts, more) = join::join_more(drafts, &join_params)?;
            // Ends that overlap only along copies of a repeat, by the reads across it; the ends
            // left may then join as before.
            let (across, after);
            (drafts, across) = join::join_across_repeats(drafts, reads, &join_params, &params.map)?;
            (drafts, after) = if across > 0 {
                join::join_more(drafts, &join_params)?
            } else {
                (drafts, 0)
            };
            let joined = joined + more + across + after;
            (drafts, crossed) = join::join_by_reads(drafts, reads, &params.map, &backend)?;
            let (paired, anchored, held);
            (drafts, paired) = join::join_by_anchors(drafts, &join_params)?;
            (drafts, anchored) = join::close_by_anchors(drafts, &join_params)?;
            let joined = joined + paired;
            // Drafts whose end runs on into their start along a stretch the joins could not
            // take (one holding an inverted segment): circles where reads run across.
            let self_closed;
            (drafts, self_closed) = circle::close_self_overlaps(drafts, reads, &params.map)?;
            // Once the last drafts are closed, what they hold of the linear ones left (a
            // variant's stretch, laid out with the bases either side of it) is left out too.
            (drafts, held) = bridge::without_contained(drafts, &params.map, bridge::FINE_CHUNK)?;
            // Linear replicons read round their hairpin ends, laid out as circles twice their
            // length or with their ends folded back, cut to what they are.
            let unfolded;
            (drafts, unfolded) = hairpin::unfold(drafts, &params.map)?;
            let left_out = left_out + held;
            // Plasmids shorter than their reads lay out as nothing: found in the reads that
            // hold them whole, once the drafts are polished and their tandems collapsed (so
            // the reads of those they hold map onto them).
            let circles = small::circles_in_reads(reads, &drafts, &params.map, &params.overlap)?;
            let found = circles.len();
            for c in circles {
                match c.replaces {
                    Some(d) => drafts[d] = c.draft,
                    None => drafts.push(c.draft),
                }
            }
            log::info!(
                "bacteria: {rewired} drafts rewired at repeats, {joined} joins of draft ends, \
                 {closed} drafts closed by another's \
                 bases, {left_out} held by others left out, {crossed} gaps crossed by reads, \
                 {anchored} closed through end anchors, {} closed and {} trimmed along their own \
                 ends, {unfolded} cut at hairpin ends, {found} \
                 small circles from whole reads; \
                 {} drafts ({} circular)",
                self_closed.circles,
                self_closed.trimmed,
                drafts.len(),
                drafts.iter().filter(|d| d.circular).count()
            );
        }
        if round > 0 {
            for d in &mut drafts {
                d.rotate(d.seq.len() / 2);
            }
        }
        let stats = polish_round(&mut drafts, reads, params, &backend)?;
        log::info!(
            "bacteria: polishing round {}: {} reads mapped, {} of {} windows polished, {:.1} s",
            round + 1,
            stats.mapped_reads,
            stats.polished,
            stats.windows,
            started.elapsed().as_secs_f64()
        );
        covered = stats.covered;
    }
    // Small plasmids read round more than once, which the first round's copies hid: one circle
    // each, polished across its new junction.
    let collapsed;
    (drafts, collapsed) = join::collapse_tandems(drafts);
    if collapsed > 0 {
        log::info!("bacteria: {collapsed} polished linear drafts are circles read round");
        covered = polish_round(&mut drafts, reads, params, &backend)?.covered;
    }
    // Once more, polished like the rest: a line laid out from a replicon's reads beside its
    // circle (*Borrelia* cp32-4, 30.3 kb twice) matches it only once both are polished.
    let kept = bridge::not_contained(&drafts, &params.map, bridge::FINE_CHUNK)?;
    if kept.contains(&false) {
        log::info!(
            "bacteria: {} polished drafts the others hold left out",
            kept.iter().filter(|&&k| !k).count()
        );
        drafts = kept_only(drafts, &kept);
        covered = kept_only(covered, &kept);
    }

    // Loci where most reads have another version; a replaced one is the consensus of a few
    // reads, so the contigs are polished once more.
    let mut mappings = None;
    if params.repair_loci {
        let started = Instant::now();
        let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
        let mut mapped = map_reads(&targets, reads, &params.map)?;
        // Tandem copies most reads have once first (the loci below are found from the reads'
        // mappings, which a left-out copy moves).
        let copies = locus::drop_minor_copies(&mut drafts, reads, &mapped, &params.map);
        // And tandem copies laid out as circles of their own put back.
        let (kept, restored) = locus::restore_copies(&mut drafts, reads, &mapped, &params.map)?;
        if restored > 0 {
            drafts = kept_only(drafts, &kept);
            covered = kept_only(covered, &kept);
            let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
            mapped = map_reads(&targets, reads, &params.map)?;
        }
        // Then a locus laid out as a line of its own beside the contig with the other version.
        let (kept, settled) = locus::settle_versions(&mut drafts, &mapped, &params.map)?;
        if settled > 0 {
            drafts = kept_only(drafts, &kept);
            covered = kept_only(covered, &kept);
        }
        let dropped = copies + restored + settled;
        if dropped > 0 {
            let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
            mapped = map_reads(&targets, reads, &params.map)?;
        }
        let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
        let circular: Vec<bool> = drafts.iter().map(|d| d.circular).collect();
        let with_parts = locus::with_parts(&targets, &circular, reads, &mapped, &params.map)?;
        let repaired = locus::repair(&mut drafts, reads, &with_parts, &params.map, &backend);
        log::info!(
            "bacteria: {copies} tandem copies most reads have once left out, {restored} put back; \
             {settled} loci laid out twice settled; {repaired} loci given the version most reads have, {:.1} s",
            started.elapsed().as_secs_f64()
        );
        if repaired + dropped > 0 {
            covered = polish_round(&mut drafts, reads, params, &backend)?.covered;
        } else {
            mappings = Some(mapped);
        }
    }

    // One more mapping of the reads, onto the polished contigs, for the base calls and the
    // contigs' checks.
    let mut mappings = match mappings {
        Some(m) => m,
        None => {
            let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
            map_reads(&targets, reads, &params.map)?
        }
    };
    // Then what few reads join to the rest of a linear contig at its ends.
    let weak = weak_ends(&drafts, &mappings, params.min_depth);
    let cut = cut_ends(&mut drafts, &weak);
    if cut > 0 {
        log::info!("bacteria: {cut} bases cut from linear contigs' ends few reads join on");
        let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
        mappings = map_reads(&targets, reads, &params.map)?;
    }
    let mappings = mappings;
    let mut checks = {
        let seqs: Vec<&[u8]> = drafts.iter().map(|d| d.seq.as_slice()).collect();
        check::check(&seqs, reads, &mappings)
    };
    report_read_divergence(&mappings, params.kind == ReadKind::Ont && !found_noisy);
    let homopolymers = call_bases(&mut drafts, &mut checks, reads, &mappings, params)?;

    let mut contigs: Vec<Contig> = drafts
        .into_iter()
        .zip(covered)
        .zip(checks)
        .map(|((d, covered), check)| {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a mean depth, not an exact count"
            )]
            let depth = covered as f64 / d.seq.len().max(1) as f64;
            Contig {
                seq: d.seq,
                circular: d.circular,
                hairpins: d.hairpins,
                depth,
                check,
            }
        })
        .collect();
    contigs.sort_by(|a, b| b.seq.len().cmp(&a.seq.len()).then(a.seq.cmp(&b.seq)));
    let (contigs, mut shallow) = drop_shallow(contigs);
    if !shallow.is_empty() {
        log::info!(
            "bacteria: {} shallow contigs left out ({} bases; depths {})",
            shallow.len(),
            shallow.iter().map(|c| c.seq.len()).sum::<usize>(),
            shallow
                .iter()
                .take(10)
                .map(|c| format!("{:.1}", c.depth))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    let (contigs, excised) = drop_excised(contigs, &params.map)?;
    if !excised.is_empty() {
        log::info!(
            "bacteria: {} shallow circles of a longer contig's bases left out ({} bases)",
            excised.len(),
            excised.iter().map(|c| c.seq.len()).sum::<usize>()
        );
        shallow.extend(excised);
    }
    Ok(Assembly {
        contigs,
        shallow,
        overlaps: kept_overlaps,
        layout,
        laid_reads,
        homopolymers,
    })
}

/// A linear contig's end is cut at a place fewer than [`BactParams::min_depth`] reads run
/// across, by this many bases either side...
const END_SPAN: u32 = 300;

/// ...within this many bases of the end...
const END_WINDOW: u32 = 10_000;

/// ...on a contig at least this deep. Reads there run across every place but a join that is
/// not in the genome.
const MIN_WEAK_END_DEPTH: f64 = 40.0;

/// The places in steps of this many bases are tried.
const END_STEP: u32 = 25;

/// Returns how many bases to cut off the start and the end of each draft: those of a linear
/// one's end outside the innermost place (within [`END_WINDOW`] of it) that fewer than
/// `min_depth` reads run across, by [`END_SPAN`] bases either side. Where few reads join a
/// line's end on, the bases past them are often another place's: with short reads (*A. tenax*,
/// N50 1.3 kb) one chimeric read laid out 1.1-1.3 kb of another place in the genome at three of
/// its four lines' ends, and on *A. loyolae* eight more. A stretch the reads cover thinly at a
/// line's end is cut as well, and the reads cannot tell it apart: on *C. difficile* NCTC 12727,
/// 6.8 kb of the genome at four lines' ends, 1-2 reads deep and held elsewhere too.
fn weak_ends(drafts: &[Draft], mappings: &[Option<Overlap>], min_depth: usize) -> Vec<[usize; 2]> {
    let mut spans: Vec<Vec<(u32, u32)>> = vec![Vec::new(); drafts.len()];
    for m in mappings.iter().flatten() {
        spans[m.target.index()].push((m.target_start, m.target_end));
    }
    drafts
        .iter()
        .zip(&spans)
        .map(|(d, spans)| {
            if d.circular || d.depth < MIN_WEAK_END_DEPTH {
                return [0; 2];
            }
            let len = u32::try_from(d.seq.len()).expect("contigs are shorter than 4 Gb");
            let window = END_WINDOW.min(len / 4);
            // Places measured from the start; the end's by flipping the spans.
            let flipped: Vec<(u32, u32)> = spans.iter().map(|&(s, e)| (len - e, len - s)).collect();
            let cut = |spans: &[(u32, u32)], hairpin: bool| -> usize {
                if hairpin {
                    return 0;
                }
                let near: Vec<(u32, u32)> = spans
                    .iter()
                    .copied()
                    .filter(|&(s, _)| s < window + END_SPAN)
                    .collect();
                let across = |at: u32| {
                    near.iter()
                        .filter(|&&(s, e)| s <= at.saturating_sub(END_SPAN) && e >= at + END_SPAN)
                        .count()
                };
                // The innermost place few reads run across, and the bases outside it.
                (0..=window / END_STEP)
                    .rev()
                    .map(|i| i * END_STEP)
                    .find(|&at| across(at) < min_depth)
                    .map_or(0, |at| at.saturating_sub(END_SPAN) as usize)
            };
            [cut(spans, d.hairpins[0]), cut(&flipped, d.hairpins[1])]
        })
        .collect()
}

/// Cuts `cuts` bases off the start and the end of each draft, unless they would take it all;
/// returns how many bases were cut.
fn cut_ends(drafts: &mut [Draft], cuts: &[[usize; 2]]) -> usize {
    let mut cut = 0;
    for (d, &[start, end]) in drafts.iter_mut().zip(cuts) {
        let start = if d.hairpins[0] { 0 } else { start };
        let end = if d.hairpins[1] { 0 } else { end };
        if start + end == 0 || start + end >= d.seq.len() {
            continue;
        }
        d.seq.truncate(d.seq.len() - end);
        d.seq.drain(..start);
        d.shallow_ends = [0; 2];
        cut += start + end;
    }
    cut
}

/// The items of `v` whose `kept` is true.
fn kept_only<T>(v: Vec<T>, kept: &[bool]) -> Vec<T> {
    v.into_iter()
        .zip(kept)
        .filter_map(|(x, &k)| k.then_some(x))
        .collect()
}

/// Reads that differ from the contigs at more than this share of their aligned columns (the
/// median read's) are noisier than the reads picoasm is tuned for: R10.4.1 reads differ at 0.5-4%,
/// R9.4.1 reads at 4-8%, and at 8-16% (R9 reads called with a fast model, 1-3 kb long) the
/// overlaps between reads are missed and the assembly falls apart.
const NOISY_READ_DIVERGENCE: f64 = 0.06;

/// Logs how much the reads differ from the contigs, and warns if they are noisy (see
/// [`NOISY_READ_DIVERGENCE`]) and were laid out raw as R10.4.1 reads (`raw_ont`).
fn report_read_divergence(mappings: &[Option<Overlap>], raw_ont: bool) {
    let mut divergence: Vec<f64> = mappings
        .iter()
        .flatten()
        .filter(|m| m.columns() > 0)
        .map(|m| {
            #[expect(clippy::cast_precision_loss, reason = "a share, not an exact count")]
            let share = m.differences() as f64 / m.columns() as f64;
            share
        })
        .collect();
    if divergence.is_empty() {
        return;
    }
    let mid = divergence.len() / 2;
    let (_, median, _) = divergence.select_nth_unstable_by(mid, f64::total_cmp);
    let median = *median;
    log::info!(
        "bacteria: the median read differs from the contigs at {:.1}% of its columns",
        100.0 * median
    );
    if raw_ont && median > NOISY_READ_DIVERGENCE {
        log::warn!(
            "bacteria: the reads are noisy (the median read differs from the contigs at {:.1}%; \
             R10.4.1 reads called with a sup or hac model differ at 0.5-4%): overlaps between reads \
             are missed, so expect gaps and a fragmented assembly; `--reads ont-noisy` corrects \
             the reads first",
            100.0 * median
        );
    }
}

/// Applies `changes` (per draft, sorted by start) to `drafts` and moves the places their
/// `checks` flag with the bases; returns how many changes were applied.
fn apply(drafts: &mut [Draft], checks: &mut [Check], changes: &[Vec<vote::Change>]) -> usize {
    let applied = vote::apply(drafts, changes);
    for (check, applied) in checks.iter_mut().zip(&applied) {
        if !applied.is_empty() && !check.flagged.is_empty() {
            let moves = vote::Moves::new(applied);
            check.lift(|at| moves.at(at));
        }
    }
    applied.iter().map(Vec::len).sum()
}

/// Calls the drafts' bases again from the reads (`mappings` onto them): homopolymer lengths and
/// per-strand votes, then versions by how the reads read each context (each round on a fresh
/// mapping). Returns the reads' lengths at the homopolymer runs if
/// [`BactParams::keep_layout`].
fn call_bases(
    drafts: &mut [Draft],
    checks: &mut [Check],
    reads: &Reads,
    mappings: &[Option<Overlap>],
    params: &BactParams,
) -> Result<Vec<Evidence>, AssembleError> {
    let mut homopolymers = Vec::new();
    if params.call_homopolymers || params.vote_margin.is_some() || params.keep_layout {
        // The votes are on bases outside runs of homopolymer::MIN_RUN or more, the homopolymer
        // calls on those runs, so their changes never meet.
        let started = Instant::now();
        let evidence = homopolymer::evidence(drafts, reads, mappings);
        let (mut changes, motifs): (Vec<Vec<vote::Change>>, usize) = match params.vote_margin {
            Some(margin) => {
                let modified = if params.learn_motifs {
                    Modified::Learn
                } else {
                    Modified::Ignore
                };
                vote::changes(drafts, reads, mappings, margin, modified)
            }
            None => (vec![Vec::new(); drafts.len()], 0),
        };
        let voted: usize = changes.iter().map(Vec::len).sum();
        if params.call_homopolymers {
            for (all, runs) in changes
                .iter_mut()
                .zip(homopolymer::changes(drafts, &evidence))
            {
                all.extend(runs);
            }
        }
        for c in &mut changes {
            c.sort_by_key(|c| c.start);
        }
        let total: usize = changes.iter().map(Vec::len).sum();
        let applied = apply(drafts, checks, &changes);
        log::info!(
            "bacteria: {} homopolymer runs measured; {motifs} misread motifs learned; {voted} \
             bases voted on, {} runs called again; {applied} changes, {:.1} s",
            evidence.len(),
            total - voted,
            started.elapsed().as_secs_f64()
        );
        if params.keep_layout {
            homopolymers = evidence;
        }
    }
    if let Some(min_gain) = params.context_gain {
        // Each round learns from the contigs the last one corrected: where most reads misread
        // a context, the first rounds learn it from the sites the contigs have right.
        for round in 0..context::ROUNDS {
            let started = Instant::now();
            let targets: Reads = drafts.iter().map(|d| d.seq.as_slice()).collect();
            let mappings = map_reads(&targets, reads, &params.map)?;
            let columns = vote::columns(drafts, reads, &mappings);
            let changes = context::changes(
                drafts,
                &columns,
                reads,
                &mappings,
                min_gain,
                params.composition,
            );
            let applied = apply(drafts, checks, &changes);
            log::info!(
                "bacteria: context round {}: {applied} bases given the version the reads favour, \
                 {:.1} s",
                round + 1,
                started.elapsed().as_secs_f64()
            );
            if applied == 0 {
                break;
            }
        }
    }
    Ok(homopolymers)
}

/// A contig shallower than this share of the chromosome's depth (the longest contig's) is left
/// out of the assembly, as Autocycler and Myloasm (`--min_depth_rel 0.1`) leave them out: reads
/// of another sample, read through a barcode, assemble at a few percent of the depth.
pub const MIN_DEPTH_SHARE: f64 = 0.1;

/// A contig that is not circular and shallower than this share of the chromosome's depth is
/// left out: a whole replicon assembles as a circle, and on the benchmark the linear contigs
/// below a quarter of the chromosome's depth were all another sample's or a chimera's.
pub const MIN_LINEAR_DEPTH_SHARE: f64 = 0.25;

/// A circular contig below this share of the chromosome's depth, most of whose bases a longer
/// contig holds, is left out: an element that excises from the chromosome (or a plasmid) and
/// closes into a circle in part of the cells, as the reference assemblies leave it out.
pub const MAX_EXCISED_DEPTH_SHARE: f64 = 0.5;

/// Splits `contigs` (longest first) into those kept and the shallow circles a longer one holds
/// (see [`MAX_EXCISED_DEPTH_SHARE`]).
fn drop_excised(
    contigs: Vec<Contig>,
    map: &OverlapParams,
) -> Result<(Vec<Contig>, Vec<Contig>), InvalidOverlapParamsError> {
    let Some(main) = contigs.first().map(|c| c.depth) else {
        return Ok((contigs, Vec::new()));
    };
    let mut excised = vec![false; contigs.len()];
    for i in 1..contigs.len() {
        let c = &contigs[i];
        if !c.circular || c.depth >= MAX_EXCISED_DEPTH_SHARE * main {
            continue;
        }
        let longer: Reads = contigs[..i]
            .iter()
            .zip(&excised)
            .filter(|&(l, &gone)| !gone && l.seq.len() > c.seq.len())
            .map(|(l, _)| l.seq.as_slice())
            .collect();
        excised[i] = bridge::held(&c.seq, bridge::Shape::Circular, &longer, map)?;
    }
    let (kept, gone): (Vec<_>, Vec<_>) = contigs.into_iter().zip(excised).partition(|(_, e)| !e);
    Ok((
        kept.into_iter().map(|(c, _)| c).collect(),
        gone.into_iter().map(|(c, _)| c).collect(),
    ))
}

/// Splits `contigs` (longest first) into those kept and those too shallow beside the first.
///
/// Depths are the contigs' median bin depths ([`Check::median_depth`]). A mean is raised by a
/// high-copy element inside the contig: on a *Streptomyces* whose reads fall from 25x to 10x
/// along the chromosome, a 3.4 Mb contig holding 250 kb at 238x had a mean of 54x, and the
/// other 4.4 Mb of the chromosome, at 9-13x, fell under a quarter of it.
fn drop_shallow(contigs: Vec<Contig>) -> (Vec<Contig>, Vec<Contig>) {
    let Some(main) = contigs.first().map(|c| c.check.median_depth) else {
        return (contigs, Vec::new());
    };
    contigs.into_iter().partition(|c| {
        let share = if c.is_whole() {
            MIN_DEPTH_SHARE
        } else {
            MIN_LINEAR_DEPTH_SHARE
        };
        c.check.median_depth >= share * main
    })
}

/// Draft ends are overlapped trying this many chains per pair: 150 kb ends can share several
/// copies of an rRNA operon, each a chain inside both, before the overlap itself.
const JOIN_CHAINS_PER_PAIR: usize = 50;

/// A contig that is not circular and whose depth is below this share of the deepest long
/// contig's (the chromosome's) is shallow: its reads are left out and the rest laid out again.
const WEAK_DEPTH_SHARE: f64 = 0.25;

/// A shallow contig's read, or a lone one, is left out only if at least this share of the reads
/// it overlaps are at home in a deeper contig. *E. coli* H10407's 5,538 bp plasmid kept its
/// reads at 0.05 to 0.5; at 0.5 the reads of *B. burgdorferi*'s cp32 plasmids, a family of near
/// copies, stayed too and two cp32s came out wrong (one three times over), but not at 0.2 or
/// below.
const WEAK_PARTNER_SHARE: f64 = 0.1;

/// Returns the reads of the layout's shallow contigs (see [`WEAK_DEPTH_SHARE`]), and those left
/// alone, sorted, that overlap (by the overlaps `used`) the deeper contigs' reads (see
/// [`WEAK_PARTNER_SHARE`]): those of a replicon of its own overlap almost only each other
/// (*E. coli* H10407's 5.5 kb plasmid at a fifth of the chromosome's depth, in reads shorter
/// than it, laid out as a linear contig and lone reads, and tied to the chromosome only by two
/// reads of a repeat they share). The laid-out reads are mapped onto the contigs, each to its
/// best place (its home), and a contig's depth is the bases mapped onto it over its length.
fn weak_reads(
    layout: &Layout,
    laid: &Reads,
    overlaps: &[Overlap],
    used: &[Usage],
    map: &OverlapParams,
) -> Result<Vec<ReadIndex>, InvalidOverlapParamsError> {
    if layout.contigs.len() < 2 {
        return Ok(Vec::new());
    }
    let targets: Reads = layout.contigs.iter().map(|c| c.seq.as_slice()).collect();
    let mut covered = vec![0_u64; targets.len()];
    let mut home = vec![None; laid.len()];
    for (r, m) in map_reads(&targets, laid, map)?.into_iter().enumerate() {
        if let Some(m) = m {
            covered[m.target.index()] += u64::from(m.target_end - m.target_start);
            home[r] = Some(m.target.index());
        }
    }
    #[expect(clippy::cast_precision_loss, reason = "a mean depth")]
    let depth: Vec<f64> = layout
        .contigs
        .iter()
        .zip(&covered)
        .map(|(c, &b)| b as f64 / c.seq.len().max(1) as f64)
        .collect();
    // The chromosome: the longest contig.
    let main = layout
        .contigs
        .iter()
        .zip(&depth)
        .max_by_key(|(c, _)| c.seq.len())
        .map_or(0.0, |(_, &d)| d);
    let shallow = |circular: bool, d: f64| !circular && d < WEAK_DEPTH_SHARE * main;
    let mut weak: Vec<ReadIndex> = layout
        .contigs
        .iter()
        .zip(&depth)
        .filter(|&(c, &d)| shallow(c.circular, d))
        .flat_map(|(c, _)| c.reads.iter().map(|p| p.read))
        .collect();
    // Reads left alone (a contig of one read) are no contig either, but may join one: a
    // single read beside a circle's junction keeps it from closing.
    weak.extend(
        layout
            .fates
            .iter()
            .zip(laid.indices())
            .filter(|&(&f, _)| f == Fate::Alone)
            .map(|(_, r)| r),
    );
    weak.sort_unstable();
    weak.dedup();
    // Each read's overlap partners, and how many of them are at home in a deeper contig.
    let deep: Vec<bool> = layout
        .contigs
        .iter()
        .zip(&depth)
        .map(|(c, &d)| !shallow(c.circular, d))
        .collect();
    let is_deep = |r: ReadIndex| home[r.index()].is_some_and(|c| deep[c]);
    let mut partners = vec![[0_u32; 2]; laid.len()];
    for (o, &u) in overlaps.iter().zip(used) {
        if u != Usage::Ignore {
            partners[o.query.index()][usize::from(is_deep(o.target))] += 1;
            partners[o.target.index()][usize::from(is_deep(o.query))] += 1;
        }
    }
    weak.retain(|r| {
        let [other, deeper] = partners[r.index()];
        f64::from(deeper) >= WEAK_PARTNER_SHARE * f64::from(other + deeper)
    });
    Ok(weak)
}

/// Logs what the layout did.
fn log_layout(layout: &Layout) {
    log::info!(
        "bacteria: layout: {} edges ({} transitive, {} short); {} chimeric, {} contained, {} in \
         tips, {} in bubbles, {} alone; {} contigs ({} circular), {} bases; ends: {} dead, {} \
         forks, {} joins",
        layout.edges,
        layout.transitive,
        layout.short,
        layout.chimeric,
        layout.contained,
        layout.tip_reads,
        layout.bubble_reads,
        layout.lone_reads,
        layout.contigs.len(),
        layout.contigs.iter().filter(|c| c.circular).count(),
        layout.contigs.iter().map(|c| c.seq.len()).sum::<usize>(),
        layout.ends.dead,
        layout.ends.forks,
        layout.ends.joins
    );
}

/// After correction, overlaps differing at more than this many times the median overlap's
/// share of columns are not edges (at most [`BactParams::edge_divergence`], at least
/// [`MIN_CORRECTED_EDGE_DIVERGENCE`]). Corrected reads of one place differ by little, and by
/// less the better the raw reads: on *S. turgidiscabies* R9 reads (7% off) the median overlap
/// differed at 0.04%, on *A. baumannii* (18% off) at 0.19%. At the raw reads' 12% the overlaps
/// that are not of one place made 26,606 errors on the first; at 1% 13,168 (Flye and Medaka:
/// 18,080). A fixed 1% on the second made 34,594 errors where 12% made 20,247; 25 times its
/// median (4.7%) made 20,237.
const CORRECTED_EDGE_FACTOR: f64 = 25.0;

/// The least cut-off on corrected reads' edges (see [`CORRECTED_EDGE_FACTOR`]): at 0.08% (twice
/// the median) the *S. turgidiscabies* chromosome fell into 18 pieces.
const MIN_CORRECTED_EDGE_DIVERGENCE: f64 = 0.01;

/// Raw ONT reads missing overlaps through their errors have fewer than this many each (a read
/// overlaps about twice the depth's worth of others): *S. turgidiscabies*, *S. stelliscabiei*,
/// *V. parahaemolyticus* and *E. coli* ATCC 25922 1.6-2.7, every other set above 7.5% (see
/// `BactParams::noisy_divergence`) 9.9 or more but *A. loyolae* (1.5, its reads too short to
/// overlap).
const NOISY_OVERLAPS_PER_READ: f64 = 5.0;

/// Reads are corrected only if their N50 (of the reads laid out) is at least this long:
/// correction from reads of about a kilobase made thousands of errors (*A. loyolae*, N50 1,134:
/// 250 -> 55,738), where at 2-2.7 kb it gained (*V. parahaemolyticus*, 1,962; *S.
/// turgidiscabies*, 2,668).
const MIN_NOISY_READ_N50: u32 = 1_500;

/// Whether raw ONT reads (`reads`, with their `overlaps`) are noisy enough to be corrected first:
/// their overlaps differ at a median above `max_divergence`, they miss overlaps (fewer than
/// [`NOISY_OVERLAPS_PER_READ`] each), and they are long enough to correct from
/// ([`MIN_NOISY_READ_N50`]).
fn looks_noisy(reads: &Reads, overlaps: &[Overlap], max_divergence: f64) -> bool {
    let median = median_divergence(overlaps);
    #[expect(clippy::cast_precision_loss, reason = "a ratio of counts")]
    let per_read = overlaps.len() as f64 / reads.len().max(1) as f64;
    let n50 = n50(reads);
    log::info!(
        "bacteria: overlaps between the reads differ at a median {:.1}% of their columns, {:.1} \
         per read; read N50 {n50}",
        100.0 * median,
        per_read
    );
    median > max_divergence && per_read < NOISY_OVERLAPS_PER_READ && n50 >= MIN_NOISY_READ_N50
}

/// The length such that reads at least as long hold half of `reads`' bases (0 if none).
fn n50(reads: &Reads) -> u32 {
    let mut lens: Vec<u32> = reads.indices().map(|r| reads.read_len(r)).collect();
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let half = reads.total_bases().div_ceil(2);
    let mut sum = 0;
    lens.into_iter()
        .find(|&len| {
            sum += len as usize;
            sum >= half
        })
        .unwrap_or(0)
}

/// The median share of columns at which `overlaps` differ (0 if there are none).
fn median_divergence(overlaps: &[Overlap]) -> f64 {
    #[expect(clippy::cast_precision_loss, reason = "shares of an overlap's columns")]
    let mut shares: Vec<f64> = overlaps
        .iter()
        .map(|o| o.differences() as f64 / o.columns().max(1) as f64)
        .collect();
    if shares.is_empty() {
        return 0.0;
    }
    let mid = shares.len() / 2;
    *shares.select_nth_unstable_by(mid, f64::total_cmp).1
}

/// Returns how the layout uses each overlap between raw reads: an edge where it differs at
/// no more than `max_divergence` of its columns.
fn edges(overlaps: &[Overlap], max_divergence: f64) -> Vec<Usage> {
    overlaps
        .iter()
        .map(|o| {
            #[expect(
                clippy::cast_precision_loss,
                reason = "a share of an overlap's columns"
            )]
            let d = o.differences() as f64 / o.columns().max(1) as f64;
            if d <= max_divergence {
                Usage::Edge
            } else {
                Usage::Ignore
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_are_noisy_only_if_their_overlaps_differ_much_are_few_and_the_reads_are_long() {
        use picoasm_overlap::{Edit, EditOp};
        use picoasm_seq::Strand;
        // An overlap differing at `differ` columns in a hundred.
        let overlap = |differ: u32| Overlap {
            query: ReadIndex(0),
            target: ReadIndex(1),
            strand: Strand::Forward,
            query_start: 0,
            query_end: 100,
            target_start: 0,
            target_end: 100,
            edits: vec![
                Edit::new(EditOp::Match, 100 - differ),
                Edit::new(EditOp::Mismatch, differ),
            ],
        };
        let long = vec![b'A'; 3_000];
        let short = vec![b'A'; 1_000];
        let reads = |read: &[u8]| -> Reads { std::iter::repeat_n(read, 10).collect() };
        let (long_reads, short_reads) = (reads(&long), reads(&short));
        let few_noisy = vec![overlap(11); 20];
        assert!(looks_noisy(&long_reads, &few_noisy, 0.075));
        // Clean overlaps, many of them, or reads too short to correct from: laid out raw.
        assert!(!looks_noisy(&long_reads, &vec![overlap(3); 20], 0.075));
        assert!(!looks_noisy(&long_reads, &vec![overlap(11); 100], 0.075));
        assert!(!looks_noisy(&short_reads, &few_noisy, 0.075));
    }

    #[test]
    fn half_the_bases_are_in_reads_at_least_the_n50_long() {
        let lens = [5, 1, 2, 2];
        let seqs: Vec<Vec<u8>> = lens.iter().map(|&n| vec![b'A'; n]).collect();
        let reads: Reads = seqs.iter().map(Vec::as_slice).collect();
        assert_eq!(n50(&reads), 5);
        assert_eq!(n50(&Reads::default()), 0);
    }

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
    fn a_tail_few_reads_join_on_is_cut_from_a_deep_line() {
        use picoasm_overlap::{Edit, EditOp};
        use picoasm_seq::Strand;
        let span = |target: u32, start: u32, end: u32| {
            Some(Overlap {
                query: ReadIndex(0),
                target: ReadIndex(target),
                strand: Strand::Forward,
                query_start: 0,
                query_end: end - start,
                target_start: start,
                target_end: end,
                edits: vec![Edit::new(EditOp::Match, end - start)],
            })
        };
        // A 20 kb line whose first 1.2 kb only reads of another copy cover (they stop there)
        // and one chimeric read joins on; the rest held by 1 kb reads every 10 bases, as is
        // its end. A circle and a shallow line with the same reads keep their bases.
        let mut mappings = Vec::new();
        for target in 0..3 {
            for s in (0..200).step_by(5) {
                mappings.push(span(target, s, 1_200));
            }
            mappings.push(span(target, 0, 2_500));
            for s in (1_200..19_000).step_by(10) {
                mappings.push(span(target, s, s + 1_000));
            }
            for s in (19_000..20_000).step_by(10) {
                mappings.push(span(target, s, 20_000));
            }
        }
        let mut line = Draft::new(random(20_000, 3), false);
        line.depth = 100.0;
        let mut circle = line.clone();
        circle.circular = true;
        let mut shallow = line.clone();
        shallow.depth = 10.0;
        let drafts = [line, circle, shallow];

        assert_eq!(
            weak_ends(&drafts, &mappings, 3),
            [[1_200, 0], [0, 0], [0, 0]]
        );
    }

    fn contig(seq: &[u8], depth: f64) -> Contig {
        Contig {
            seq: seq.to_vec(),
            circular: true,
            hairpins: [false; 2],
            depth,
            check: Check::default(),
        }
    }

    #[test]
    fn a_high_copy_element_does_not_make_the_rest_of_the_chromosome_shallow() {
        // Lines of one chromosome whose depth falls from 25x to 10x: the longest holds a
        // high-copy element that raises its mean to 54x, but its median is 20x.
        let line = |len: usize, mean: f64, median: f64| Contig {
            circular: false,
            check: Check {
                median_depth: median,
                ..Check::default()
            },
            ..contig(&random(len, 3), mean)
        };
        let contigs = vec![
            line(3_000, 54.0, 20.0),
            line(2_000, 12.0, 12.0),
            line(1_000, 10.0, 9.0),
            // Another sample's reads.
            line(500, 3.0, 3.0),
        ];
        let (kept, shallow) = drop_shallow(contigs);
        assert_eq!(kept.len(), 3);
        assert_eq!(shallow.len(), 1);
        assert_eq!(shallow[0].seq.len(), 500);
    }

    #[test]
    fn shallow_circles_of_the_chromosomes_bases_are_left_out() {
        let chromosome = random(100_000, 1);
        let contigs = vec![
            contig(&chromosome, 40.0),
            // Shallow and the chromosome's: an excised element.
            contig(&chromosome[20_000..30_000], 10.0),
            // Another, 13 kb, cut where most of a 5 kb piece runs across its junction.
            contig(
                &[&chromosome[60_516..63_005], &chromosome[50_000..60_516]].concat(),
                10.0,
            ),
            // As deep as half the chromosome: kept.
            contig(&chromosome[50_000..60_000], 30.0),
            // Shallow but its own bases: a plasmid.
            contig(&random(10_000, 2), 10.0),
        ];
        let map = BactParams::new(ReadKind::Ont).map;
        let (kept, excised) = drop_excised(contigs, &map).expect("valid settings");
        assert_eq!(kept.len(), 3);
        assert_eq!(excised.len(), 2);
        assert_eq!(excised[0].seq, chromosome[20_000..30_000]);
    }
}
