//! Settings for [`assemble`](crate::assemble), by kind of read.

use crate::consensus::PoaParams;
use crate::overlap::{Extension, LayoutParams, OverlapParams};

/// The kind of long read assembled.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReadKind {
    /// Oxford Nanopore R10.4.1 reads (Dorado sup or hac): 1-3% errors, mostly homopolymer
    /// lengths. Reads found noisier are corrected first, as [`ReadKind::OntNoisy`]'s (see
    /// [`BactParams::noisy_divergence`]).
    Ont,
    /// Noisier Oxford Nanopore reads (R9.4.1 or older, or a fast model): 5-18% errors. The
    /// reads are corrected from each other before they are laid out.
    OntNoisy,
    /// PacBio HiFi reads: about one error per thousand bases.
    Hifi,
}

/// Pieces of an Oxford Nanopore overlap between seeds may differ at this share of their columns
/// (see [`OverlapParams::max_piece_divergence`]): stretches of a read with errors at 20% or
/// more lost overlaps of tens of kilobases, and with them contigs' joins.
const MAX_PIECE_DIVERGENCE: f64 = 0.35;

/// Overlaps between Oxford Nanopore reads chain across indels of up to about this many bases
/// (see [`OverlapParams::chain_slack`]): a burst of errors in one read breaks a chain at 50,
/// and the overlap with it.
const OVERLAP_CHAIN_SLACK: u32 = 200;

/// Reads mapped onto a draft chain across indels of up to about this many bases (see
/// [`OverlapParams::chain_slack`]): a draft spelled from raw reads can lack, or hold twice, a few
/// hundred bases where one read hands over to the next, and polishing mends that only if the
/// reads' alignments run across it.
const MAP_CHAIN_SLACK: u32 = 2_000;

/// Noisy Oxford Nanopore reads are corrected in this many rounds. On *A. baumannii* R9 reads
/// (82% identity) one round left 65,629 errors in the assembly and two 20,247 (Flye and Medaka:
/// 21,235, and one of the three replicons complete where these complete all three).
const NOISY_CORRECTION_ROUNDS: u32 = 2;

/// Oxford Nanopore reads whose overlaps differ at a median above this share of their columns may
/// be noisy (see [`BactParams::noisy_divergence`]; the reads must also miss overlaps and be long
/// enough to correct from). Over 120 benchmark sets the median is 1.9-5% for R10.4.1 reads
/// called with sup or hac, 6-9% for R9.4.1 reads and R10.4.1 reads of heavily methylated
/// genomes, and 8-11.5% for reads 8-20% off; it cannot go much higher, as overlaps over 12% are
/// not found. Corrected first, *A. baumannii* R9 went 0/3 -> 3/3 complete, *E. coli* ATCC 25922
/// (8.0%) 4/5 -> 5/5 and 757 -> 300 errors, *V. parahaemolyticus* 13,029 -> 1,613 errors, *S.
/// turgidiscabies* 27,092 -> 12,526 and *S. stelliscabiei* 1.05 -> 0.52 Mb missing; above 7.5%
/// the sets that still find many overlaps per read changed little or lost (*A. mictus*, 10.0%:
/// 7,844 -> 13,044 errors).
const NOISY_OVERLAP_DIVERGENCE: f64 = 0.075;

/// Settings for [`assemble`](crate::assemble).
///
/// Start from [`BactParams::new`] for the kind of read, then adjust fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BactParams {
    /// The kind of read.
    pub kind: ReadKind,
    /// Overlaps between reads, for the layout.
    pub overlap: OverlapParams,
    /// Overlaps differing at no more than this share of their columns are the string graph's
    /// edges (ONT; HiFi reads' edges are those the correction finds consistent).
    pub edge_divergence: f64,
    /// The string graph.
    pub layout: LayoutParams,
    /// Rounds of read correction before the layout (HiFi and noisy ONT; R10.4.1 reads are laid
    /// out raw).
    pub correction_rounds: u32,
    /// ONT reads laid out raw are corrected first after all (in as many rounds as
    /// [`ReadKind::OntNoisy`]'s) if their overlaps differ at a median above this share of their
    /// columns (7.5% for [`ReadKind::Ont`]), they find few overlaps each and they are not too
    /// short (`None`: never).
    pub noisy_divergence: Option<f64>,
    /// Reads mapped onto the contigs, for polishing.
    pub map: OverlapParams,
    /// Rounds of polishing; at least two are run (draft ends are joined after the first).
    pub polish_rounds: u32,
    /// Polishing windows are about this many contig bases long.
    pub window_len: usize,
    /// A window is polished from at most this many read segments.
    pub max_depth: usize,
    /// A window with fewer read segments keeps its bases.
    pub min_depth: usize,
    /// Partial-order alignment scoring.
    pub poa: PoaParams,
    /// Reads deeper than this many times the genome's depth (a high-copy plasmid's) are
    /// thinned to it before the layout (see `select.rs`).
    pub max_depth_ratio: f64,
    /// Where the genome is deeper than this (as `select.rs` counts it), only the longest reads
    /// enough for this depth are laid out.
    pub max_layout_depth: f64,
    /// Homopolymer runs are called again from the reads after polishing, by how the reads
    /// give runs of each length (see `homopolymer.rs`); `false` leaves them as polishing called
    /// them.
    pub call_homopolymers: bool,
    /// After polishing, loci where most reads have another version than the contig (a
    /// segment turned round, or another sequence) get that version, and are polished again
    /// (see `locus.rs`).
    pub repair_loci: bool,
    /// After polishing, the reads vote on each base outside runs, per strand, and a base is
    /// replaced, removed or added when they favour that by this margin (see `vote.rs`); `None`
    /// leaves the bases as polishing called them.
    pub vote_margin: Option<f64>,
    /// With the vote, learn the motifs one strand's reads misread (modified bases) and take the
    /// other strand's base where a contig has the misread one (see `motif.rs`).
    pub learn_motifs: bool,
    /// After the votes, take another version of a base where it makes the reads at least this
    /// much likelier (a natural log), by how they read each context (see `context.rs`); `None`
    /// to skip.
    pub context_gain: Option<f64>,
    /// What the context step charges a substitution that makes the genome's 6-mers round it
    /// rarer: times the fall in their log counts (see `context.rs`); 0 for nothing.
    pub composition: f64,
    /// Keep the overlaps, the layout and the reads' homopolymer lengths in the
    /// [`Assembly`](crate::Assembly) (for development).
    pub keep_layout: bool,
}

impl BactParams {
    /// Returns the settings for `kind` of reads.
    #[must_use]
    pub fn new(kind: ReadKind) -> Self {
        // Raw reads' overlaps differ at a few percent, and their snapped overhangs put a read's
        // place off by up to hang::MAX_HANG bases.
        let layout = LayoutParams {
            max_divergence: 1.0,
            fuzz: 1_000,
            max_tip_reads: 6,
            ..LayoutParams::default()
        };
        match kind {
            ReadKind::Ont | ReadKind::OntNoisy => Self {
                kind,
                overlap: OverlapParams {
                    kmer_len: 19,
                    density: 0.05,
                    max_divergence: 0.12,
                    extension: Extension::Clipped,
                    max_piece_divergence: Some(MAX_PIECE_DIVERGENCE),
                    chain_slack: OVERLAP_CHAIN_SLACK,
                    ..OverlapParams::default()
                },
                edge_divergence: 0.12,
                layout,
                correction_rounds: if kind == ReadKind::OntNoisy {
                    NOISY_CORRECTION_ROUNDS
                } else {
                    0
                },
                noisy_divergence: (kind == ReadKind::Ont).then_some(NOISY_OVERLAP_DIVERGENCE),
                map: OverlapParams {
                    kmer_len: 19,
                    density: 0.05,
                    max_divergence: 0.15,
                    max_piece_divergence: Some(MAX_PIECE_DIVERGENCE),
                    chain_slack: MAP_CHAIN_SLACK,
                    min_overlap: 500,
                    max_candidates: 20,
                    ..OverlapParams::default()
                },
                polish_rounds: 2,
                window_len: 500,
                max_depth: 40,
                min_depth: 3,
                poa: PoaParams::default(),
                max_depth_ratio: crate::select::MAX_DEPTH_RATIO,
                max_layout_depth: crate::select::MAX_LAYOUT_DEPTH,
                repair_loci: true,
                call_homopolymers: true,
                vote_margin: Some(crate::vote::MARGIN),
                learn_motifs: true,
                context_gain: Some(crate::context::MIN_GAIN),
                composition: crate::context::COMPOSITION_WEIGHT,
                keep_layout: false,
            },
            ReadKind::Hifi => Self {
                kind,
                overlap: OverlapParams::default(),
                edge_divergence: crate::overlap::MAX_DIVERGENCE,
                layout: LayoutParams::default(),
                correction_rounds: 2,
                noisy_divergence: None,
                map: OverlapParams {
                    min_overlap: 500,
                    max_candidates: 20,
                    chain_slack: MAP_CHAIN_SLACK,
                    ..OverlapParams::default()
                },
                polish_rounds: 2,
                window_len: 500,
                // As deep as ONT's, and runs called again from the reads: partial-order
                // alignment spreads a run's extra base over several columns, so HiFi reads
                // of 6 to 9 bases mostly right still came out a base off (L. monocytogenes'
                // T6 from 8 of 9 reads). Over the 25 HiFi sets, 8 such errors fewer, none
                // more, for 8% more time.
                max_depth: 40,
                min_depth: 3,
                poa: PoaParams::default(),
                max_depth_ratio: crate::select::MAX_DEPTH_RATIO,
                max_layout_depth: crate::select::MAX_LAYOUT_DEPTH,
                repair_loci: true,
                call_homopolymers: true,
                vote_margin: None,
                learn_motifs: false,
                context_gain: None,
                composition: 0.0,
                keep_layout: false,
            },
        }
    }
}
