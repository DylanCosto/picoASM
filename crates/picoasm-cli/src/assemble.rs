//! `picoasm assemble` and `picoasm polish-short`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Args, ValueEnum};
use picoasm::{
    Assembly, BactParams, Check, Contig, Flag, FlaggedBin, ReadKind, ShortParams, assemble,
    polish_short,
};
use picoasm_io::{
    FastaWriter, GfaLink, GfaSegment, GfaTag, GfaWriter, OrientedSegment, Overlap as GfaOverlap,
    ReadRecord, ReadSource,
};
use picoasm_overlap::Reads;
use picoasm_seq::Strand;

#[derive(Debug, Args)]
pub(crate) struct AssembleArgs {
    /// Input read files (FASTA, FASTQ, gzip or BGZF, or unaligned BAM).
    #[arg(required = true, value_name = "READS")]
    reads: Vec<PathBuf>,

    /// Output FASTA of the contigs.
    #[arg(short, long, value_name = "FASTA")]
    output: PathBuf,

    /// The kind of reads.
    #[arg(long = "reads", value_enum, value_name = "KIND", default_value_t = Kind::Ont)]
    kind: Kind,

    /// Worker threads [default: all available cores].
    #[arg(short, long, value_name = "N")]
    threads: Option<NonZeroUsize>,

    /// Rounds of polishing [default: 2].
    #[arg(long, value_name = "N")]
    polish_rounds: Option<u32>,

    /// Also write the contigs left out as too shallow beside the chromosome (below 10% of its
    /// depth, or 25% if not circular: another sample's reads, or a chimera's) to this FASTA.
    #[arg(long, value_name = "FASTA")]
    shallow: Option<PathBuf>,

    /// Also write the places the reads show the contigs may be wrong to this BED file: stretches
    /// few reads cover (`shallow`), places many reads stop aligning (`stops`), and contigs not
    /// closed into a whole replicon (`linear`).
    #[arg(long, value_name = "BED")]
    suspect: Option<PathBuf>,

    /// Also write the contigs as a GFA graph (for Bandage and other graph tools): circular
    /// contigs link end to start, and hairpin ends link a contig to its reverse complement.
    #[arg(long, value_name = "GFA")]
    gfa: Option<PathBuf>,

    /// Short reads (Illumina) of the same isolate to polish the contigs with at the end; give
    /// the option once per file (both files of a pair, in any order).
    #[arg(long, value_name = "READS")]
    short_reads: Vec<PathBuf>,

    /// Settings for development, each overriding the default for the kind of reads.
    #[command(flatten)]
    dev: DevArgs,
}

/// Settings for tuning on the benchmark; hidden from `--help`.
#[derive(Debug, Args)]
struct DevArgs {
    /// Two paths to a read agree on where it starts within this many bases.
    #[arg(long, hide = true)]
    fuzz: Option<u32>,
    /// Overlaps differing at no more than this share of their columns are edges.
    #[arg(long, hide = true)]
    edge_divergence: Option<f64>,
    /// Overlaps differing at more than this share of their columns are not found.
    #[arg(long, hide = true)]
    max_divergence: Option<f64>,
    /// Rounds of read correction before the layout.
    #[arg(long, hide = true)]
    correction_rounds: Option<u32>,
    /// Minimizer length for the overlaps, in homopolymer-compressed bases.
    #[arg(long, hide = true)]
    kmer_len: Option<u8>,
    /// Share of compressed k-mers kept as minimizers, for the overlaps.
    #[arg(long, hide = true)]
    density: Option<f64>,
    /// Edges shorter than this share of a read's longest are dropped.
    #[arg(long, hide = true)]
    min_overlap_ratio: Option<f64>,
    /// Dead-end paths of at most this many reads are tips.
    #[arg(long, hide = true)]
    max_tip_reads: Option<usize>,
    /// Bubbles closing within this many bases are popped.
    #[arg(long, hide = true)]
    max_bubble_len: Option<u64>,
    /// Chained seeds of an overlap may leave their diagonal by this many bases.
    #[arg(long, hide = true)]
    chain_slack: Option<u32>,
    /// Leave homopolymer runs as polishing called them.
    #[arg(long, hide = true)]
    no_homopolymers: bool,
    /// Call homopolymer runs' lengths from the reads (the default for ONT reads).
    #[arg(long, hide = true, conflicts_with = "no_homopolymers")]
    homopolymers: bool,
    /// Thin reads deeper than this many times the genome's depth.
    #[arg(long, hide = true)]
    max_depth_ratio: Option<f64>,
    /// Lay out only the longest reads up to this depth.
    #[arg(long, hide = true)]
    max_layout_depth: Option<f64>,
    /// Leave loci where most reads have another version as the contigs have them.
    #[arg(long, hide = true)]
    no_repair_loci: bool,
    /// Do not learn the motifs one strand misreads.
    #[arg(long, hide = true)]
    no_motifs: bool,
    /// The least log-likelihood gain for the context step to change a base (negative: skip).
    #[arg(long, hide = true)]
    context_gain: Option<f64>,
    /// Weight of the genome's 6-mer composition at the context step's substitutions.
    #[arg(long, hide = true)]
    composition: Option<f64>,
    /// Margin of the reads' vote on each base after polishing (negative: no vote).
    #[arg(long, hide = true, allow_negative_numbers = true)]
    vote_margin: Option<f64>,
    /// Polishing window length.
    #[arg(long, hide = true)]
    window_len: Option<usize>,
    /// Read segments per polishing window, at most.
    #[arg(long, hide = true)]
    max_depth: Option<usize>,
    /// Write the overlaps, the draft contigs' reads and each read's fate to this directory.
    #[arg(long, hide = true, value_name = "DIR")]
    dump: Option<PathBuf>,
}

impl DevArgs {
    fn apply(&self, params: &mut BactParams) {
        if let Some(x) = self.fuzz {
            params.layout.fuzz = x;
        }
        if let Some(x) = self.edge_divergence {
            params.edge_divergence = x;
        }
        if let Some(x) = self.max_divergence {
            params.overlap.max_divergence = x;
        }
        if let Some(x) = self.correction_rounds {
            params.correction_rounds = x;
        }
        if let Some(x) = self.kmer_len {
            params.overlap.kmer_len = x;
        }
        if let Some(x) = self.density {
            params.overlap.density = x;
        }
        if let Some(x) = self.min_overlap_ratio {
            params.layout.min_overlap_ratio = x;
        }
        if let Some(x) = self.max_tip_reads {
            params.layout.max_tip_reads = x;
        }
        if let Some(x) = self.max_bubble_len {
            params.layout.max_bubble_len = x;
        }
        if let Some(x) = self.chain_slack {
            params.overlap.chain_slack = x;
        }
        if self.no_homopolymers {
            params.call_homopolymers = false;
        }
        if self.homopolymers {
            params.call_homopolymers = true;
        }
        if let Some(x) = self.max_depth_ratio {
            params.max_depth_ratio = x;
        }
        if let Some(x) = self.max_layout_depth {
            params.max_layout_depth = x;
        }
        if self.no_repair_loci {
            params.repair_loci = false;
        }
        if self.no_motifs {
            params.learn_motifs = false;
        }
        if let Some(x) = self.context_gain {
            params.context_gain = (x >= 0.0).then_some(x);
        }
        if let Some(x) = self.composition {
            params.composition = x;
        }
        if let Some(x) = self.vote_margin {
            params.vote_margin = (x >= 0.0).then_some(x);
        }
        if let Some(x) = self.window_len {
            params.window_len = x;
        }
        if let Some(x) = self.max_depth {
            params.max_depth = x;
        }
        params.keep_layout = self.dump.is_some();
    }
}

/// Returns `ok` for a contig the reads show nothing wrong with, else what they show, comma
/// separated: `linear`, `shallow:N` (bins of [`picoasm::CHECK_BIN`] bases under a quarter of its median
/// depth), `stops:N` (bins where many reads stop aligning).
fn check_tag(contig: &Contig) -> String {
    let check = &contig.check;
    if check.is_clean(contig.is_whole()) {
        return "ok".to_owned();
    }
    let mut tags = Vec::new();
    if !contig.is_whole() {
        tags.push("linear".to_owned());
    }
    if check.shallow_bins() > 0 {
        tags.push(format!("shallow:{}", check.shallow_bins()));
    }
    if check.stop_bins() > 0 {
        tags.push(format!("stops:{}", check.stop_bins()));
    }
    tags.join(",")
}

/// One region of [`suspect_regions`]: a run of adjacent flagged bins of one kind.
#[derive(Clone, Debug, PartialEq)]
struct Suspect {
    start: u32,
    end: u32,
    reason: String,
}

/// The places the reads show `contig` may be wrong, in order: adjacent bins of one kind
/// merged, a shallow stretch with its lowest depth, a place reads stop with its bin where most
/// reads stop (and that bin's depth); and the whole contig if it is not whole. A bin the last
/// base calls took every base out of has no place, and is left out.
fn suspect_regions(contig: &Contig) -> Vec<Suspect> {
    let len = u32::try_from(contig.seq.len()).expect("contigs are shorter than 4 Gb");
    let median = contig.check.median_depth;
    let mut out = Vec::new();
    if !contig.is_whole() {
        out.push(Suspect {
            start: 0,
            end: len,
            reason: "linear".to_owned(),
        });
    }
    for flag in [Flag::Shallow, Flag::Stop] {
        // Runs of adjacent bins, each with its shallowest bin and its bin where most reads stop.
        let mut runs: Vec<(u32, u32, FlaggedBin, FlaggedBin)> = Vec::new();
        for &bin in contig.check.flagged.iter().filter(|b| b.flag == flag) {
            match runs.last_mut() {
                Some(run) if run.1 == bin.start => {
                    run.1 = bin.end;
                    if bin.depth < run.2.depth {
                        run.2 = bin;
                    }
                    if bin.stops > run.3.stops {
                        run.3 = bin;
                    }
                }
                _ => runs.push((bin.start, bin.end, bin, bin)),
            }
        }
        for (start, end, shallowest, most_stops) in runs {
            let (start, end) = (start.min(len), end.min(len));
            if start >= end {
                continue;
            }
            let reason = match flag {
                Flag::Shallow => {
                    format!("shallow:depth={:.1},median={median:.1}", shallowest.depth)
                }
                Flag::Stop => format!(
                    "stops:reads={},depth={:.1}",
                    most_stops.stops, most_stops.depth
                ),
                _ => format!("{flag:?}"),
            };
            out.push(Suspect { start, end, reason });
        }
    }
    out.sort_by_key(|s| (s.start, s.end));
    out
}

/// Writes the places the reads show `contigs` may be wrong (see [`suspect_regions`]) as BED,
/// contigs named as in [`write_contigs`].
fn write_suspect(path: &Path, prefix: &str, contigs: &[Contig]) -> Result<()> {
    let write = || -> std::io::Result<()> {
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(out, "#contig\tstart\tend\treason")?;
        for (i, contig) in contigs.iter().enumerate() {
            for s in suspect_regions(contig) {
                writeln!(
                    out,
                    "{prefix}_{}\t{}\t{}\t{}",
                    i + 1,
                    s.start,
                    s.end,
                    s.reason
                )?;
            }
        }
        out.flush()
    };
    write().with_context(|| format!("writing {}", path.display()))
}

/// Writes `contigs` as GFA 1, named as in [`write_contigs`], each with its length and depth: a
/// circular contig links its end to its start, and a hairpin end links the contig to its own
/// reverse complement there.
fn write_gfa(path: &Path, prefix: &str, contigs: &[Contig]) -> Result<()> {
    let names: Vec<String> = (1..=contigs.len())
        .map(|i| format!("{prefix}_{i}"))
        .collect();
    let write = || -> Result<(), picoasm_io::WriteError> {
        let mut gfa = GfaWriter::create(path)?;
        gfa.write_header(&[])?;
        for (name, contig) in names.iter().zip(contigs) {
            let tags = [
                GfaTag::length(contig.seq.len() as u64),
                GfaTag::depth(contig.depth),
            ];
            gfa.write_segment(
                &GfaSegment::new(name)
                    .with_sequence(&contig.seq)
                    .with_tags(&tags),
            )?;
        }
        for (name, contig) in names.iter().zip(contigs) {
            let side = |strand| OrientedSegment::new(name, strand);
            let mut links = Vec::new();
            if contig.circular {
                links.push((Strand::Forward, Strand::Forward));
            }
            if contig.hairpins[1] {
                links.push((Strand::Forward, Strand::Reverse));
            }
            if contig.hairpins[0] {
                links.push((Strand::Reverse, Strand::Forward));
            }
            for (from, to) in links {
                gfa.write_link(&GfaLink::new(side(from), side(to), GfaOverlap::Matches(0)))?;
            }
        }
        gfa.finish()?;
        Ok(())
    };
    write().with_context(|| format!("writing {}", path.display()))
}

/// Writes `contigs` to a FASTA file at `path`, named `{prefix}_1` and on, with each one's
/// length, whether it is circular, its hairpin ends if it has any (`telomeres=hairpin` for
/// both, `hairpin=start` or `hairpin=end` for one) and its depth in the header.
fn write_contigs(path: &Path, prefix: &str, contigs: &[Contig]) -> Result<()> {
    let mut fasta =
        FastaWriter::create(path).with_context(|| format!("creating {}", path.display()))?;
    for (i, contig) in contigs.iter().enumerate() {
        let ends = match contig.hairpins {
            [true, true] => " telomeres=hairpin",
            [true, false] => " hairpin=start",
            [false, true] => " hairpin=end",
            [false, false] => "",
        };
        let header = format!(
            "{prefix}_{} length={} circular={}{ends} depth={:.1} check={}",
            i + 1,
            contig.seq.len(),
            contig.circular,
            contig.depth,
            check_tag(contig)
        );
        fasta
            .write_record(header, &contig.seq)
            .with_context(|| format!("writing {}", path.display()))?;
    }
    fasta
        .finish()
        .with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

/// Writes the overlaps, the draft contigs' reads, each read's fate and the reads' homopolymer
/// lengths to `dir`, as `overlaps.tsv`, `layout.tsv`, `fates.tsv` and `homopolymers.tsv`. Reads
/// are named by their index among the reads laid out; `laid.tsv` gives each one's index among
/// the input reads and its length once trimmed.
fn dump(dir: &Path, assembly: &Assembly) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    let len = |r: picoasm_overlap::ReadIndex| assembly.laid_reads[r.index()].1;
    let path = dir.join("laid.tsv");
    let mut out = BufWriter::new(
        File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    for (i, (input, trimmed)) in assembly.laid_reads.iter().enumerate() {
        writeln!(out, "{i}\t{}\t{trimmed}", input.0)
            .with_context(|| format!("writing {}", path.display()))?;
    }
    out.flush()
        .with_context(|| format!("writing {}", path.display()))?;
    let path = dir.join("overlaps.tsv");
    let mut out = BufWriter::new(
        File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    for (o, usage) in &assembly.overlaps {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{usage:?}",
            o.query.0,
            len(o.query),
            o.query_start,
            o.query_end,
            if o.strand.is_reverse() { '-' } else { '+' },
            o.target.0,
            len(o.target),
            o.target_start,
            o.target_end,
            o.columns(),
            o.differences(),
        )
        .with_context(|| format!("writing {}", path.display()))?;
    }
    out.flush()
        .with_context(|| format!("writing {}", path.display()))?;
    let path = dir.join("layout.tsv");
    let mut out = BufWriter::new(
        File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    for (i, c) in assembly.layout.contigs.iter().enumerate() {
        for p in &c.reads {
            writeln!(
                out,
                "{i}\t{}\t{}\t{}\t{:?}\t{}",
                p.read.0,
                if p.strand.is_reverse() { '-' } else { '+' },
                p.start,
                c.stops,
                c.seq.len()
            )
            .with_context(|| format!("writing {}", path.display()))?;
        }
    }
    out.flush()
        .with_context(|| format!("writing {}", path.display()))?;
    let path = dir.join("homopolymers.tsv");
    let mut out = BufWriter::new(
        File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    for e in &assembly.homopolymers {
        let counts: Vec<String> = e
            .counts
            .iter()
            .enumerate()
            .filter(|&(_, &n)| n > 0)
            .map(|(k, n)| format!("{k}:{n}"))
            .collect();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            e.run.contig,
            e.run.start,
            e.run.len,
            char::from(e.run.base),
            counts.join(",")
        )
        .with_context(|| format!("writing {}", path.display()))?;
    }
    out.flush()
        .with_context(|| format!("writing {}", path.display()))?;
    let path = dir.join("fates.tsv");
    let mut out = BufWriter::new(
        File::create(&path).with_context(|| format!("creating {}", path.display()))?,
    );
    for (read, fate) in assembly.layout.fates.iter().enumerate() {
        writeln!(out, "{read}\t{fate:?}").with_context(|| format!("writing {}", path.display()))?;
    }
    out.flush()
        .with_context(|| format!("writing {}", path.display()))
}

/// The kinds of reads, as the command line names them.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Kind {
    /// Oxford Nanopore R10.4.1 (Dorado sup or hac); reads found noisier are corrected first, as
    /// for ont-noisy.
    Ont,
    /// Oxford Nanopore R9.4.1 or older, or a fast model: reads 5-18% off, corrected first.
    OntNoisy,
    /// Oxford Nanopore reads laid out raw, even if they look noisy.
    OntClean,
    /// PacBio HiFi.
    Hifi,
}

/// Polish a bacterial assembly with short reads (Illumina) of the same isolate.
#[derive(Debug, Args)]
pub(crate) struct PolishShortArgs {
    /// The assembly (FASTA); a contig whose header has `circular=true` is a circle.
    #[arg(value_name = "ASSEMBLY")]
    assembly: PathBuf,

    /// The short reads (FASTQ or FASTA, gzip or not; both files of a pair, in any order).
    #[arg(required = true, value_name = "READS")]
    reads: Vec<PathBuf>,

    /// Output FASTA: the contigs polished, under their headers (`length=` updated).
    #[arg(short, long, value_name = "FASTA")]
    output: PathBuf,

    /// Worker threads [default: all available cores].
    #[arg(short, long, value_name = "N")]
    threads: Option<NonZeroUsize>,
}

/// A thread pool of `threads` workers, or one per core.
fn pool(threads: Option<NonZeroUsize>) -> Result<rayon::ThreadPool> {
    let mut builder = rayon::ThreadPoolBuilder::new();
    if let Some(threads) = threads {
        builder = builder.num_threads(threads.get());
    }
    builder.build().context("starting the worker threads")
}

/// Reads every read of `paths`.
fn read_all(paths: &[PathBuf]) -> Result<Reads> {
    let mut stream = ReadSource::new(paths).open().context("opening the reads")?;
    let mut reads = Reads::new();
    let mut record = ReadRecord::default();
    while stream
        .read_record(&mut record)
        .context("reading the reads")?
    {
        reads.push(&record.seq);
    }
    reads.shrink_to_fit();
    Ok(reads)
}

/// Polishes `contigs` with the short reads in `paths` and logs what changed.
fn polish_with(contigs: &mut [Contig], paths: &[PathBuf]) -> Result<()> {
    let reads = read_all(paths)?;
    log::info!(
        "bacteria: {} short reads, {} bases",
        reads.len(),
        reads.total_bases()
    );
    let stats = polish_short(contigs, &reads, &ShortParams::new());
    log::info!(
        "bacteria: short reads: {} aligned, {} of {} positions covered, {} changed",
        stats.aligned,
        stats.covered,
        stats.positions,
        stats.changed
    );
    Ok(())
}

/// Reads an assembly's contigs, with their headers (`circular=true` marks a circle, and
/// `telomeres=hairpin`, `hairpin=start` or `hairpin=end` a linear replicon's hairpin ends).
fn read_assembly(path: &Path) -> Result<(Vec<String>, Vec<Contig>)> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut headers = Vec::new();
    let mut contigs = Vec::new();
    for record in text.split('>').filter(|r| !r.trim().is_empty()) {
        let (header, seq) = record.split_once('\n').unwrap_or((record, ""));
        let seq: Vec<u8> = seq
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let has = |tag: &str| {
            header
                .split_whitespace()
                .any(|w| w.eq_ignore_ascii_case(tag))
        };
        let both = has("telomeres=hairpin");
        headers.push(header.trim().to_owned());
        contigs.push(Contig {
            seq,
            circular: has("circular=true"),
            hairpins: [both || has("hairpin=start"), both || has("hairpin=end")],
            depth: 0.0,
            check: Check::default(),
        });
    }
    Ok((headers, contigs))
}

/// Writes `contigs` under their `headers`, each `length=` updated.
fn write_assembly(path: &Path, headers: &[String], contigs: &[Contig]) -> Result<()> {
    let mut fasta =
        FastaWriter::create(path).with_context(|| format!("creating {}", path.display()))?;
    for (header, contig) in headers.iter().zip(contigs) {
        let header: Vec<String> = header
            .split_whitespace()
            .map(|w| {
                if w.starts_with("length=") {
                    format!("length={}", contig.seq.len())
                } else {
                    w.to_owned()
                }
            })
            .collect();
        fasta
            .write_record(header.join(" "), &contig.seq)
            .with_context(|| format!("writing {}", path.display()))?;
    }
    fasta
        .finish()
        .with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

pub(crate) fn run_polish_short(args: &PolishShortArgs) -> Result<()> {
    let pool = pool(args.threads)?;
    let (headers, mut contigs) = read_assembly(&args.assembly)?;
    pool.install(|| polish_with(&mut contigs, &args.reads))?;
    write_assembly(&args.output, &headers, &contigs)
}

pub(crate) fn run(args: &AssembleArgs) -> Result<()> {
    let pool = pool(args.threads)?;
    let reads = read_all(&args.reads)?;
    log::info!(
        "bacteria: {} reads, {} bases",
        reads.len(),
        reads.total_bases()
    );

    let mut params = BactParams::new(match args.kind {
        Kind::Ont => ReadKind::Ont,
        Kind::OntNoisy => ReadKind::OntNoisy,
        Kind::OntClean => ReadKind::Ont,
        Kind::Hifi => ReadKind::Hifi,
    });
    if args.kind == Kind::OntClean {
        params.noisy_divergence = None;
    }
    if let Some(rounds) = args.polish_rounds {
        params.polish_rounds = rounds;
    }
    args.dev.apply(&mut params);
    let mut assembly = pool
        .install(|| assemble(&reads, &params))
        .context("assembling")?;
    drop(reads);
    if !args.short_reads.is_empty() {
        pool.install(|| polish_with(&mut assembly.contigs, &args.short_reads))?;
    }

    if let Some(dir) = &args.dev.dump {
        dump(dir, &assembly)?;
    }
    write_contigs(&args.output, "contig", &assembly.contigs)?;
    if let Some(path) = &args.shallow {
        write_contigs(path, "shallow", &assembly.shallow)?;
    }
    if let Some(path) = &args.suspect {
        write_suspect(path, "contig", &assembly.contigs)?;
    }
    if let Some(path) = &args.gfa {
        write_gfa(path, "contig", &assembly.contigs)?;
    }
    for (i, contig) in assembly.contigs.iter().enumerate() {
        if !contig.check.is_clean(contig.is_whole()) {
            log::warn!(
                "bacteria: contig_{} ({} bases) needs a look: {} (min depth {:.1})",
                i + 1,
                contig.seq.len(),
                check_tag(contig),
                contig.check.min_depth
            );
        }
    }
    log::info!(
        "bacteria: {} contigs ({} circular), {} bases",
        assembly.contigs.len(),
        assembly.contigs.iter().filter(|c| c.circular).count(),
        assembly.contigs.iter().map(|c| c.seq.len()).sum::<usize>()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bin(start: u32, flag: Flag, depth: f64, stops: u64) -> FlaggedBin {
        FlaggedBin {
            start,
            end: start + picoasm::CHECK_BIN,
            flag,
            depth,
            stops,
        }
    }

    fn contig(len: usize, circular: bool, hairpins: [bool; 2], flagged: Vec<FlaggedBin>) -> Contig {
        Contig {
            seq: vec![b'A'; len],
            circular,
            hairpins,
            depth: 40.0,
            check: Check {
                min_depth: 2.0,
                median_depth: 40.0,
                flagged,
            },
        }
    }

    #[test]
    fn adjacent_flagged_bins_are_one_region_clipped_to_the_contig() {
        // Shallow over 2-4 kb, stops over 3-5 kb (most at 4 kb), and a shallow last bin past
        // the contig's 9,990 bases.
        let c = contig(
            9_990,
            true,
            [false; 2],
            vec![
                bin(2_000, Flag::Shallow, 6.0, 0),
                bin(3_000, Flag::Shallow, 2.5, 1),
                bin(3_000, Flag::Stop, 2.5, 7),
                bin(4_000, Flag::Stop, 30.0, 12),
                bin(9_000, Flag::Shallow, 9.0, 0),
            ],
        );
        let got: Vec<(u32, u32, String)> = suspect_regions(&c)
            .into_iter()
            .map(|s| (s.start, s.end, s.reason))
            .collect();
        assert_eq!(
            got,
            [
                (2_000, 4_000, "shallow:depth=2.5,median=40.0".to_owned()),
                (3_000, 5_000, "stops:reads=12,depth=30.0".to_owned()),
                (9_000, 9_990, "shallow:depth=9.0,median=40.0".to_owned()),
            ]
        );
    }

    #[test]
    fn a_contig_not_closed_is_suspect_over_its_length() {
        let regions = suspect_regions(&contig(5_000, false, [true, false], Vec::new()));
        assert_eq!(
            regions,
            [Suspect {
                start: 0,
                end: 5_000,
                reason: "linear".to_owned()
            }]
        );
        assert!(suspect_regions(&contig(5_000, true, [false; 2], Vec::new())).is_empty());
    }

    #[test]
    fn the_graph_links_circles_end_to_start_and_hairpins_to_their_reverse() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("asm.gfa");
        let contigs = [
            contig(8, true, [false; 2], Vec::new()),
            contig(6, false, [true, true], Vec::new()),
            contig(4, false, [false; 2], Vec::new()),
        ];
        write_gfa(&path, "contig", &contigs).expect("the graph is written");
        let text = std::fs::read_to_string(&path).expect("the graph is read back");
        assert_eq!(
            text,
            "H\tVN:Z:1.0\n\
             S\tcontig_1\tAAAAAAAA\tLN:i:8\tdp:f:40\n\
             S\tcontig_2\tAAAAAA\tLN:i:6\tdp:f:40\n\
             S\tcontig_3\tAAAA\tLN:i:4\tdp:f:40\n\
             L\tcontig_1\t+\tcontig_1\t+\t0M\n\
             L\tcontig_2\t+\tcontig_2\t-\t0M\n\
             L\tcontig_2\t-\tcontig_2\t+\t0M\n"
        );
    }

    #[test]
    fn a_bin_with_no_bases_left_has_no_region() {
        let mut gone = bin(1_000, Flag::Shallow, 1.0, 0);
        gone.end = gone.start;
        let c = contig(5_000, true, [false; 2], vec![gone]);
        assert!(suspect_regions(&c).is_empty());
    }

    #[test]
    fn a_line_with_hairpins_at_both_ends_is_whole_and_one_hairpin_is_not() {
        let both = contig(5_000, false, [true, true], Vec::new());
        assert!(suspect_regions(&both).is_empty());
        let start = contig(5_000, false, [true, false], Vec::new());
        assert_eq!(suspect_regions(&start)[0].reason, "linear");
    }

    #[test]
    fn the_bed_names_contigs_as_the_fasta_does_and_has_a_header() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("suspect.bed");
        let contigs = [
            contig(3_000, true, [false; 2], Vec::new()),
            contig(
                2_500,
                false,
                [false; 2],
                vec![bin(1_000, Flag::Stop, 20.0, 6)],
            ),
        ];
        write_suspect(&path, "contig", &contigs).expect("the BED is written");
        let text = std::fs::read_to_string(&path).expect("the BED is read back");
        assert_eq!(
            text,
            "#contig\tstart\tend\treason\n\
             contig_2\t0\t2500\tlinear\n\
             contig_2\t1000\t2000\tstops:reads=6,depth=20.0\n"
        );
        write_suspect(&path, "contig", &[]).expect("an empty BED is written");
        let text = std::fs::read_to_string(&path).expect("the BED is read back");
        assert_eq!(text, "#contig\tstart\tend\treason\n");
    }

    #[test]
    fn an_empty_assembly_is_a_graph_with_only_its_header() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("asm.gfa");
        write_gfa(&path, "contig", &[]).expect("the graph is written");
        let text = std::fs::read_to_string(&path).expect("the graph is read back");
        assert_eq!(text, "H\tVN:Z:1.0\n");
    }
}
