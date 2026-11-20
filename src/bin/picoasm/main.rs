//! The picoasm command-line interface: `picoasm assemble` and `picoasm polish-short`.

mod assemble;

use std::process::ExitCode;

use anyhow::Result;
use clap::{Parser, Subcommand};

/// De novo assembler for bacterial isolates from Oxford Nanopore or PacBio HiFi reads.
///
/// picoasm overlaps the reads base by base, lays them out in a string graph and polishes the
/// contigs from all the reads. Inputs may be FASTA, FASTQ (plain, gzip or BGZF) or unaligned
/// BAM, in any mix; formats are detected from the file contents.
#[derive(Debug, Parser)]
#[command(name = "picoasm", version, propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Log more: -v for debug messages, -vv for trace messages (RUST_LOG overrides).
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    verbose: u8,

    /// Log only warnings and errors.
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Assemble a bacterial isolate from Oxford Nanopore or PacBio HiFi reads.
    ///
    /// Writes the contigs to a FASTA file; each header gives the contig's length, whether it
    /// is a whole circular replicon, and the depth of the reads on it.
    Assemble(Box<assemble::AssembleArgs>),
    /// Polish a bacterial assembly with short reads (Illumina) of the same isolate.
    ///
    /// For an assembly polished after `picoasm assemble` (with Dorado or Medaka): writes its
    /// contigs with each base the short reads clearly show otherwise changed.
    PolishShort(assemble::PolishShortArgs),
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let level = match (cli.quiet, cli.verbose) {
        (true, _) => log::LevelFilter::Warn,
        (false, 0) => log::LevelFilter::Info,
        (false, 1) => log::LevelFilter::Debug,
        (false, _) => log::LevelFilter::Trace,
    };
    env_logger::Builder::new()
        .filter_level(level)
        .parse_default_env()
        .format_timestamp_secs()
        .init();
    match run(cli.command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(command: Command) -> Result<()> {
    match command {
        Command::Assemble(args) => assemble::run(&args),
        Command::PolishShort(args) => assemble::run_polish_short(&args),
    }
}
