"""Runs bacterial assemblers on the bacterial benchmark and scores them against the references.

Usage: python3 scripts/bact_bench.py --tool NAME[,NAME...] [--only DATASET,...] [--threads N]
                                     [--jobs N] [--label LABEL] [--reads ont|hifi] [--keep-work]

For each dataset in scripts/bact_datasets.tsv ($BACT_DATASETS, if set: one list or several,
comma-separated; a dataset named twice runs once) with reads of the chosen kind (paths relative
to the repository root), runs each tool, `--jobs` at a time with `--threads` each, under `nice -n 19` and `/usr/bin/time -v` (wall time, and
the peak memory of its largest process), scores its contigs with `bact_score.py`, and appends
one row to data/bacteria/bench/results.tsv. Outputs go to
data/bacteria/bench/<dataset>/<label or tool>/: the assembly as asm.fa, the score as
score.json, the tool's log; the tool's own work directory is deleted unless --keep-work is
given.

Tools (each run as a shell command; the assemblers in the micromamba environment $ASM_ENV,
default `autocycler`, Medaka in $MEDAKA_ENV, default `medaka`; see the README for the installs):
- `flye`: Flye, `--nano-hq` or `--pacbio-hifi`.
- `flye_medaka`: Flye `--nano-hq` and then Medaka's bacterial model, the core of ONT's EPI2ME
  `wf-bacterial-genomes`; ONT only.
- `autocycler`: Autocycler's fully automated pipeline: four read subsets, each assembled by
  every assembler in $AUTOCYCLER_ASSEMBLERS (default the official script's nine: raven myloasm
  miniasm flye metamdbg necat nextdenovo plassembler canu), then compress, cluster, trim,
  resolve and combine.
- `autocycler_medaka`: the same, then Medaka; ONT only.
- `raven`, `miniasm` (miniasm and Minipolish), through `autocycler helper`; `hifiasm` from
  PATH (HiFi only, primary contigs).
- `picoasm`: `picoasm assemble` from target/release (or $PICOASM), with any extra settings in
  $PICOASM_ARGS.
- `picoasm_medaka`: the same, then Medaka; ONT only.
- `picoasm_dorado`: the same, then `dorado polish --bacteria` (on the GPU) with the reads
  aligned by `dorado aligner`; Dorado from $DORADO (default `dorado` on PATH); ONT only.
- `picoasm_short`, `picoasm_dorado_short`, `flye_medaka_short`: `picoasm`,
  `picoasm_dorado` and `flye_medaka` with the contigs then polished by the isolate's Illumina
  reads (`picoasm polish-short`): the files
  `illumina/*_1.fastq.gz` and `*_2.fastq.gz` beside the reference; datasets without them are
  skipped.
- `picoasm_multi`: Autocycler's pipeline with picoasm as its only assembler, on its four
  read subsets and on all the reads; `picoasm_multi_dorado`: the same, then Dorado.

A dataset list's optional fifth column holds comma-separated flags:
- `structure-only`: the reference no longer matches the isolate base for base (an old finished
  genome the culture has drifted from), so the row's errors and Q are written as NA;
  completeness, missing bases, extra contigs and misassemblies still count, and score.json keeps
  the error count.
"""
import argparse
import concurrent.futures
import csv
import datetime
import json
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bact_score  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATASETS = [pathlib.Path(p) for p in
            (os.environ.get('BACT_DATASETS') or str(ROOT / 'scripts' / 'bact_datasets.tsv')).split(',')]
BENCH = ROOT / 'data' / 'bacteria' / 'bench'
RESULTS = BENCH / 'results.tsv'
COLUMNS = ['date', 'commit', 'label', 'tool', 'dataset', 'reads', 'threads', 'wall_s', 'cpu_s',
           'rss_gb', 'contigs', 'total_mb', 'replicons', 'complete', 'chromosome_complete',
           'errors', 'missing_bases', 'extra_contigs', 'extra_bases', 'misassemblies', 'qscore',
           'perfect', 'per_replicon']

AUTOCYCLER_ASSEMBLERS = 'raven myloasm miniasm flye metamdbg necat nextdenovo plassembler canu'


def datasets():
    """Rows of the dataset lists: name, ONT reads, HiFi reads, reference, flags (each name once)."""
    rows = {}
    for listed in DATASETS:
        with open(listed) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    r = line.rstrip('\n').split('\t')
                    rows.setdefault(r[0], r)
    path = lambda p: ROOT / p if p else None
    return [{'name': r[0], 'ont': path(r[1]), 'hifi': path(r[2]), 'reference': ROOT / r[3],
             'flags': set(r[4].split(',')) if len(r) > 4 and r[4] else set()}
            for r in rows.values()]


def in_env(env_var, default, command):
    """`command` run in the micromamba environment named by $env_var (else `default`)."""
    env = os.environ.get(env_var, default)
    return f'micromamba run -n {shlex.quote(env)} {command}'


def asm(command):
    return in_env('ASM_ENV', 'autocycler', command)


# The model `--bacteria` picks for R10.4.1 sup or hac reads (v4.2 to v5.2), named outright:
# reads from the SRA lack the basecaller's model in their headers, and `--bacteria` refuses them.
MEDAKA_MODEL = 'r1041_e82_400bps_bacterial_methylation'


def medaka(reads, draft, out, threads):
    q = lambda p: shlex.quote(str(p))
    return in_env('MEDAKA_ENV', 'medaka', f'medaka_consensus -i {q(reads)} -d {q(draft)} '
                                          f'-o {q(out)} -t {threads} -m {MEDAKA_MODEL}')


# The public read sets carry no basecaller model in a BAM header for Dorado to pick its model by,
# so the bacterial one (the same as Medaka's) is named: $DORADO_MODEL, a folder `dorado download
# --model dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0` wrote.
DORADO_MODEL = os.environ.get(
    'DORADO_MODEL',
    'dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0')


def dorado_polish(reads, draft, out, threads):
    """Shell lines that polish `draft` with `dorado polish --bacteria` into `out`/polished.fasta."""
    q = lambda p: shlex.quote(str(p))
    dorado = q(os.environ.get('DORADO', 'dorado'))
    bam = out / 'aligned.bam'
    return (f'mkdir -p {q(out)}\n'
            f'{dorado} aligner -t {threads} {q(draft)} {q(reads)} | samtools sort -@ 4 -o {q(bam)} -\n'
            f'samtools index {q(bam)}\n'
            f'{dorado} polish --bacteria -t {threads} --model-override {q(DORADO_MODEL)} '
            f'{q(bam)} {q(draft)} > {q(out / "polished.fasta")}')


def autocycler_script(kind, reads, work, threads, genome_size, polish, picoasm=None):
    """The shell commands of Autocycler's fully automated pipeline; with `picoasm` (its path),
    picoasm assembles each read subset and all the reads instead of Autocycler's assemblers."""
    ac = asm('autocycler')
    assemblers = os.environ.get('AUTOCYCLER_ASSEMBLERS', AUTOCYCLER_ASSEMBLERS).split()
    read_type = 'ont_r10' if kind == 'ont' else 'pacbio_hifi'
    lines = ['set -euo pipefail', f'cd {shlex.quote(str(work))}',
             f'{ac} subsample --reads {shlex.quote(str(reads))} --out_dir subsampled_reads '
             f'--genome_size {genome_size}',
             'mkdir -p assemblies']
    if picoasm:
        samples = [f'subsampled_reads/sample_{i}.fastq' for i in ('01', '02', '03', '04')]
        for i, sample in enumerate(samples + [shlex.quote(str(reads))]):
            lines.append(f'{shlex.quote(picoasm)} assemble --reads {kind} -t {threads} '
                         f'-o assemblies/picoasm_{i}.fasta {sample} 2> assemblies/picoasm_{i}.log')
        assemblers = []
    for assembler in assemblers:
        for i in ('01', '02', '03', '04'):
            lines.append(
                f'{ac} helper {assembler} --reads subsampled_reads/sample_{i}.fastq '
                f'--out_prefix assemblies/{assembler}_{i} --threads {threads} '
                f'--genome_size {genome_size} --read_type {read_type} || true')
    lines += ['rm -f assemblies/*.gfa',
              f'{ac} compress -i assemblies -a autocycler_out --threads {threads}',
              f'{ac} cluster -a autocycler_out',
              'for c in autocycler_out/clustering/qc_pass/cluster_*; do',
              f'  {ac} trim -c "$c" --threads {threads}',
              f'  {ac} resolve -c "$c"',
              'done',
              f'{ac} combine -a autocycler_out '
              '-i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa']
    out = 'autocycler_out/consensus_assembly.fasta'
    if polish == 'dorado':
        lines.append(dorado_polish(reads, out, pathlib.Path('dorado'), threads))
        out = 'dorado/polished.fasta'
    elif polish:
        lines.append(medaka(reads, out, 'medaka', threads))
        out = 'medaka/consensus.fasta'
    return '\n'.join(lines), work / out


SHORT_TOOLS = {'picoasm_short': 'picoasm', 'picoasm_dorado_short': 'picoasm_dorado',
               'flye_medaka_short': 'flye_medaka'}


def short_reads(reference):
    """The Illumina read files beside `reference` (in illumina/), or [] if there are none."""
    return sorted((reference.parent / 'illumina').glob('*_[12].fastq.gz'))


def tool_script(tool, kind, reads, work, threads, genome_size, short=()):
    """Returns (shell script, path of the assembly it writes) for one tool."""
    q = lambda p: shlex.quote(str(p))
    picoasm = os.environ.get('PICOASM') or str(ROOT / 'target' / 'release' / 'picoasm')
    if tool in SHORT_TOOLS:
        script, assembly = tool_script(SHORT_TOOLS[tool], kind, reads, work, threads, genome_size)
        polished = work / 'short_polished.fasta'
        files = ' '.join(q(f) for f in short)
        return (f'set -euo pipefail\n{script}\n'
                f'{q(picoasm)} polish-short -t {threads} -o {q(polished)} {q(assembly)} {files}',
                polished)
    flye_mode = '--nano-hq' if kind == 'ont' else '--pacbio-hifi'
    flye = asm(f'flye {flye_mode} {q(reads)} --out-dir {q(work / "flye")} --threads {threads}')
    if tool == 'flye':
        return flye, work / 'flye' / 'assembly.fasta'
    if tool == 'flye_medaka':
        polish = medaka(reads, work / 'flye' / 'assembly.fasta', work / 'medaka', threads)
        return f'set -euo pipefail\n{flye}\n{polish}', work / 'medaka' / 'consensus.fasta'
    if tool in ('autocycler', 'autocycler_medaka'):
        return autocycler_script(kind, reads, work, threads, genome_size,
                                 tool == 'autocycler_medaka')
    if tool in ('raven', 'miniasm'):
        read_type = 'ont_r10' if kind == 'ont' else 'pacbio_hifi'
        helper = asm(f'autocycler helper {tool} --reads {q(reads)} --out_prefix {q(work / tool)} '
                     f'--threads {threads} --genome_size {genome_size} --read_type {read_type}')
        return helper, work / f'{tool}.fasta'
    if tool == 'hifiasm':
        # hifiasm from PATH, not the assemblers' environment.
        prefix = work / 'hifiasm'
        gfa_to_fasta = f"awk '/^S/{{print \">\"$2\"\\n\"$3}}' {q(prefix)}.bp.p_ctg.gfa"
        return (f'set -euo pipefail\nhifiasm -t {threads} -o {q(prefix)} {q(reads)}\n'
                f'{gfa_to_fasta} > {q(work / "p_ctg.fa")}', work / 'p_ctg.fa')
    if tool in ('picoasm_multi', 'picoasm_multi_dorado'):
        polish = 'dorado' if tool == 'picoasm_multi_dorado' else None
        return autocycler_script(kind, reads, work, threads, genome_size, polish, picoasm)
    if tool in ('picoasm', 'picoasm_medaka', 'picoasm_dorado'):
        extra = os.environ.get('PICOASM_ARGS', '')
        bacteria = (f'{q(picoasm)} assemble --reads {kind} -t {threads} -o {q(work / "asm.fa")} '
                    f'{extra} {q(reads)}')
        if tool == 'picoasm':
            return bacteria, work / 'asm.fa'
        if tool == 'picoasm_dorado':
            polish = dorado_polish(reads, work / 'asm.fa', work / 'dorado', threads)
            return f'set -euo pipefail\n{bacteria}\n{polish}', work / 'dorado' / 'polished.fasta'
        polish = medaka(reads, work / 'asm.fa', work / 'medaka', threads)
        return f'set -euo pipefail\n{bacteria}\n{polish}', work / 'medaka' / 'consensus.fasta'
    sys.exit(f'unknown tool {tool}')


def time_used(stderr):
    """Returns wall seconds, CPU seconds and peak GB that `/usr/bin/time -v` wrote."""
    log = stderr.read_text()
    wall = re.findall(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', log)[-1]
    seconds = sum(float(x) * 60 ** i for i, x in enumerate(reversed(wall.split(':'))))
    user = float(re.findall(r'User time \(seconds\): (\S+)', log)[-1])
    system = float(re.findall(r'System time \(seconds\): (\S+)', log)[-1])
    rss_kb = int(re.findall(r'Maximum resident set size \(kbytes\): (\d+)', log)[-1])
    return seconds, user + system, rss_kb / 1e6


def git_commit():
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT,
                            capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.strip()
    return commit + ('+dirty' if dirty else '')


def run(tool, kind, dataset, threads, label, keep_work):
    reads = dataset[kind]
    genome_size = sum(len(s) for s in bact_score.read_fasta(dataset['reference']).values())
    out = BENCH / dataset['name'] / label
    work = out / 'work'
    if out.exists():
        shutil.rmtree(out)
    work.mkdir(parents=True)
    script, assembly = tool_script(tool, kind, reads, work, threads, genome_size,
                                   short_reads(dataset['reference']))
    (out / 'run.sh').write_text(script + '\n')
    command = ['/usr/bin/time', '-v', 'nice', '-n', '19', 'bash', str(out / 'run.sh')]
    with open(out / 'log.txt', 'w') as log, open(out / 'time.txt', 'w') as err:
        # `time -v` writes to its own stderr after the tool's; keep the two apart.
        status = subprocess.run(command, stdout=log, stderr=err, cwd=work).returncode
    wall, cpu, rss = time_used(out / 'time.txt')
    if status != 0 or not assembly.exists():
        print(f'{dataset["name"]}\t{label}\tfailed (exit {status}); see {out}/time.txt',
              file=sys.stderr)
        return None
    shutil.copy(assembly, out / 'asm.fa')
    m = bact_score.score(out / 'asm.fa', dataset['reference'], threads, out / 'score')
    shutil.rmtree(out / 'score', ignore_errors=True)
    (out / 'score.json').write_text(json.dumps(m, indent=1))
    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)
    structure_only = 'structure-only' in dataset['flags']
    per_replicon = ';'.join(f"{r['replicon']}:{r['length']}:{r['status']}:"
                            f"{'NA' if structure_only else r['errors']}"
                            for r in m['replicon_rows'])
    if structure_only:
        m = dict(m, errors='NA', qscore='NA')
    return {'date': datetime.date.today().isoformat(), 'commit': git_commit(), 'label': label,
            'tool': tool, 'dataset': dataset['name'], 'reads': kind, 'threads': threads,
            'wall_s': round(wall), 'cpu_s': round(cpu), 'rss_gb': round(rss, 2),
            'contigs': m['contigs'], 'total_mb': round(m['total_bases'] / 1e6, 3),
            'replicons': m['replicons'], 'complete': m['complete'],
            'chromosome_complete': m['chromosome_complete'], 'errors': m['errors'],
            'missing_bases': m['missing_bases'], 'extra_contigs': m['extra_contigs'],
            'extra_bases': m['extra_bases'], 'misassemblies': m['misassemblies'],
            'qscore': m['qscore'], 'perfect': m['perfect'], 'per_replicon': per_replicon}


def append(row):
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    new = not RESULTS.exists()
    with open(RESULTS, 'a', newline='') as f:
        writer = csv.DictWriter(f, COLUMNS, delimiter='\t')
        if new:
            writer.writeheader()
        writer.writerow(row)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--tool', required=True, help='comma-separated tool names')
    p.add_argument('--only', help='comma-separated dataset names')
    p.add_argument('--reads', choices=['ont', 'hifi'], default='ont')
    p.add_argument('--threads', type=int, default=16, help='threads per run')
    p.add_argument('--jobs', type=int, default=1, help='runs at once')
    p.add_argument('--label', help='output label [default: TOOL-READS]')
    p.add_argument('--keep-work', action='store_true')
    args = p.parse_args()
    only = set(args.only.split(',')) if args.only else None
    runs = []
    for tool in args.tool.split(','):
        if args.reads == 'hifi' and (tool.endswith('medaka') or tool.endswith('dorado')):
            print(f'{tool} is for ONT reads only; skipped', file=sys.stderr)
            continue
        for d in datasets():
            if (only and d['name'] not in only) or d[args.reads] is None:
                continue
            if tool in SHORT_TOOLS and not short_reads(d['reference']):
                continue
            runs.append((tool, d, args.label or f'{tool}-{args.reads}'))
    # Runs go on in threads (each waits on its tool's process); rows are written here only.
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(args.jobs, 1)) as pool:
        futures = {pool.submit(run, tool, args.reads, d, args.threads, label, args.keep_work):
                   (d, label) for tool, d, label in runs}
        for future in concurrent.futures.as_completed(futures):
            d, label = futures[future]
            row = future.result()
            if row is None:
                continue
            append(row)
            print(f"{d['name']}\t{label}\t{row['wall_s']} s\t{row['rss_gb']} GB\t"
                  f"complete {row['complete']}/{row['replicons']}\terrors {row['errors']}\t"
                  f"missing {row['missing_bases']}\textra {row['extra_contigs']}\t"
                  f"misassemblies {row['misassemblies']}\tQ {row['qscore']}", flush=True)


if __name__ == '__main__':
    main()
