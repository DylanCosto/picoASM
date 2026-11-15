"""Scores the Autocycler paper's published assemblies with bact_score.py, for the baseline table.

Usage: python3 scripts/bact_score_paper.py ASSEMBLIES_DIR [--jobs N] [--out results.tsv]

ASSEMBLIES_DIR is the `assemblies/` folder of the figshare archive "Autocycler: assemblies and
reference genomes" (doi 10.26188/28910099): one folder per genome holding `reference.fasta` and
`<tool>_<subset>.fasta` for every tool and read subset (1-6). Writes one row per assembly to
data/bacteria/bench/paper.tsv (or --out), then prints, per tool, the medians and totals over all
genomes and subsets. These assemblies were made on the paper's hardware, so there is no time or
memory column; bact_bench.py measures those on this machine.
"""
import argparse
import csv
import multiprocessing
import pathlib
import re
import statistics
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bact_score  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
COLUMNS = ['genome', 'tool', 'subset', 'contigs', 'replicons', 'complete', 'chromosome_complete',
           'errors', 'missing_bases', 'extra_contigs', 'extra_bases', 'misassemblies', 'qscore',
           'perfect', 'per_replicon']
THREADS_PER_JOB = 2


def score_one(job):
    genome, tool, subset, fasta, reference = job
    with tempfile.TemporaryDirectory(prefix='bact_paper.') as work:
        m = bact_score.score(fasta, reference, THREADS_PER_JOB, work)
    per_replicon = ';'.join(f"{r['replicon']}:{r['status']}:{r['errors']}"
                            for r in m['replicon_rows'])
    return {'genome': genome, 'tool': tool, 'subset': subset, 'contigs': m['contigs'],
            'replicons': m['replicons'], 'complete': m['complete'],
            'chromosome_complete': m['chromosome_complete'], 'errors': m['errors'],
            'missing_bases': m['missing_bases'], 'extra_contigs': m['extra_contigs'],
            'extra_bases': m['extra_bases'], 'misassemblies': m['misassemblies'],
            'qscore': m['qscore'], 'perfect': m['perfect'], 'per_replicon': per_replicon}


def summary(rows):
    """Prints per tool: perfect, all replicons complete, median errors, total errors, etc."""
    by_tool = {}
    for r in rows:
        by_tool.setdefault(r['tool'], []).append(r)
    print('tool\tassemblies\tperfect\tall_complete\tchromosome_complete\tmedian_errors\t'
          'total_errors\ttotal_missing\textra_contigs\tmisassemblies')
    ranked = sorted(by_tool.items(), key=lambda kv: statistics.median(
        int(r['errors']) + int(r['missing_bases']) for r in kv[1]))
    for tool, rs in ranked:
        n = len(rs)
        print(f"{tool}\t{n}\t{sum(str(r['perfect']) == 'True' for r in rs)}\t"
              f"{sum(int(r['complete']) == int(r['replicons']) for r in rs)}\t"
              f"{sum(str(r['chromosome_complete']) == 'True' for r in rs)}\t"
              f"{statistics.median(int(r['errors']) for r in rs):g}\t"
              f"{sum(int(r['errors']) for r in rs)}\t"
              f"{sum(int(r['missing_bases']) for r in rs)}\t"
              f"{sum(int(r['extra_contigs']) for r in rs)}\t"
              f"{sum(int(r['misassemblies']) for r in rs)}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('assemblies', type=pathlib.Path)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--out', type=pathlib.Path,
                   default=ROOT / 'data' / 'bacteria' / 'bench' / 'paper.tsv')
    args = p.parse_args()
    jobs = []
    for genome in sorted(d for d in args.assemblies.iterdir() if d.is_dir()):
        reference = genome / 'reference.fasta'
        for fasta in sorted(genome.glob('*.fasta')):
            m = re.fullmatch(r'(.+)_(\d+)\.fasta', fasta.name)
            if m:
                jobs.append((genome.name, m.group(1), int(m.group(2)), fasta, reference))
    with multiprocessing.Pool(args.jobs) as pool:
        rows = pool.map(score_one, jobs, chunksize=1)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    summary(rows)


if __name__ == '__main__':
    main()
