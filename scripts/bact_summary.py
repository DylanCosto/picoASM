"""Summarises the bacterial benchmark per tool, as the Autocycler paper's tables do.

Usage: python3 scripts/bact_summary.py LABEL [LABEL ...]

Reads data/bacteria/bench/results.tsv (rows whose dataset is one of the paper's 30 read
subsets, `<genome>_<1-6>`, for each LABEL) and data/bacteria/bench/paper.tsv and
paper_medaka.tsv (the paper's own assemblies, scored by bact_score_paper.py), and prints for
each: assemblies, perfect ones, those with every replicon complete, median and total errors,
missing bases, extra contigs, misassemblies, and for our runs the median wall and CPU time and
the largest peak memory. Then the errors per genome (median over its six subsets).
"""
import collections
import csv
import pathlib
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
BENCH = ROOT / 'data' / 'bacteria' / 'bench'
GENOMES = {'sflex': 'Shigella_flexneri', 'kpneu': 'Klebsiella_pneumoniae',
           'ehorm': 'Enterobacter_hormaechei', 'linno': 'Listeria_innocua',
           'prett': 'Providencia_rettgeri'}


def ours(label):
    """The last row per subset dataset for `label`, as dicts like paper.tsv's."""
    rows = {}
    with open(BENCH / 'results.tsv') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            name = r['dataset']
            if r['label'] == label and name[:5] in GENOMES and name[-1].isdigit() and '_' in name:
                rows[name] = dict(r, genome=GENOMES[name[:5]], tool=label)
    return list(rows.values())


def paper(path):
    if not path.exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def summarise(name, rows):
    if not rows:
        return
    errors = [int(r['errors']) for r in rows]
    line = (f"{name:22} {len(rows):3} {sum(str(r['perfect']) == 'True' for r in rows):7} "
            f"{sum(int(r['complete']) == int(r['replicons']) for r in rows):12} "
            f"{statistics.median(errors):8g} {sum(errors):9} "
            f"{sum(int(r['missing_bases']) for r in rows):8} "
            f"{sum(int(r['extra_contigs']) for r in rows):6} "
            f"{sum(int(r['misassemblies']) for r in rows):6}")
    if 'wall_s' in rows[0]:
        line += (f"   {statistics.median(int(r['wall_s']) for r in rows):5g} s "
                 f"{statistics.median(int(r['cpu_s']) for r in rows):6g} CPU s "
                 f"{max(float(r['rss_gb']) for r in rows):5.2f} GB")
    print(line)


def main():
    tools = collections.defaultdict(list)
    for r in paper(BENCH / 'paper.tsv') + paper(BENCH / 'paper_medaka.tsv'):
        tools[r['tool']].append(r)
    for label in sys.argv[1:]:
        tools[label] = ours(label)
    print(f"{'tool':22} {'n':>3} {'perfect':>7} {'all complete':>12} {'median':>8} "
          f"{'total':>9} {'missing':>8} {'extra':>6} {'misasm':>6}")
    order = sorted(tools, key=lambda t: statistics.median(int(r['errors']) for r in tools[t])
                   if tools[t] else 1e9)
    for tool in order:
        summarise(tool, tools[tool])
    print()
    print('median errors per genome (over its six read subsets):')
    shown = [t for t in order if t in sys.argv[1:] or t in (
        'autocycler', 'autocycler_medaka', 'flye', 'flye_medaka', 'hybracter', 'canu')]
    print(f"{'genome':26}" + ''.join(f'{t[:18]:>19}' for t in shown))
    for genome in GENOMES.values():
        cells = []
        for t in shown:
            errs = [int(r['errors']) for r in tools[t] if r['genome'] == genome]
            cells.append(f'{statistics.median(errs):19g}' if errs else f"{'-':>19}")
        print(f'{genome:26}' + ''.join(cells))


if __name__ == '__main__':
    main()
