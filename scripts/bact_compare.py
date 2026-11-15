#!/usr/bin/env python3
"""Tabulates each dataset of a list against every tool run on it (the last row per dataset and label in
data/bacteria/bench/results.tsv): complete replicons, errors, missing bases, misassemblies, extra contigs and
wall time, picoasm first.

Usage: bact_compare.py DATASETS.tsv PICOASM_LABEL [OTHER_LABEL ...]
"""
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent


def main():
    datasets = [l.split('\t')[0] for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
    labels = sys.argv[2:]
    rows = {}
    for r in csv.DictReader(open(ROOT / 'data/bacteria/bench/results.tsv'), delimiter='\t'):
        rows[(r['dataset'], r['label'])] = r
    print('dataset\t' + '\t'.join(labels))
    for d in datasets:
        cells = []
        for label in labels:
            r = rows.get((d, label))
            cells.append('-' if r is None else
                         f"{r['complete']}/{r['replicons']} e{r['errors']} m{r['missing_bases']} "
                         f"x{r['misassemblies']} +{r['extra_contigs']} {r['wall_s']}s")
        print(d + '\t' + '\t'.join(cells))


if __name__ == '__main__':
    main()
