#!/usr/bin/env python3
"""Counts each tool's differences from a dataset's reference that no other tool shares.

Usage: bact_unshared.py DATASETS.tsv LABEL [LABEL ...]

For each dataset of the list, the assemblies data/bacteria/bench/<dataset>/<label>/asm.fa are compared with
the reference by bact_errors.py (substitutions and indels within aligned stretches). A difference another
tool's assembly has at the same place too (within 10 bases, same kind) may well be the reference's: a
reference built from the same reads, or from another culture. Each tool's count of differences no other
tool has (its own), and their bases (indels by length), is printed beside its total.
"""
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
BENCH = ROOT / 'data' / 'bacteria' / 'bench'
SLACK = 10


def differences(reference, assembly):
    """(replicon, position, kind, bases) of each difference, by bact_errors.py."""
    out = subprocess.run([sys.executable, str(ROOT / 'scripts' / 'bact_errors.py'), str(reference),
                          str(assembly), 'x'], capture_output=True, text=True, check=True).stdout
    found = []
    for line in out.splitlines():
        _, replicon, pos, kind, edit, *_ = line.split('\t')
        found.append((replicon, int(pos), kind, 1 if kind == 'sub' else len(edit)))
    return found


def in_any(d, others):
    """Whether any of the other tools' differences is `d` (same replicon and kind, within SLACK)."""
    replicon, pos, kind, _ = d
    return any(o[0] == replicon and o[2] == kind and abs(o[1] - pos) <= SLACK
               for other in others for o in other)


def main():
    rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
    labels = sys.argv[2:]
    print('dataset\t' + '\t'.join(f'{l} (own / all)' for l in labels))
    for row in rows:
        name, reference = row[0], ROOT / row[3]
        found = {}
        for label in labels:
            assembly = BENCH / name / label / 'asm.fa'
            if assembly.exists():
                found[label] = differences(reference, assembly)
        cells = []
        for label in labels:
            if label not in found:
                cells.append('-')
                continue
            others = [found[o] for o in found if o != label]
            own = [d for d in found[label] if not in_any(d, others)]
            cells.append(f'{len(own)} ev {sum(d[3] for d in own)} bp / {len(found[label])} ev')
        print(name + '\t' + '\t'.join(cells), flush=True)


if __name__ == '__main__':
    main()
