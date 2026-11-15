#!/usr/bin/env python3
"""Lists and classifies an assembly's small differences from its reference.

Usage: bact_errors.py reference.fa asm.fa [label]

Maps each contig with minimap2 (asm5, --cs) and prints one line per difference: label,
replicon, position, kind (sub/ins/del), the edit, the reference homopolymer run it sits in
(base and length), and the reference context (12 bases either side).
"""
import re
import subprocess
import sys

import bact_score


def runs_at(seq, pos):
    """Returns (base, length) of the homopolymer containing or just right of `pos`."""
    if pos >= len(seq):
        return ('-', 0)
    b = seq[pos]
    i = j = pos
    while i > 0 and seq[i - 1] == b:
        i -= 1
    while j < len(seq) and seq[j] == b:
        j += 1
    return (b, j - i)


def main():
    ref_path, asm_path = sys.argv[1], sys.argv[2]
    label = sys.argv[3] if len(sys.argv) > 3 else asm_path
    ref = {k.split()[0]: v.upper() for k, v in bact_score.read_fasta(ref_path).items()}
    out = subprocess.run(['minimap2', '-cx', 'asm5', '--cs', '-t', '4', ref_path, asm_path],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        f = line.split('\t')
        if int(f[11]) < 30 or int(f[10]) < 10000:
            continue
        cs = next(x for x in f if x.startswith('cs:Z:'))[5:]
        name, p = f[5], int(f[7])
        seq = ref[name]
        for op in re.findall(r'(:\d+|\*[a-z]{2}|[+-][a-z]+)', cs):
            if op[0] == ':':
                p += int(op[1:])
                continue
            kind = {'*': 'sub', '+': 'ins', '-': 'del'}[op[0]]
            edit = op[1:].upper()
            # An insertion or deletion of one repeated base sits in the run of that base.
            if kind == 'sub':
                b, n = runs_at(seq, p)
            else:
                b, n = runs_at(seq, p if kind == 'del' else p)
                if kind == 'ins' and (n == 0 or b != edit[0]):
                    b2, n2 = runs_at(seq, p - 1) if p > 0 else ('-', 0)
                    if b2 == edit[0]:
                        b, n = b2, n2
            ctx = seq[max(0, p - 12):p].lower() + seq[p:p + max(1, len(edit) if kind != 'ins' else 0)] + seq[p + max(1, len(edit) if kind != 'ins' else 0):p + 12 + max(1, len(edit))].lower()
            print(label, name, p, kind, edit, f'{b}{n}', ctx, sep='\t')
            if kind in ('sub',):
                p += 1
            elif kind == 'del':
                p += len(edit)


if __name__ == '__main__':
    main()
