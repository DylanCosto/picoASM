#!/usr/bin/env python3
"""Counts how reads report homopolymer runs of known length.

Usage: bact_hp_truth.py reference.fa reads.fastq out.tsv [threads]

Maps the reads to the reference (minimap2 --cs=long) and, for every run of at least 3 of one
base in the reference with a different base either side, records the longest stretch of that
base in each read's piece between the flanking bases. Writes rows of: true length, run base in
the read's orientation, read length called, count.
"""
import collections
import re
import subprocess
import sys

import bact_score

COMP = str.maketrans('ACGT', 'TGCA')


def runs(seq, min_len=3):
    out = []
    i = 0
    while i < len(seq):
        j = i + 1
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len and i > 0 and j < len(seq):
            out.append((i, j - i))
        i = j
    return out


def main():
    ref_path, reads, out_path = sys.argv[1:4]
    threads = sys.argv[4] if len(sys.argv) > 4 else '8'
    refs = {k.split()[0]: v.upper() for k, v in bact_score.read_fasta(ref_path).items()}
    run_starts = {n: runs(s) for n, s in refs.items()}
    counts = collections.Counter()
    p = subprocess.Popen(['minimap2', '-cx', 'map-ont', '--cs=long', '--secondary=no', '-t', threads,
                          ref_path, reads], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    for line in p.stdout:
        f = line.split('\t')
        if int(f[11]) < 30:
            continue
        name, strand, ts, te = f[5], f[4], int(f[7]), int(f[8])
        cs = next(x for x in f if x.startswith('cs:Z:'))[5:].rstrip()
        # For each reference position, the read bases from it up to the next one.
        bases = [''] * (te - ts + 1)
        t = 0
        for op, s in re.findall(r'([=*+\-~])([A-Za-z]+)', cs):
            s = s.upper()
            if op == '=':
                for c in s:
                    bases[t] += c
                    t += 1
            elif op == '*':
                bases[t] += s[1]
                t += 1
            elif op == '+':
                if t > 0:
                    bases[t - 1] += s
            elif op == '-':
                t += len(s)
        seq = refs[name]
        for start, length in run_starts[name]:
            if start - 1 < ts or start + length + 1 > te:
                continue
            piece = ''.join(bases[start - 1 - ts:start + length + 1 - ts])
            b = seq[start]
            k = max((len(m) for m in re.findall(b + '+', piece)), default=0)
            ob = b if strand == '+' else b.translate(COMP)
            counts[(length, ob, k)] += 1
    with open(out_path, 'w') as out:
        for (length, b, k), n in sorted(counts.items()):
            print(length, b, k, n, sep='\t', file=out)


if __name__ == '__main__':
    main()
