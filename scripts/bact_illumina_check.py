"""Usage: python3 scripts/bact_illumina_check.py reference.fa illumina_on_reference.bam asm.fa [asm.fa ...]

For each tool's differences from the reference (minimap2 --cs, chromosome-length alignments),
asks the Illumina reads mapped to the reference which side they take: the share of reads that
show the tool's version. >= 0.5: the reference is wrong there (the tool right); <= 0.1: the
tool's error; between: mixed (a variable site, or reads unsure)."""
import re, subprocess, sys
from collections import Counter
ref, bam, *asms = sys.argv[1:]
WINDOW = 15

def diffs(asm):
    paf = subprocess.run(['minimap2', '-cx', 'asm5', '--cs', '-t', '4', ref, asm],
                         capture_output=True, text=True, check=True).stdout
    out = set()
    for line in paf.splitlines():
        f = line.split('\t')
        if int(f[10]) < 100000: continue
        tpos = int(f[7]); cs = [x for x in f if x.startswith('cs:Z:')][0][5:]
        for op, val in re.findall(r'([:*+\-~])([0-9a-z]+)', cs):
            if op == ':': tpos += int(val)
            elif op == '*': out.add((f[5], tpos, 'sub', val[1].upper())); tpos += 1
            elif op == '+': out.add((f[5], tpos, 'ins', len(val)))
            elif op == '-': out.add((f[5], tpos, 'del', len(val))); tpos += len(val)
    return out

def pileup(chrom, a, b):
    out = subprocess.run(['samtools', 'mpileup', '-B', '-Q', '0', '-q', '0', '-d', '0',
                          '-f', ref, '-r', f'{chrom}:{a}-{b}', bam],
                         capture_output=True, text=True, check=True).stdout
    rows = {}
    for line in out.splitlines():
        f = line.split('\t')
        rows[int(f[1])] = (int(f[3]), f[4])
    return rows

def events(s):
    """(mismatch letters, indels as ('+'|'-', len)) of one pileup column."""
    letters, indels, i = Counter(), Counter(), 0
    while i < len(s):
        c = s[i]
        if c == '^': i += 2; continue
        if c in '+-':
            m = re.match(r'\d+', s[i + 1:]); n = int(m.group())
            indels[(c, n)] += 1; i += 1 + len(m.group()) + n; continue
        if c.upper() in 'ACGT': letters[c.upper()] += 1
        i += 1
    return letters, indels

def support(d):
    chrom, tpos, kind, x = d
    rows = pileup(chrom, max(1, tpos - WINDOW), tpos + WINDOW + 1)
    if kind == 'sub':
        depth, s = rows.get(tpos + 1, (0, ''))
        return (events(s)[0][x] / depth if depth else None), depth
    depth = rows.get(tpos, (0, ''))[0]
    sign = '+' if kind == 'ins' else '-'
    n = sum(events(s)[1][(sign, x)] for _, s in rows.values())
    return (min(1.0, n / depth) if depth else None), depth

cache = {}
for asm in asms:
    name = asm.split('/')[-2]
    verdicts = Counter()
    for d in sorted(diffs(asm)):
        if d not in cache: cache[d] = support(d)
        share, depth = cache[d]
        v = 'no reads' if share is None else 'ref wrong' if share >= 0.5 else 'tool error' if share <= 0.1 else 'mixed'
        verdicts[v] += 1
        if v in ('ref wrong', 'mixed'): verdicts[f'  {d[2]} {v}'] += 1
    print(f'{name:22} ' + ', '.join(f'{k}: {v}' for k, v in sorted(verdicts.items())))
print('reference-wrong sites:', sorted(d for d, (s, _) in cache.items() if s is not None and s >= 0.5))
