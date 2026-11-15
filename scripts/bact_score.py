"""Scores a bacterial assembly against a trusted complete reference, replicon by replicon.

Usage: python3 scripts/bact_score.py assembly.fa reference.fa [--threads N] [--json out.json]
                                     [--work DIR] [--keep]

Needs minimap2 on PATH (or $MINIMAP2). Both FASTA files may be gzipped.

Each contig of at least 500 bp (MIN_CONTIG) is aligned to the reference
(`minimap2 -cx asm5 --eqx`) and assigned to the replicon that takes most of its aligned bases;
contigs that align nowhere, and those beside the one that holds a complete replicon, are
"extra". A replicon is **complete** when one contig holds it whole: the contig that takes most
of the replicon's aligned bases is within 1% of its length, and once both are rotated to a
shared base and aligned end to end the alignment covers 99% of each. The errors of a complete
replicon are those of that end-to-end alignment, so they include the ones where the assembler
closed the circle: mismatches, inserted and deleted bases, and bases of either sequence left
out of the alignment. A contig that holds a small replicon two (or more) times over, as
assemblers make of small plasmids, is "doubled", not complete.

For a replicon that is not complete, the errors are those inside its contigs' alignments, and
the reference bases no contig covers are "missing". A **misassembly** is a pair of consecutive
alignments (each at least 1 kb, MIN_BREAK_ALIGNMENT) of one contig that land on different
replicons or strands, or that jump more than 1 kb (MAX_JUMP) on the replicon, other than the
step from a circular replicon's end to its start.

A `mask.bed` beside the reference (replicon, 0-based start, end) lists windows the reference
cannot call, such as long homopolymers or sites mixed in the culture: differences inside them are
not counted as errors (they are reported as "masked"). References without one are scored whole.

Prints a table per replicon and a one-line summary; --json writes all of it. The genome's Q is
-10 log10((errors + missing bases) / reference length). An assembly is **perfect** when every
replicon is complete with no error, and there is no extra contig and no misassembly.
"""
import argparse
import bisect
import collections
import gzip
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

MIN_CONTIG = 500
MIN_ASSIGN_ALIGNMENT = 200
MIN_BREAK_ALIGNMENT = 1_000
MAX_JUMP = 1_000
COMPLETE_SHARE = 0.99
LENGTH_TOLERANCE = 0.01
WRAP_SLACK = 1_000
ANCHOR_RUN = 50

Alignment = collections.namedtuple(
    'Alignment', 'query qlen qstart qend strand target tlen tstart tend matches block nm cigar')


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_fasta(path):
    """Returns {name: sequence} (upper case) in file order."""
    seqs = {}
    name, parts = None, []
    with open_text(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(parts).upper()
    return seqs


def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for name, seq in seqs.items():
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + '\n')


_COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def minimap2(target, query, out, threads):
    exe = os.environ.get('MINIMAP2') or shutil.which('minimap2')
    if not exe:
        sys.exit('minimap2 not found: put it on PATH or set $MINIMAP2')
    with open(out, 'w') as f:
        subprocess.run([exe, '-cx', 'asm5', '--eqx', '-t', str(threads), str(target), str(query)],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)


def read_paf(path):
    out = []
    with open(path) as f:
        for line in f:
            r = line.rstrip('\n').split('\t')
            tags = {t[:2]: t[5:] for t in r[12:]}
            if tags.get('tp', 'P') != 'P':
                continue
            out.append(Alignment(r[0], int(r[1]), int(r[2]), int(r[3]), r[4], r[5], int(r[6]),
                                 int(r[7]), int(r[8]), int(r[9]), int(r[10]),
                                 int(tags.get('NM', 0)), tags.get('cg', '')))
    return out


class Mask:
    """The windows of one replicon left out of scoring, as sorted (start, end) pairs."""

    def __init__(self, windows=()):
        self.windows = sorted(windows)
        self.starts = [s for s, _ in self.windows]

    def __bool__(self):
        return bool(self.windows)

    def __contains__(self, pos):
        i = bisect.bisect_right(self.starts, pos) - 1
        return i >= 0 and pos < self.windows[i][1]


def read_mask(reference):
    """{replicon: Mask} from the `mask.bed` beside `reference`, or {} if there is none."""
    path = pathlib.Path(reference).parent / 'mask.bed'
    windows = collections.defaultdict(list)
    if path.exists():
        with open(path) as f:
            for line in f:
                r = line.split('\t')
                if len(r) >= 3 and not line.startswith('#'):
                    windows[r[0]].append((int(r[1]), int(r[2])))
    return {name: Mask(w) for name, w in windows.items()}


def cigar_errors(cigar, tstart=0, mask=None, shift=0, length=0):
    """Returns (mismatches, inserted bases, deleted bases, masked) of an `--eqx` CIGAR whose
    target starts at `tstart`. With a mask, differences at target positions it holds are counted
    as masked instead; target position p is replicon position (p + shift) % length (the target
    may be the replicon rotated by `shift`)."""
    if not mask:
        counts = collections.Counter()
        for n, op in re.findall(r'(\d+)([=XID])', cigar):
            counts[op] += int(n)
        return counts['X'], counts['I'], counts['D'], 0
    held = lambda p: (p + shift) % length in mask
    x = i = d = masked = 0
    t = tstart
    for n, op in re.findall(r'(\d+)([=XID])', cigar):
        n = int(n)
        if op == 'X':
            m = sum(held(t + k) for k in range(n))
            x, masked = x + n - m, masked + m
            t += n
        elif op == 'D':
            m = sum(held(t + k) for k in range(n))
            d, masked = d + n - m, masked + m
            t += n
        elif op == 'I':
            # An insertion sits between t - 1 and t: masked if either side is.
            if held(t) or held(t - 1):
                masked += n
            else:
                i += n
        else:
            t += n
    return x, i, d, masked


def covered(intervals):
    total, start, end = 0, -1, -1
    for a, b in sorted(intervals):
        if a > end:
            total += max(0, end - start)
            start, end = a, b
        else:
            end = max(end, b)
    return total + max(0, end - start)


def n50(lengths):
    total, seen = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        seen += x
        if seen * 2 >= total:
            return x
    return 0


def anchor(alignment):
    """Returns (query pos, target pos) at the middle of the longest exact match run of a '+'
    alignment, so both sequences can be rotated to start at the same base."""
    q, t = alignment.qstart, alignment.tstart
    best = (-1, q, t)
    for n, op in re.findall(r'(\d+)([=XID])', alignment.cigar):
        n = int(n)
        if op == '=':
            if n > best[0]:
                best = (n, q + n // 2, t + n // 2)
            q, t = q + n, t + n
        elif op == 'X':
            q, t = q + n, t + n
        elif op == 'I':
            q += n
        else:
            t += n
    return best[1], best[2]


def misassemblies(alignments, replicon_len):
    """Counts the breaks between consecutive alignments of each contig (see the module doc)."""
    by_contig = collections.defaultdict(list)
    for a in alignments:
        if a.block >= MIN_BREAK_ALIGNMENT:
            by_contig[a.query].append(a)
    breaks = []
    for contig, blocks in by_contig.items():
        blocks.sort(key=lambda a: a.qstart)
        for a, b in zip(blocks, blocks[1:]):
            if a.target != b.target or a.strand != b.strand:
                breaks.append((contig, 'replicon' if a.target != b.target else 'strand'))
                continue
            dq = b.qstart - a.qend
            dt = (b.tstart - a.tend) if a.strand == '+' else (a.tstart - b.tend)
            if abs(dt - dq) <= MAX_JUMP:
                continue
            length = replicon_len[a.target]
            # Across a circular replicon's origin: the step wraps from one end to the other.
            if a.strand == '+':
                wrapped = a.tend > length - WRAP_SLACK and b.tstart < WRAP_SLACK
            else:
                wrapped = a.tstart < WRAP_SLACK and b.tend > length - WRAP_SLACK
            if not wrapped:
                breaks.append((contig, f'jump {abs(dt - dq)} bp'))
    return breaks


def end_to_end(contig_seq, replicon_seq, alignment, work, threads, tag, mask=None):
    """Rotates an oriented contig and its replicon to one shared base and aligns them end to
    end; returns (errors, share of each covered, details). Differences in `mask` (replicon
    coordinates) are counted as masked, not as errors."""
    q, t = anchor(alignment)
    contig_rot = contig_seq[q:] + contig_seq[:q]
    replicon_rot = replicon_seq[t:] + replicon_seq[:t]
    qa, ta = work / f'{tag}.contig.fa', work / f'{tag}.replicon.fa'
    write_fasta(qa, {'contig': contig_rot})
    write_fasta(ta, {'replicon': replicon_rot})
    paf = work / f'{tag}.paf'
    minimap2(ta, qa, paf, threads)
    blocks = sorted((a for a in read_paf(paf) if a.strand == '+'), key=lambda a: a.qstart)
    mism = ins = dels = masked = 0
    qcov = covered([(a.qstart, a.qend) for a in blocks])
    tcov = covered([(a.tstart, a.tend) for a in blocks])
    for a in blocks:
        x, i, d, m = cigar_errors(a.cigar, a.tstart, mask, t, len(replicon_seq))
        mism, ins, dels, masked = mism + x, ins + i, dels + d, masked + m
    # Bases of either sequence outside every alignment are errors too (a missing or doubled
    # piece at the circle's junction, or a large indel that split the alignment).
    unaligned = (len(contig_rot) - qcov) + (len(replicon_rot) - tcov)
    errors = mism + ins + dels + unaligned
    details = {'mismatches': mism, 'insertions': ins, 'deletions': dels,
               'unaligned_bases': unaligned, 'pieces': len(blocks), 'masked': masked}
    return errors, qcov / len(contig_rot), tcov / len(replicon_rot), details


def score(assembly, reference, threads=4, work=None):
    """Returns the metrics dict of one assembly against one reference."""
    ref = read_fasta(reference)
    masks = read_mask(reference)
    asm = {n: s for n, s in read_fasta(assembly).items() if len(s) >= MIN_CONTIG}
    replicon_len = {n: len(s) for n, s in ref.items()}
    genome = sum(replicon_len.values())
    work = pathlib.Path(work)
    work.mkdir(parents=True, exist_ok=True)

    raw = work / 'contigs.fa'
    write_fasta(raw, asm)
    paf = work / 'contigs.paf'
    minimap2(reference, raw, paf, threads)
    alignments = [a for a in read_paf(paf) if a.block >= MIN_ASSIGN_ALIGNMENT]
    breaks = misassemblies(alignments, replicon_len)

    # Assign each contig to the replicon (and strand) with most of its aligned bases.
    votes = collections.defaultdict(collections.Counter)
    for a in alignments:
        votes[a.query][(a.target, a.strand)] += a.matches
    assigned = collections.defaultdict(list)
    for contig, v in votes.items():
        (replicon, strand), _ = v.most_common(1)[0]
        assigned[replicon].append((contig, strand))
    extra = [c for c in asm if c not in votes]

    replicons = []
    total_errors = total_missing = total_masked = 0
    for name, seq in ref.items():
        mine = assigned.get(name, [])
        mask = masks.get(name)
        row = {'replicon': name, 'length': len(seq), 'contigs': [c for c, _ in mine],
               'status': 'missing', 'errors': 0, 'missing_bases': len(seq), 'masked': 0}
        if mine:
            own = [a for a in alignments if a.target == name and
                   any(a.query == c for c, _ in mine)]
            row['missing_bases'] = len(seq) - covered([(a.tstart, a.tend) for a in own])
            by_contig = collections.Counter()
            for a in own:
                by_contig[a.query] += a.matches
            top, top_strand = max(mine, key=lambda cs: by_contig[cs[0]])
            ratio = len(asm[top]) / len(seq)
            row['status'] = (f'fragmented ({len(mine)} contigs)' if len(mine) > 1
                             else f'wrong length (x{ratio:.3f})')
            if ratio >= 2 - 2 * LENGTH_TOLERANCE and by_contig[top] >= COMPLETE_SHARE * len(seq):
                row['status'] = f'doubled (x{ratio:.2f})'
            if abs(ratio - 1) <= LENGTH_TOLERANCE:
                contig_seq = asm[top] if top_strand == '+' else revcomp(asm[top])
                oriented = work / f'{len(replicons)}.oriented.fa'
                write_fasta(oriented, {top: contig_seq})
                opaf = work / f'{len(replicons)}.oriented.paf'
                minimap2(reference, oriented, opaf, threads)
                best = max((a for a in read_paf(opaf) if a.target == name and a.strand == '+'),
                           key=lambda a: a.matches, default=None)
                if best is not None:
                    errors, qshare, tshare, details = end_to_end(
                        contig_seq, seq, best, work, threads, str(len(replicons)), mask)
                    if qshare >= COMPLETE_SHARE and tshare >= COMPLETE_SHARE:
                        row.update(status='complete', errors=errors, missing_bases=0,
                                   contig=top, contig_length=len(contig_seq), **details)
            if row['status'] != 'complete':
                # Errors inside the alignments of the replicon's contigs only.
                errs = masked = 0
                for a in own:
                    x, i, d, m = cigar_errors(a.cigar, a.tstart, mask, 0, len(seq))
                    errs, masked = errs + x + i + d, masked + m
                row.update(errors=errs, masked=masked)
        total_errors += row['errors']
        total_missing += row['missing_bases']
        total_masked += row['masked']
        replicons.append(row)

    # Contigs beyond the one a complete replicon needs (duplicates, leftover fragments) are
    # extra too.
    for row in replicons:
        if row['status'] == 'complete':
            extra += [c for c in row['contigs'] if c != row['contig']]
    extra_bases = sum(len(asm[c]) for c in extra)
    lengths = [len(s) for s in asm.values()]
    wrong = total_errors + total_missing
    return {
        'assembly': str(assembly), 'reference': str(reference),
        'contigs': len(asm), 'total_bases': sum(lengths), 'n50': n50(lengths),
        'replicons': len(ref),
        'complete': sum(r['status'] == 'complete' for r in replicons),
        'chromosome_complete': replicons[0]['status'] == 'complete' if replicons else False,
        'errors': total_errors, 'masked': total_masked,
        'missing_bases': total_missing, 'extra_contigs': len(extra), 'extra_bases': extra_bases,
        'misassemblies': len(breaks), 'breaks': breaks,
        'qscore': round(-10 * math.log10(max(wrong, 1) / genome), 2) if genome else 0.0,
        'perfect': wrong == 0 and not extra and not breaks and all(
            r['status'] == 'complete' for r in replicons),
        'replicon_rows': replicons,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('assembly')
    p.add_argument('reference')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--json')
    p.add_argument('--work', help='directory for alignments [default: a temporary one]')
    p.add_argument('--keep', action='store_true', help='keep the work directory')
    args = p.parse_args()
    work = args.work or tempfile.mkdtemp(prefix='bact_score.')
    try:
        m = score(args.assembly, args.reference, args.threads, work)
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
    print('replicon\tlength\tstatus\terrors\tmissing\tcontigs')
    for r in m['replicon_rows']:
        contigs = ','.join(r['contigs'][:3]) + (f",+{len(r['contigs']) - 3}"
                                                if len(r['contigs']) > 3 else '')
        print(f"{r['replicon']}\t{r['length']}\t{r['status']}\t{r['errors']}\t"
              f"{r['missing_bases']}\t{contigs}")
    for contig, kind in m['breaks']:
        print(f'misassembly\t{contig}\t{kind}')
    print(f"summary\tcomplete {m['complete']}/{m['replicons']}\terrors {m['errors']}\t"
          + (f"masked {m['masked']}\t" if m['masked'] else '') +
          f"missing {m['missing_bases']}\t"
          f"extra contigs {m['extra_contigs']} ({m['extra_bases']} bp)\t"
          f"misassemblies {m['misassemblies']}\tQ {m['qscore']}\t"
          f"{'perfect' if m['perfect'] else ''}")
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(m, f, indent=1)


if __name__ == '__main__':
    main()
