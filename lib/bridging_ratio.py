"""Classify gene-bridging junctions by read support relative to the bridged genes.

usage: python bridging_ratio.py <annotation.gtf> <maker_ref.gtf> <junction_counts_dir> <out.tsv>

A bridging intron has its donor in one MAKER gene and its acceptor in another (same strand).
ratio = reads(bridging intron) / min(median reads of the annotation's introns inside gene A,
inside gene B), reads summed over all samples (unique reads, from the BAMs).
Low ratio -> minor read-through product; ratio near 1 -> the two MAKER models are probably
fragments of one gene.
"""
import sys, re, os, bisect
from collections import defaultdict
import numpy as np

gtf, ref, jdir, out = sys.argv[1:5]
attr = lambda s, k: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]

reads = defaultdict(int)
for fn in os.listdir(jdir):
    if fn.endswith('.junc.tsv'):
        for l in open(os.path.join(jdir, fn)):
            c, s, e, n = l.split('\t'); reads[(c, int(s), int(e))] += int(n)

rg = defaultdict(list)
for l in open(ref):
    f = l.split('\t')
    if len(f) > 8 and f[1] == 'maker' and f[2] == 'gene':
        rg[(f[0], f[6])].append((int(f[3]), int(f[4]), attr(f[8], 'gene_id')))
for k in rg: rg[k].sort()
starts = {k: [g[0] for g in v] for k, v in rg.items()}
maxlen = max(e - s for v in rg.values() for s, e, _ in v)
def genes_at(key, p):
    if key not in rg: return set()
    lst, ss = rg[key], starts[key]; i = bisect.bisect_left(ss, p - maxlen); h = set()
    while i < len(lst) and lst[i][0] <= p:
        if lst[i][1] >= p: h.add(lst[i][2])
        i += 1
    return h

ex, info = defaultdict(list), {}
for l in open(gtf):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[2] == 'exon':
        t = attr(f[8], 'transcript_id'); ex[t].append((int(f[3]), int(f[4]))); info[t] = (f[0], f[6], attr(f[8], 'gene_id'))
introns = defaultdict(set)   # intron -> transcripts
for t, v in ex.items():
    v.sort(); c, st, g = info[t]
    for (a, b), (x, y) in zip(v, v[1:]): introns[(c, st, b + 1, x - 1)].add(t)

# introns inside each single MAKER gene
inside = defaultdict(list)
bridging = []
for (c, st, s, e), ts in introns.items():
    d, a = genes_at((c, st), s - 1), genes_at((c, st), e + 1)
    if d and a and not (d & a): bridging.append(((c, st, s, e), ts, d, a))
    elif d & a:
        for g in d & a: inside[g].append(reads[(c, s, e)])
med = {g: float(np.median(v)) for g, v in inside.items()}

rows = []
for (c, st, s, e), ts, d, a in bridging:
    md = max((med.get(g, 0) for g in d), default=0); ma = max((med.get(g, 0) for g in a), default=0)
    base = min(md, ma)
    r = reads[(c, s, e)]
    ratio = r / base if base > 0 else float('inf')
    rows.append((c, st, s, e, r, md, ma, ratio, ','.join(sorted(d)), ','.join(sorted(a)), len(ts), ','.join(sorted(ts))[:200]))

with open(out, 'w') as o:
    o.write('chrom\tstrand\tstart\tend\tbridge_reads\tmedian_reads_in_donor_gene\tmedian_reads_in_acceptor_gene\tratio\tdonor_maker_gene\tacceptor_maker_gene\tn_tx\ttranscripts\n')
    for r in sorted(rows): o.write('\t'.join(f'{x:.3f}' if isinstance(x, float) else str(x) for x in r) + '\n')

R = np.array([r[7] for r in rows])
print(f"bridging introns: {len(rows)}")
for lo, hi, lab in [(0, 0.05, '< 0.05 (minor read-through)'), (0.05, 0.2, '0.05-0.2'), (0.2, 0.5, '0.2-0.5'), (0.5, np.inf, '>= 0.5 (comparable: likely split MAKER model)')]:
    m = (R >= lo) & (R < hi); print(f"  ratio {lab:<48} {int(m.sum()):>6}")
print(f"  no internal introns to compare (inf)               {int(np.isinf(R).sum()):>6}")
