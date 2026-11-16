"""Junction QC: intron length distribution and junctions/transcripts bridging genes.
usage: python junction_qc.py merged.gtf reference.gtf outprefix
"""
import sys, re, bisect
from collections import defaultdict
import numpy as np

gtf, ref, out = sys.argv[1:4]
attr = lambda s, k: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]

# --- merged transcripts
ex = defaultdict(list); info = {}
for l in open(gtf):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[2] != 'exon': continue
    t = attr(f[8], 'transcript_id')
    ex[t].append((int(f[3]), int(f[4])))
    info[t] = (f[0], f[6], attr(f[8], 'gene_id'))

# --- reference (MAKER) genes and their exons, per chrom/strand
rg = defaultdict(list); rex = defaultdict(list)
for l in open(ref):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[1] != 'maker': continue
    if f[2] == 'gene':
        rg[(f[0], f[6])].append((int(f[3]), int(f[4]), attr(f[8], 'gene_id')))
    elif f[2] == 'exon':
        rex[(f[0], f[6])].append((int(f[3]), int(f[4]), attr(f[8], 'gene_id')))
for d in (rg, rex):
    for k in d: d[k].sort()
starts = {k: [g[0] for g in v] for k, v in rg.items()}
exstarts = {k: [g[0] for g in v] for k, v in rex.items()}
maxlen = max(e - s for v in rg.values() for s, e, _ in v)
maxexlen = max(e - s for v in rex.values() for s, e, _ in v)

def genes_overlapping(key, s, e, use_exons=False):
    arr, st, ml = (rex, exstarts, maxexlen) if use_exons else (rg, starts, maxlen)
    if key not in arr: return set()
    lst, ss = arr[key], st[key]
    i = bisect.bisect_left(ss, s - ml)
    hits = set()
    while i < len(lst) and lst[i][0] <= e:
        if lst[i][1] >= s: hits.add(lst[i][2])
        i += 1
    return hits

def genes_contained(key, s, e):
    if key not in rg: return set()
    lst, ss = rg[key], starts[key]
    i = bisect.bisect_left(ss, s)
    hits = set()
    while i < len(lst) and lst[i][0] <= e:
        if lst[i][1] <= e: hits.add(lst[i][2])
        i += 1
    return hits

introns = {}  # (chr,strand,s,e) -> list of tx
tx_multi = []
for t, exons in ex.items():
    exons.sort(); c, st, g = info[t]; key = (c, st)
    for (a, b), (x, y) in zip(exons, exons[1:]):
        introns.setdefault((c, st, b + 1, x - 1), []).append(t)
    # reference genes whose EXONS this transcript's exons overlap
    hit = set()
    for a, b in exons: hit |= genes_overlapping(key, a, b, use_exons=True)
    if len(hit) >= 2:
        tx_multi.append((t, g, c, st, exons[0][0], exons[-1][1], len(exons), ','.join(sorted(hit))))

L = np.array([k[3] - k[2] + 1 for k in introns])
print(f"merged transcripts: {len(ex)}   unique introns: {len(L)}")
print("intron length quantiles (bp): " + "  ".join(f"q{q}={int(np.percentile(L, q))}" for q in (50, 90, 99, 99.9)) + f"  max={L.max()}")
for thr in (5000, 10000, 20000, 50000, 100000):
    print(f"  introns > {thr:>6} bp: {int((L > thr).sum())}")

# reference intron lengths for comparison
rL = []
mrna_ex = defaultdict(list)
for l in open(ref):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[1] == 'maker' and f[2] == 'exon':
        mrna_ex[attr(f[8], 'transcript_id')].append((int(f[3]), int(f[4])))
for v in mrna_ex.values():
    v.sort(); rL += [x - b - 1 for (a, b), (x, y) in zip(v, v[1:])]
rL = np.array(rL)
print("MAKER intron length quantiles: " + "  ".join(f"q{q}={int(np.percentile(rL, q))}" for q in (50, 90, 99, 99.9)) + f"  max={rL.max()}")

# junctions bridging genes: donor & acceptor in different ref genes, or intron containing whole ref gene(s)
bridge = []
for (c, st, s, e), ts in introns.items():
    key = (c, st)
    don = genes_overlapping(key, s - 1, s - 1); acc = genes_overlapping(key, e + 1, e + 1)
    inside = genes_contained(key, s, e)
    cross = bool(don and acc and not (don & acc))
    if cross or inside:
        bridge.append((c, st, s, e, e - s + 1, len(ts), ','.join(sorted(don)) or '-', ','.join(sorted(acc)) or '-',
                       ','.join(sorted(inside)) or '-', ','.join(sorted({info[t][2] for t in ts})), ','.join(ts[:5])))
ncross = sum(1 for b in bridge if b[6] != '-' and b[7] != '-' and not set(b[6].split(',')) & set(b[7].split(',')))
ncont = sum(1 for b in bridge if b[8] != '-')
print(f"\nintrons whose donor and acceptor lie in different MAKER genes: {ncross}")
print(f"introns that contain an entire MAKER gene (same strand): {ncont}")
print(f"union: {len(bridge)} introns in {len({g for b in bridge for g in b[9].split(',')})} merged genes")
print(f"\ntranscripts whose exons overlap exons of >=2 MAKER genes (read-through/fusion candidates): {len(tx_multi)} "
      f"in {len({x[1] for x in tx_multi})} merged genes")

with open(out + '.bridging_introns.tsv', 'w') as o:
    o.write('chrom\tstrand\tintron_start\tintron_end\tlength\tn_tx\tdonor_gene\tacceptor_gene\tcontained_genes\tmerged_genes\texample_tx\n')
    for b in sorted(bridge): o.write('\t'.join(map(str, b)) + '\n')
with open(out + '.multigene_transcripts.tsv', 'w') as o:
    o.write('transcript_id\tgene_id\tchrom\tstrand\tstart\tend\tn_exons\tmaker_genes\n')
    for x in sorted(tx_multi, key=lambda x: (x[2], x[4])): o.write('\t'.join(map(str, x)) + '\n')
with open(out + '.long_introns.tsv', 'w') as o:
    o.write('chrom\tstrand\tintron_start\tintron_end\tlength\tn_tx\tmerged_genes\n')
    for (c, st, s, e), ts in sorted(introns.items()):
        if e - s + 1 > 20000:
            o.write(f"{c}\t{st}\t{s}\t{e}\t{e-s+1}\t{len(ts)}\t{','.join(sorted({info[t][2] for t in ts}))}\n")
