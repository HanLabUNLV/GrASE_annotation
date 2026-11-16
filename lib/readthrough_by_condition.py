"""Per-condition read-through ratios for minor gene-bridging junctions.

usage: python readthrough_by_condition.py <annotation.gtf> <maker_ref.gtf> <bridging_ratio.tsv>
                                          <sample_conditions.tsv> <junction_counts_dir> <out.tsv>
                                          [max_ratio=0.2] [fold=3] [min_reads=10]

Candidates: bridging introns (donor and acceptor in different MAKER genes) with pooled
ratio < max_ratio (see bridging_ratio.py). For each candidate and each taxon x condition
group, ratio = bridge reads / min(median reads of the annotation's introns inside the
donor gene, inside the acceptor gene), reads summed over the group's samples.
A candidate is kept as stress-elevated if, in some taxon, a stress group (drought or
waterlogging) has ratio >= max_ratio AND >= fold x the same taxon's control ratio
(control ratio floored at 0.01) with >= min_reads bridging reads summed over that group. Output: one row per candidate junction with the ratio
per taxon:condition group, the stress groups where it is elevated, the transcripts using
it, and decision keep (stress-elevated) / drop. Apply with readthrough_filter.py.
"""
import sys, re, os
from collections import defaultdict
import numpy as np

gtf, ref, bridge_tsv, cond_tsv, jdir, out_tab = sys.argv[1:7]
max_ratio = float(sys.argv[7]) if len(sys.argv) > 7 else 0.2
fold = float(sys.argv[8]) if len(sys.argv) > 8 else 3.0
min_reads = int(sys.argv[9]) if len(sys.argv) > 9 else 10
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]

group = {}
for l in open(cond_tsv):
    s, t, c = l.rstrip('\n').split('\t'); group[s] = (t, c)
groups = sorted(set(group.values()))

# candidates and the genes they bridge
cand = {}
for l in open(bridge_tsv):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'chrom' or float(f[7]) >= max_ratio: continue
    cand[(f[0], f[1], int(f[2]), int(f[3]))] = (f[8].split(','), f[9].split(','))

# annotation introns -> transcripts; introns inside each MAKER gene (reuse bridging table genes)
ex = defaultdict(list); strand = {}
for l in open(gtf):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'exon': t = tid(f[8]); ex[t].append((int(f[3]), int(f[4]), f[0])); strand[t] = f[6]
intr_tx = defaultdict(set)
for t, v in ex.items():
    v.sort()
    for (a, b, c), (x, y, _) in zip(v, v[1:]): intr_tx[(c, strand[t], b + 1, x - 1)].add(t)

# MAKER gene spans, for genes named in candidates
need = {g for d, a in cand.values() for g in d + a}
span = {}
for l in open(ref):
    f = l.split('\t')
    if len(f) > 8 and f[1] == 'maker' and f[2] == 'gene':
        g = re.search(r'gene_id "([^"]+)"', f[8])[1]
        if g in need: span[g] = (f[0], f[6], int(f[3]), int(f[4]))
inside = defaultdict(list)   # gene -> its internal annotation introns
by_cs = defaultdict(list)
for k in intr_tx: by_cs[(k[0], k[1])].append(k)
for g, (gc, gs, a, b) in span.items():
    inside[g] = [k for k in by_cs[(gc, gs)] if a <= k[2] - 1 and k[3] + 1 <= b and k not in cand]

# per-group read counts for the introns we need
want = set(cand) | {k for v in inside.values() for k in v}
want_key = {(k[0], k[2], k[3]) for k in want}
reads = defaultdict(lambda: defaultdict(int))   # (c,s,e) -> group -> reads
for s, gk in group.items():
    p = os.path.join(jdir, f'{s}.junc.tsv')
    for l in open(p):
        c, a, b, n = l.split('\t'); k = (c, int(a), int(b))
        if k in want_key: reads[k][gk] += int(n)

def ratio(k, d, a, gk):
    meds = []
    for genes in (d, a):
        m = max((np.median([reads[(i[0], i[2], i[3])][gk] for i in inside[g]]) if inside[g] else 0 for g in genes), default=0)
        meds.append(m)
    base = min(meds)
    r = reads[(k[0], k[2], k[3])][gk]
    return r / base if base > 0 else (np.inf if r > 0 else np.nan)

rows, keep = [], set()
for k, (d, a) in cand.items():
    rs = {gk: ratio(k, d, a, gk) for gk in groups}
    elevated = []
    for (t, c) in groups:
        if c == 'control': continue
        r, rc = rs[(t, c)], rs.get((t, 'control'), np.nan)
        if np.isnan(r) or np.isnan(rc): continue
        if r >= max_ratio and r >= fold * max(rc, 0.01) and reads[(k[0], k[2], k[3])][(t, c)] >= min_reads:
            elevated.append(f'{t}:{c}')
    if elevated: keep.add(k)
    rows.append((k, rs, elevated))

with open(out_tab, 'w') as o:
    o.write('chrom\tstrand\tstart\tend\t' + '\t'.join(f'{t}:{c}' for t, c in groups)
            + '\tstress_elevated_in\tdecision\ttranscripts\n')
    for k, rs, el in sorted(rows):
        o.write('\t'.join(map(str, k)) + '\t' + '\t'.join(f'{rs[g]:.3f}' for g in groups) + '\t'
                + (','.join(el) or '-') + '\t' + ('keep' if k in keep else 'drop') + '\t'
                + ','.join(sorted(intr_tx.get(k, ()))) + '\n')

print(f"read-through candidate junctions (pooled ratio < {max_ratio}): {len(cand)}")
print(f"  stress-elevated (keep): {len(keep)}   drop: {len(cand) - len(keep)}")
multi = sum(1 for k, rs, el in rows if len(el) >= 2)
print(f"  elevated in >= 2 taxon:stress groups: {multi}")
