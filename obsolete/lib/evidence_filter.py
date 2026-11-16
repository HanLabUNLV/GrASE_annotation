"""Evidence rule: drop transcripts never assembled as such, unless needed for a junction.

usage: python evidence_filter.py <in.gtf> <out.gtf> <support_table.tsv>...

A transcript (any source) is a candidate when no sample's own assembly contains its exact
intron chain (gffcompare "=", column `total` = 0 in the step 01/03 support tables).
Candidates are removed lowest coverage first ("cov" attribute; 0 if absent), each only if
every one of its junctions is still present in another remaining transcript of the same
gene, so no junction is lost. Single-exon candidates have no junctions and are removed.
Rationale: short reads support junctions, rarely whole intron chains; Aletsch reports many
combinations of the same junctions from its pooled graph (42% of the redundant ones were
never assembled by any sample), which multiply GrASE's paths without adding junctions.
"""
import sys, re
from collections import defaultdict

src, out, *tables = sys.argv[1:]
total = {}
for p in tables:
    for l in open(p):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'transcript_id': total[f[0]] = int(f[2])

tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]
ex, gene, cov = defaultdict(list), {}, {}
for l in open(src):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t'); t = tid(f[8])
    if f[2] == 'transcript':
        gene[t] = re.search(r'gene_id "([^"]+)"', f[8])[1]
        m = re.search(r'cov "([^"]+)"', f[8]); cov[t] = float(m[1]) if m else 0.0
    elif f[2] == 'exon':
        ex[t].append((int(f[3]), int(f[4])))
missing = [t for t in gene if t not in total]
if missing: sys.exit(f'{len(missing)} transcripts have no support record, e.g. {missing[:3]}')
J = {}
for t, v in ex.items():
    v.sort(); J[t] = [(b + 1, x - 1) for (a, b), (x, y) in zip(v, v[1:])]
bygene = defaultdict(list)
for t, g in gene.items(): bygene[g].append(t)

drop = set()
for g, ts in bygene.items():
    n = defaultdict(int)
    for t in ts:
        for j in J[t]: n[j] += 1
    for t in sorted([t for t in ts if total[t] == 0], key=lambda t: cov[t]):
        if all(n[j] > 1 for j in J[t]):
            drop.add(t)
            for j in J[t]: n[j] -= 1

with open(out, 'w') as o:
    for l in open(src):
        if l[0] != '#' and tid(l.split('\t')[8]) in drop: continue
        o.write(l)
cand = sum(1 for t in gene if total[t] == 0)
print(f"evidence rule: {len(gene)} transcripts, {cand} never assembled exactly in any sample, "
      f"{len(drop)} removed ({sum(len(J[t]) == 0 for t in drop)} single-exon), {len(gene) - len(drop)} kept; no junction lost")
