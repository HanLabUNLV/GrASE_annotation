"""Why are original significant SE events not tested with a new annotation?

usage: python diagnose_missing_se.py <combined.gtf> <junction_counts_dir> <orig_filtered_SE.csv>... -- <new_SE.MATS.JC.txt>...
For each original significant SE event not tested in the new run (pooled over files), check its 3
junctions (inclusion up/down, skipping): present in the new annotation? reads (unique, all samples)?
"""
import sys, os, re
from collections import defaultdict
import pandas as pd
gtf, jdir, *rest = sys.argv[1:]
i = rest.index('--'); orig, new = rest[:i], rest[i+1:]
ann = set()
ex = defaultdict(list)
for l in open(gtf):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'exon': ex[re.search(r'transcript_id "([^"]+)"', f[8])[1]].append((f[0], int(f[3]), int(f[4])))
for v in ex.values():
    v.sort(key=lambda x: x[1])
    for (c, a, b), (_, x, y) in zip(v, v[1:]): ann.add((c, b + 1, x - 1))
reads = defaultdict(int); nsamp = defaultdict(int)
for fn in os.listdir(jdir):
    if fn.endswith('.junc.tsv'):
        for l in open(os.path.join(jdir, fn)):
            c, s, e, n = l.split('\t'); k = (c, int(s), int(e)); reads[k] += int(n); nsamp[k] += int(n) >= 3
K = ['chr', 'strand', 'exonStart_0base', 'exonEnd', 'upstreamEE', 'downstreamES']
tested = set()
for p in new: tested |= set(map(tuple, pd.read_csv(p, sep='\t')[K].values))
o = pd.concat([pd.read_csv(p) for p in orig]).drop_duplicates(K)
miss = o[[tuple(r) not in tested for r in o[K].values]]
cat = defaultdict(int)
for r in miss.itertuples():
    js = [(r.chr, r.upstreamEE + 1, r.exonStart_0base), (r.chr, r.exonEnd + 1, r.downstreamES - 0), (r.chr, r.upstreamEE + 1, r.downstreamES)]
    js[1] = (r.chr, r.exonEnd + 1, r.downstreamES)   # rMATS: downstreamES is 0-based start -> intron end = downstreamES
    inann = [j in ann for j in js]
    weak = [nsamp[j] < 3 for j in js]
    if all(inann): cat['all 3 junctions annotated (event structure differs)'] += 1
    elif any((not a) and w for a, w in zip(inann, weak)): cat['missing junction is weakly supported (<3 samples >=3 reads)'] += 1
    else: cat['missing junction is well supported (>=3 samples >=3 reads)'] += 1
# where was each well-supported missing junction lost?
def introns_of(path):
    ex = defaultdict(list); out = set()
    for l in open(path):
        if l[0] == '#': continue
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[re.search(r'transcript_id "([^"]+)"', f[8])[1]].append((f[0], int(f[3]), int(f[4])))
    for v in ex.values():
        v.sort(key=lambda x: x[1])
        for (c, a, b), (_, x, y) in zip(v, v[1:]): out.add((c, b + 1, x - 1))
    return out
stages = [(n, introns_of(p)) for n, p in [('f0.01 merge', 'merged.f0.01.gtf'), ('support3pertaxon', 'merged.support3pertaxon.gtf'),
          ('StringTie final', 'merged.support3pertaxon.canonical.sj3x3.gtf'), ('Aletsch raw', 'assembly_scallop2_aletsch/aletsch/aletsch.gtf'),
          ('Aletsch support3', 'assembly_scallop2_aletsch/aletsch.support3pertaxon.gtf'), ('combined', 'merged.combined.gtf')]]
lost = defaultdict(int); seen = set()
for r in miss.itertuples():
    for j in [(r.chr, r.upstreamEE + 1, r.exonStart_0base), (r.chr, r.exonEnd + 1, r.downstreamES), (r.chr, r.upstreamEE + 1, r.downstreamES)]:
        if j in ann or nsamp[j] < 3 or j in seen: continue
        seen.add(j)
        lost[' / '.join(n for n, st in stages if j in st) or 'in none of them'] += 1
print("well-supported missing junctions, by which annotations contain them:")
for k, v in sorted(lost.items(), key=lambda x: -x[1]): print(f"  {v:5d}  {k}")
print(f"original significant SE (unique across contrasts): {len(o)}; not tested in new run: {len(miss)}")
for k, v in sorted(cat.items(), key=lambda x: -x[1]): print(f"  {v:5d}  {k}")
