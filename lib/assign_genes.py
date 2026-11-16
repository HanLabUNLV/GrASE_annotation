"""Assign gene_ids by exon overlap and attach reference gene ids (step 03).

usage: python assign_genes.py <in.gtf> <out.gtf> [ref.gtf] [prefix=SG]
Genes: transcripts on the same chrom and strand whose EXONS overlap (>= 1 bp) are joined,
transitively, except that with ref.gtf two groups touching different reference genes are not
joined by exon overlap alone (only a transcript touching both, e.g. a bridge kept by step 03 D2,
joins them). Spans that merely overlap (a gene inside another's intron, or
neighbours) stay separate. Replaces gffcompare's XLOC loci, which group by span overlap.
gene_id "<prefix><n>", numbered by position. Any existing gene_id/ref_gene_id is replaced.
With ref.gtf (e.g. MAKER models): ref_gene_id "<id>[,<id>...]" on each transcript (and on
its exons) = reference genes whose exons overlap the transcript's exons on the same strand,
or "novel" if none; printed: genes whose transcripts touch >= 2 reference genes.
"""
import sys, re, bisect, collections

src, out = sys.argv[1:3]
ref = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] else None
prefix = sys.argv[4] if len(sys.argv) > 4 else 'SG'
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]

ex, info = collections.defaultdict(list), {}
for l in open(src):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'exon':
        t = tid(f[8]); ex[t].append((int(f[3]), int(f[4]))); info[t] = (f[0], f[6])

refhit, rset = {}, {t: set() for t in ex}
if ref:
    rx = collections.defaultdict(list)
    for l in open(ref):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            rx[(f[0], f[6])].append((int(f[3]), int(f[4]), re.search(r'gene_id "([^"]+)"', f[8])[1]))
    idx = {}
    for k, v in rx.items():
        v.sort(); idx[k] = ([a for a, _, _ in v], v, max(b - a for a, b, _ in v))
    for t, e in ex.items():
        h = set()
        if info[t] in idx:
            st, v, mx = idx[info[t]]
            for a, b in e:
                i = bisect.bisect_left(st, a - mx)
                while i < len(v) and v[i][0] <= b:
                    if v[i][1] >= a: h.add(v[i][2])
                    i += 1
        rset[t] = h; refhit[t] = ','.join(sorted(h)) or 'novel'

# Group by exon overlap. With a reference, two groups touching DIFFERENT reference genes are not
# joined by exon overlap alone (e.g. neighbours whose UTR extensions overlap); a transcript that
# itself touches both genes (a bridge kept by the read-through filter) shares a gene with each
# side and so still joins them. Order: bridges first, then by position.
parent = {t: t for t in ex}; gref = {t: set(rset[t]) for t in ex}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
blocked = 0
def join(x, y):
    global blocked
    rx_, ry = find(x), find(y)
    if rx_ == ry: return
    if gref[rx_] and gref[ry] and not (gref[rx_] & gref[ry]): blocked += 1; return
    parent[ry] = rx_; gref[rx_] |= gref[ry]
by = collections.defaultdict(list)          # (chrom, strand) -> exons (start, end, t)
for t, e in ex.items():
    for a, b in e: by[info[t]].append((a, b, t))
pairs = []
for k, v in by.items():                     # sweep: all pairs of transcripts with overlapping exons
    v.sort(); active = []
    for a, b, t in v:
        active = [x for x in active if x[0] >= a]
        for _, u in active:
            if u != t: pairs.append((t, u))
        active.append((b, t))
pairs.sort(key=lambda p: -max(len(rset[p[0]]), len(rset[p[1]])))   # bridges (>= 2 genes) first
for x, y in pairs: join(x, y)
groups = collections.defaultdict(list)
for t in ex: groups[find(t)].append(t)
order = sorted(groups, key=lambda r: (info[r][0], min(a for t in groups[r] for a, _ in ex[t])))
gid = {}
for i, r in enumerate(order, 1):
    for t in groups[r]: gid[t] = f'{prefix}{i:06d}'
if ref:
    print(f'exon-overlap joins blocked because the two sides touch different reference genes: {blocked}')
    g2r = collections.defaultdict(set)
    for t, r in refhit.items():
        if r != 'novel': g2r[gid[t]] |= set(r.split(','))
    print(f'genes touching >= 2 reference genes: {sum(1 for v in g2r.values() if len(v) >= 2)} of {len(groups)}; '
          f'reference genes hit: {len(set().union(*g2r.values())) if g2r else 0}')

strip = re.compile(r'\s*(gene_id|ref_gene_id) "[^"]*";')
with open(out, 'w') as o:
    for l in open(src):
        if l[0] == '#': o.write(l); continue
        f = l.rstrip('\n').split('\t'); t = tid(f[8])
        attrs = strip.sub('', f[8]).strip()
        extra = f' ref_gene_id "{refhit[t]}";' if ref else ''
        f[8] = f'gene_id "{gid[t]}"; {attrs}{extra}'
        o.write('\t'.join(f) + '\n')
print(f'{len(ex)} transcripts in {len(groups)} genes ({len(ex) / len(groups):.2f} per gene)')
