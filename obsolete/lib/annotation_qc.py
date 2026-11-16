"""Side-by-side QC of transcript annotations (GTFs).

usage: python annotation_qc.py <genome.fa> <maker_ref.gtf> <junction_counts_dir> <label=gtf> [<label=gtf> ...] > table.tsv

Per GTF: transcripts, distinct intron chains, gene_ids, loci (same-strand
span-overlap clusters), transcripts per locus, single-exon transcripts, unique introns,
intron length quantiles, introns >20/100 kb, non-canonical introns (not GT-AG/GC-AG/AT-AC
on the transcript strand), introns with weak read support (<3 samples with >=3 unique
spliced reads, counted from the BAMs; see support/junction_counts), gene-bridging introns
(donor and acceptor in different MAKER genes, or containing a whole MAKER gene), and
transcripts overlapping exons of >=2 MAKER genes.
"""
import sys, re, os, bisect
from collections import defaultdict
import numpy as np

genome, ref, jdir, *specs = sys.argv[1:]
attr = lambda s, k: (re.search(k + r' "([^"]+)"', s) or [None, None])[1]

def read_gtf(path):
    ex, info = defaultdict(list), {}
    for l in open(path):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        t = attr(f[8], 'transcript_id')
        ex[t].append((int(f[3]), int(f[4])))
        info[t] = (f[0], f[6], attr(f[8], 'gene_id'))
    for v in ex.values(): v.sort()
    return ex, info

# genome
seq, name, buf = {}, None, []
for l in open(genome):
    if l[0] == '>':
        if name: seq[name] = ''.join(buf)
        name, buf = l[1:].split()[0], []
    else: buf.append(l.strip().upper())
seq[name] = ''.join(buf)
comp = str.maketrans('ACGT', 'TGCA')
def motif(c, st, s, e):
    m = seq[c][s-1:s+1] + '-' + seq[c][e-2:e]
    return m[::-1].translate(comp) if st == '-' else m
CANON = {'GT-AG', 'GC-AG', 'AT-AC'}

# junction read support from the BAM-derived per-sample counts
nsamp = defaultdict(int)
for fn in os.listdir(jdir):
    if not fn.endswith('.junc.tsv'): continue
    for l in open(os.path.join(jdir, fn)):
        c, s, e, n = l.split('\t')
        if int(n) >= 3: nsamp[(c, int(s), int(e))] += 1

# MAKER genes and exons
rg, rex = defaultdict(list), defaultdict(list)
for l in open(ref):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[1] != 'maker': continue
    if f[2] == 'gene': rg[(f[0], f[6])].append((int(f[3]), int(f[4]), attr(f[8], 'gene_id')))
    elif f[2] == 'exon': rex[(f[0], f[6])].append((int(f[3]), int(f[4]), attr(f[8], 'gene_id')))
for d in (rg, rex):
    for k in d: d[k].sort()
ml = {id(rg): max(e - s for v in rg.values() for s, e, _ in v), id(rex): max(e - s for v in rex.values() for s, e, _ in v)}
starts = {id(d): {k: [x[0] for x in v] for k, v in d.items()} for d in (rg, rex)}
def overlapping(d, key, s, e):
    if key not in d: return set()
    lst, ss = d[key], starts[id(d)][key]
    i, hits = bisect.bisect_left(ss, s - ml[id(d)]), set()
    while i < len(lst) and lst[i][0] <= e:
        if lst[i][1] >= s: hits.add(lst[i][2])
        i += 1
    return hits
def contained(key, s, e):
    if key not in rg: return set()
    lst, ss = rg[key], starts[id(rg)][key]
    i, hits = bisect.bisect_left(ss, s), set()
    while i < len(lst) and lst[i][0] <= e:
        if lst[i][1] <= e: hits.add(lst[i][2])
        i += 1
    return hits

rows = []
for spec in specs:
    label, path = spec.split('=', 1)
    ex, info = read_gtf(path)
    genes = {v[2] for v in info.values()}
    # loci: same-strand transcripts with overlapping genomic spans (independent of gene_id)
    spans = sorted((info[t][0], info[t][1], v[0][0], v[-1][1]) for t, v in ex.items())
    loci, sizes, cur, end, n_in = 0, [], None, -1, 0
    for c, st, s0, e0 in spans:
        if (c, st) != cur or s0 > end:
            if n_in: sizes.append(n_in)
            loci += 1; cur, end, n_in = (c, st), e0, 0
        end = max(end, e0); n_in += 1
    sizes.append(n_in)
    chains = {(info[t][0], info[t][1], tuple(v)) if len(v) == 1 else (info[t][0], info[t][1], tuple((b, x) for (a, b), (x, y) in zip(v, v[1:]))) for t, v in ex.items()}
    introns, multi = set(), 0
    for t, exons in ex.items():
        c, st, g = info[t]
        for (a, b), (x, y) in zip(exons, exons[1:]): introns.add((c, st, b + 1, x - 1))
        hit = set()
        for a, b in exons: hit |= overlapping(rex, (c, st), a, b)
        multi += len(hit) >= 2
    L = np.array([e - s + 1 for c, st, s, e in introns])
    nc = weak = bridge = 0
    for c, st, s, e in introns:
        nc += motif(c, st, s, e) not in CANON
        weak += nsamp.get((c, s, e), 0) < 3
        don, acc = overlapping(rg, (c, st), s - 1, s - 1), overlapping(rg, (c, st), e + 1, e + 1)
        bridge += bool((don and acc and not (don & acc)) or contained((c, st), s, e))
    n, ni = len(ex), len(introns)
    pct = lambda x, d: f"{x} ({100 * x / d:.1f}%)" if d else "0"
    rows.append((label, {
        'transcripts': n, 'distinct intron chains': len(chains), 'gene_ids': len(genes),
        'loci (same-strand span overlap)': loci, 'transcripts/locus mean / max': f"{n / loci:.2f} / {max(sizes)}",
        'single-exon transcripts': pct(sum(len(v) == 1 for v in ex.values()), n),
        'unique introns': ni,
        'intron length median / q99 / max': f"{int(np.median(L))} / {int(np.percentile(L, 99))} / {L.max()}",
        'introns > 20 kb': pct(int((L > 20000).sum()), ni), 'introns > 100 kb': pct(int((L > 100000).sum()), ni),
        'non-canonical introns': pct(nc, ni),
        'weakly supported introns (<3 samples >=3 reads)': pct(weak, ni),
        'gene-bridging introns (MAKER)': pct(bridge, ni),
        'transcripts spanning >=2 MAKER genes': pct(multi, n)}))

keys = list(rows[0][1])
print('metric\t' + '\t'.join(r[0] for r in rows))
for k in keys: print(k + '\t' + '\t'.join(str(r[1][k]) for r in rows))
