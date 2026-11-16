"""Combined annotation = StringTie final + Aletsch-final transcripts that add a junction.

usage: python build_combined_gtf.py <stringtie_final.gtf> <aletsch_final.gtf> <AL_vs_ST.tmap> <out.gtf>

Added: every Aletsch transcript with >= 1 intron (chrom, strand, start, end) absent from the
StringTie annotation, except intron-retention variants of StringTie transcripts (gffcompare
class "m"). gene_id: the overlapping StringTie gene (tmap ref_gene_id) for same-strand
overlapping classes; otherwise a new id (ALG.<n>) per cluster of added transcripts with
same-strand overlapping spans. Aletsch's own gene_ids are graph-instance ids, not genes.
"""
import sys, re
from collections import defaultdict

st_gtf, al_gtf, tmap, out = sys.argv[1:5]
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]
SAME_GENE = {'=', 'c', 'j', 'k', 'm', 'n', 'o'}   # same-strand overlap with a StringTie transcript

def read(path):
    tx, ex = {}, defaultdict(list)
    for l in open(path):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        t = tid(f[8])
        if f[2] == 'transcript': tx[t] = f
        elif f[2] == 'exon': ex[t].append(f)
    for v in ex.values(): v.sort(key=lambda f: int(f[3]))
    return tx, ex

def introns(exs):
    return {(e[0], e[6], int(a[4]) + 1, int(e[3]) - 1) for a, e in zip(exs, exs[1:])}

st_tx, st_ex = read(st_gtf)
al_tx, al_ex = read(al_gtf)
st_introns = set().union(*(introns(v) for v in st_ex.values()))

cls, ref_gene = {}, {}
for l in open(tmap):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'ref_gene_id': continue
    cls[f[4]], ref_gene[f[4]] = f[2], f[0]

add = [t for t in al_ex if (introns(al_ex[t]) - st_introns) and cls.get(t) != 'm']
n_m = sum(1 for t in al_ex if (introns(al_ex[t]) - st_introns) and cls.get(t) == 'm')

# gene ids for added transcripts
gene, new = {}, []
for t in add:
    if cls.get(t) in SAME_GENE and ref_gene.get(t, '-') != '-': gene[t] = ref_gene[t]
    else: new.append(t)
spans = sorted((al_tx[t][0], al_tx[t][6], int(al_tx[t][3]), int(al_tx[t][4]), t) for t in new)
k, cur, end = 0, None, -1
for c, s, a, b, t in spans:
    if (c, s) != cur or a > end: k += 1; cur, end = (c, s), b
    end = max(end, b); gene[t] = f'ALG.{k}'

def rewrite(f, g, t, extra=''):
    f = f[:]; f[8] = f'gene_id "{g}"; transcript_id "{t}";{extra}'
    return '\t'.join(f)

records = []
for t, f in st_tx.items():
    records.append((f[0], int(f[3]), t, ['\t'.join(f)] + ['\t'.join(e) for e in st_ex[t]]))
for t in add:
    g, f = gene[t], al_tx[t]
    cov = re.search(r'cov "[^"]+"', f[8])
    lines = [rewrite(f, g, t, f' source "aletsch"; aletsch_class "{cls.get(t, "?")}";' + (f' {cov[0]};' if cov else ''))]
    for i, e in enumerate(al_ex[t], 1): lines.append(rewrite(e, g, t, f' exon_number "{i}";'))
    records.append((f[0], int(f[3]), t, lines))
records.sort()
with open(out, 'w') as o:
    o.write(f'# combined: {st_gtf} ({len(st_tx)} tx) + {len(add)} Aletsch tx with novel introns from {al_gtf}\n')
    for *_, lines in records: o.write('\n'.join(lines) + '\n')

print(f"StringTie transcripts: {len(st_tx)}")
print(f"Aletsch transcripts added: {len(add)} (skipped {n_m} class-m intron-retention variants)")
print(f"  assigned to StringTie genes: {len(add) - len(new)}; in new loci (ALG.*): {len(new)} transcripts, {k} loci")
print(f"total: {len(st_tx) + len(add)} transcripts")
