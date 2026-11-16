"""Give each disjoint piece of a gene its own gene_id.

usage: python split_disjoint_genes.py <in.gtf> <out.gtf>

Transcripts sharing a gene_id but not connected by overlapping spans (same chrom/strand)
are split into pieces ordered by position; the first keeps the gene_id, later pieces get
"<gene_id>.p2", ".p3", ... (left behind when filters removed the transcripts linking them).
"""
import sys, re
from collections import defaultdict

src, out = sys.argv[1:3]
gid = lambda s: re.search(r'gene_id "([^"]+)"', s)[1]
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]
spans = defaultdict(list)
for l in open(src):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'transcript': spans[gid(f[8])].append((int(f[3]), int(f[4]), tid(f[8])))
new, nsplit = {}, 0
for g, v in spans.items():
    v.sort(); piece, end = 1, -1
    for a, b, t in v:
        if end >= 0 and a > end: piece += 1
        end = max(end, b); new[t] = g if piece == 1 else f'{g}.p{piece}'
    nsplit += piece > 1
with open(out, 'w') as o:
    for l in open(src):
        if l[0] != '#':
            f = l.split('\t'); t = tid(f[8]); g = gid(f[8])
            if new[t] != g: f[8] = f[8].replace(f'gene_id "{g}"', f'gene_id "{new[t]}"', 1); l = '\t'.join(f)
        o.write(l)
print(f"genes split: {nsplit}; gene_ids {len(spans)} -> {len(set(new.values()))}")
