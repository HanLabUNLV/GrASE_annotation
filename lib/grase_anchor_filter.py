"""Overhang (anchor) filter for GrASE bipartition events (step 11 D; post-hoc filter for
annotations built before the step 01 overhang rule, 2026-11-12).

A junction passes when >= min_samples samples have >= min_reads unique reads AND >= 1 read
with >= anchor_min aligned bases on both sides (per-sample counts from
lib/junction_anchor_counts.sh). Introns next to an internal exon shorter than anchor_min in
some transcript of the GTF (micro-exon) are exempt. An event (row of a bipartition file) is
dropped when any intron edge on any path of either side (path1, path2) is a failing junction.
Annotations built with the step 01 overhang rule contain no failing junctions, so nothing is
dropped for them.

usage: python grase_anchor_filter.py <grase_workdir> <annotation.gtf> <anchor_dir>
                                     [min_samples=3] [min_reads=3] [anchor_min=12]
reads    <wd>/bipartition.filtered/*.txt, <wd>/graphml/<gene>.graphml
writes   <wd>/bipartition.filtered.anc<anchor_min>/*.txt  (same files, failing events removed)
         <wd>/bipartition.filtered.anc<anchor_min>.removed.tsv  (removed events + failing junctions)
Graph intron edges: junction = min(from_pos, to_pos) .. max(from_pos, to_pos) - 1 (both strands).
"""
import sys, os, glob, re
import xml.etree.ElementTree as ET
from collections import defaultdict

wd, gtf, adir = sys.argv[1:4]
ms = int(sys.argv[4]) if len(sys.argv) > 4 else 3
mr = int(sys.argv[5]) if len(sys.argv) > 5 else 3
amin = int(sys.argv[6]) if len(sys.argv) > 6 else 12
NS = '{http://graphml.graphdrawing.org/xmlns}'

# samples with >= mr reads and >= 1 anchored read, per junction
nanch = defaultdict(int)
fns = glob.glob(os.path.join(adir, '*.junc_anchor.tsv'))
if not fns: sys.exit(f'no *.junc_anchor.tsv in {adir}')
for fn in fns:
    for l in open(fn):
        c, s, e, n, a = l.split('\t')
        if int(n) >= mr and int(a) >= 1: nanch[(c, int(s), int(e))] += 1

# introns next to an internal micro-exon (any transcript): exempt
tx = defaultdict(list)
for l in open(gtf):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] != 'exon': continue
    tx[re.search(r'transcript_id "([^"]+)"', f[8])[1]].append((f[0], int(f[3]), int(f[4])))
exempt = set()
for ex in tx.values():
    ex.sort(key=lambda x: x[1])
    for i in range(1, len(ex) - 1):
        c, s, e = ex[i]
        if e - s + 1 < amin:
            exempt.add((c, ex[i-1][2] + 1, s - 1)); exempt.add((c, e + 1, ex[i+1][1] - 1))
del tx

def fails(j): return j not in exempt and nanch[j] < ms

graphs = {}
def intron_edges(gene):
    """(source name, target name) -> junction for every intron edge of the gene's graph."""
    if gene not in graphs:
        g = ET.parse(os.path.join(wd, 'graphml', f'{gene}.graphml')).getroot().find(NS + 'graph')
        chrom = {d.get('key'): d.text for d in g.findall(NS + 'data')}['g_chrom']
        name = {n.get('id'): {d.get('key'): d.text for d in n.findall(NS + 'data')}['v_name'] for n in g.findall(NS + 'node')}
        ed = {}
        for e in g.findall(NS + 'edge'):
            d = {x.get('key'): x.text for x in e.findall(NS + 'data')}
            if d.get('e_ex_or_in') != 'in': continue
            a, b = int(d['e_from_pos']), int(d['e_to_pos'])
            ed[(name[e.get('source')], name[e.get('target')])] = (chrom, min(a, b), max(a, b) - 1)
        graphs.clear(); graphs[gene] = ed   # files are per gene: keep one graph in memory
    return graphs[gene]

src = os.path.join(wd, 'bipartition.filtered')
dst = os.path.join(wd, f'bipartition.filtered.anc{amin}'); os.makedirs(dst, exist_ok=True)
nrow = nrm = nfile = 0
with open(f'{dst}.removed.tsv', 'w') as rm:
    rm.write('file\tgene\tsource\tsink\tfailing_junctions\n')
    for fn in sorted(glob.glob(os.path.join(src, '*.txt'))):
        lines = open(fn).read().splitlines(keepends=True)
        if not lines: continue
        h = lines[0].rstrip('\n').split('\t'); ig, ip1, ip2 = h.index('gene'), h.index('path1'), h.index('path2')
        out = [lines[0]]
        for line in lines[1:]:
            f = line.rstrip('\n').split('\t'); nrow += 1
            ed = intron_edges(f[ig]); bad = set()
            for p in f[ip1].split(',') + f[ip2].split(','):
                nodes = p.strip().split('-')
                for u, v in zip(nodes, nodes[1:]):
                    j = ed.get((u, v))
                    if j and fails(j): bad.add(j)
            if bad:
                nrm += 1
                rm.write(f"{os.path.basename(fn)}\t{f[ig]}\t{f[h.index('source')]}\t{f[h.index('sink')]}\t"
                         + ','.join(f'{c}:{s}-{e}' for c, s, e in sorted(bad)) + '\n')
            else:
                out.append(line)
        open(os.path.join(dst, os.path.basename(fn)), 'w').writelines(out); nfile += 1
print(f'{nfile} files, {nrow} events, {nrm} removed (junction with < {ms} samples having >= {mr} reads and '
      f'>= 1 read anchored >= {amin} bp each side; {len(exempt)} micro-exon introns exempt) -> {dst}')
