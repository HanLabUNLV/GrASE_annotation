"""Remove read-through transcripts flagged by readthrough_by_condition.py.

usage: python readthrough_filter.py <annotation.gtf> <readthrough_by_condition.tsv> <out.gtf>

Drops every transcript that uses a junction with decision "drop", unless it also uses a
junction with decision "keep" (stress-elevated read-through). Edit the decision column to
override individual junctions.
"""
import sys, re

gtf, tab, out = sys.argv[1:4]
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]
drop, keep = set(), set()
for l in open(tab):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'chrom': hdr = f; continue
    tx = set(filter(None, f[hdr.index('transcripts')].split(',')))
    (keep if f[hdr.index('decision')] == 'keep' else drop).update(tx)
drop -= keep
n = 0
with open(out, 'w') as o:
    for l in open(gtf):
        if l[0] != '#' and tid(l.split('\t')[8]) in drop: continue
        o.write(l)
print(f"dropped {len(drop)} transcripts ({len(keep)} kept for stress-elevated read-through)")
