"""Match significant rMATS events between two runs by event coordinates.

usage: python compare_rmats_events.py <filtered_dir_A> <rmats_out_B> <filtered_dir_B>
For each event type: significant in A; of those, tested in B (same coordinates); significant in B.
"""
import sys, os, pandas as pd
fa, ob, fb = sys.argv[1:4]
COORD = {'SE': ['chr','strand','exonStart_0base','exonEnd','upstreamEE','downstreamES'],
         'A5SS': ['chr','strand','longExonStart_0base','longExonEnd','shortES','shortEE','flankingES','flankingEE'],
         'A3SS': ['chr','strand','longExonStart_0base','longExonEnd','shortES','shortEE','flankingES','flankingEE'],
         'MXE': ['chr','strand','1stExonStart_0base','1stExonEnd','2ndExonStart_0base','2ndExonEnd','upstreamEE','downstreamES'],
         'RI': ['chr','strand','riExonStart_0base','riExonEnd','upstreamEE','downstreamES']}
def key(df, c, ev):
    """Splice-site key: events defined by the same junction coordinates match even if
    flanking exons' far ends (which depend on transcript ends) differ."""
    if ev in ('A5SS', 'A3SS'):   # alt exon ends + the flanking exon's splice site
        fl = df.flankingES.where(df.flankingES > df.longExonEnd, df.flankingEE)
        d = df[['chr','strand','longExonStart_0base','longExonEnd','shortES','shortEE']].copy(); d['fl'] = fl
        return set(map(tuple, d.astype(str).values))
    if ev == 'RI':               # the retained intron
        return set(map(tuple, df[['chr','strand','upstreamEE','downstreamES']].astype(str).values))
    return set(map(tuple, df[c].astype(str).values))
for ev, c in COORD.items():
    a = key(pd.read_csv(os.path.join(fa, f'{ev}.MATS.JC.filtered.csv')), c, ev)
    t = key(pd.read_csv(os.path.join(ob, f'{ev}.MATS.JC.txt'), sep='\t'), c, ev)
    b = key(pd.read_csv(os.path.join(fb, f'{ev}.MATS.JC.filtered.csv')), c, ev)
    print(f"{ev}\t{len(a)}\t{len(a & t)}\t{len(a & b)}\t{len(b)}\t{len(b - a)}")
