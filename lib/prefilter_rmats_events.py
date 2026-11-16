"""Filter an rMATS event set (fromGTF.<event>.txt) by the junction rules before testing.

usage: python prefilter_rmats_events.py <rmats_out_dir> <events_out_dir> <junc_dir> <genome.fa>
                                        [max_intron=20000] [min_samples=3] [min_reads=3] [anchor_dir]
anchor_dir (optional; "" = off): per-sample overhang counts (lib/junction_anchor_counts.sh) for the
overhang rule in rmats_junction_rules.py (>= 1 read with >= 12 bp each side in each supporting sample).
Reads <rmats_out_dir>/fromGTF.{SE,MXE,A3SS,A5SS,RI}.txt (all events detected, incl. --novelSS),
writes the passing rows, unchanged, to <events_out_dir>/fromGTF.<event>.txt for
rmats.py --fixed-event-set.
"""
import sys, os, csv
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rmats_junction_rules import JunctionRules

src, dst, jdir, genome = sys.argv[1:5]
mx = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
ms = int(sys.argv[6]) if len(sys.argv) > 6 else 3
mr = int(sys.argv[7]) if len(sys.argv) > 7 else 3
anc = sys.argv[8] if len(sys.argv) > 8 and sys.argv[8] else None
rules = JunctionRules(jdir, genome, mx, ms, mr, anc)
os.makedirs(dst, exist_ok=True)
for ev in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']:
    path = os.path.join(src, f'fromGTF.{ev}.txt')
    lines = open(path).read().splitlines(keepends=True)
    # parse only to decide; rows are written back verbatim (rMATS keeps quotes, e.g.
    # GeneID "MSTRG.1", and needs them unchanged to match its own records)
    df = pd.read_csv(path, sep='\t', quoting=csv.QUOTE_NONE, dtype=str, keep_default_na=False)
    num = df.copy()
    for c in num.columns:
        if c not in ('GeneID', 'geneSymbol', 'chr', 'strand'):
            num[c] = pd.to_numeric(num[c])
    keep = rules.passes(ev, num)
    with open(os.path.join(dst, f'fromGTF.{ev}.txt'), 'w') as o:
        o.write(lines[0])
        for line, k in zip(lines[1:], keep):
            if k: o.write(line)
    print(f'{ev}\t{len(df)}\t{sum(keep)}')
