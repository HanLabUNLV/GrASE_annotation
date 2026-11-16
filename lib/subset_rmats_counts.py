"""Copy an rMATS output's event and count tables, restricted to a kept event set,
so that `rmats.py --task stat` tests only those events.

usage: python subset_rmats_counts.py <rmats_out_dir> <kept_events_dir> <stat_dir>
  <kept_events_dir>/fromGTF.<event>.txt  (from prefilter_rmats_events.py)
  writes <stat_dir>/fromGTF.<event>.txt, JC.raw.input.<event>.txt, JCEC.raw.input.<event>.txt
  (lines copied verbatim; counts as produced by the original post step)
"""
import sys, os

src, kept_dir, dst = sys.argv[1:4]
os.makedirs(dst, exist_ok=True)
for ev in ['SE', 'MXE', 'A3SS', 'A5SS', 'RI']:
    kept = open(os.path.join(kept_dir, f'fromGTF.{ev}.txt')).read().splitlines(keepends=True)
    ids = {l.split('\t', 1)[0] for l in kept[1:]}
    with open(os.path.join(dst, f'fromGTF.{ev}.txt'), 'w') as o: o.writelines(kept)
    n = {}
    for t in ('JC', 'JCEC'):
        lines = open(os.path.join(src, f'{t}.raw.input.{ev}.txt')).read().splitlines(keepends=True)
        sel = [l for l in lines[1:] if l.split('\t', 1)[0] in ids]
        with open(os.path.join(dst, f'{t}.raw.input.{ev}.txt'), 'w') as o: o.writelines(lines[:1] + sel)
        n[t] = len(sel)
    print(f"{ev}\tevents {len(ids)}\tJC rows {n['JC']}\tJCEC rows {n['JCEC']}")
