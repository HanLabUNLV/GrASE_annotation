"""Map rMATS events onto GrASE bipartitions once for all contrasts (step 14).

rMATS detects events separately in every contrast, but the event lists overlap heavily (in the water
study only 17% of the 1.17 M event rows of 8 contrasts were distinct), and GrASE's matching of an event
(grase:::match_rmats_bipartitions) depends only on the event's own coordinates and its gene's graph and
bipartitions. So the distinct events are mapped once and the result is handed to each contrast.

usage:
  python rmats_event_union.py union <union_dir> <contrast_out_dir> [<contrast_out_dir> ...]
      writes <union_dir>/fromGTF.<type>.txt: every distinct event (by all coordinate columns), with a new ID
  python rmats_event_union.py split <union_map_dir> <union_dir> <contrast_out_dir> <contrast_map_dir>
      writes <contrast_map_dir>/combined.fromGTF.<type>.txt and combined.dexseq.<type>.mapped.txt as
      map_rmats_bipartitions.R would have written them for that contrast's own events and IDs
Types: SE A5SS A3SS MXE RI (those present). Files are read and written as text, so quoting stays as rMATS wrote it.
"""
import sys, os, csv

TYPES = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']
ID_COLS = {'ID', 'GeneID', 'geneSymbol'}

def read(path):
    with open(path) as f:
        r = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        h = next(r)
        return h, [row for row in r]

def write(path, header, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar=None, escapechar=None, lineterminator='\n')
        w.writerow(header); w.writerows(rows)

def key_of(h, row):
    return tuple(v for c, v in zip(h, row) if c not in ID_COLS)

def union(udir, cdirs):
    os.makedirs(udir, exist_ok=True)
    for t in TYPES:
        fs = [os.path.join(d, f'fromGTF.{t}.txt') for d in cdirs]
        fs = [f for f in fs if os.path.exists(f)]
        if not fs: continue
        seen, out, header, total = {}, [], None, 0
        for f in fs:
            h, rows = read(f)
            if header is None: header = h
            elif h != header: sys.exit(f'{f}: columns differ from {fs[0]}')
            total += len(rows)
            for row in rows:
                k = key_of(h, row)
                if k not in seen:
                    seen[k] = len(out); new = list(row); new[h.index('ID')] = str(len(out)); out.append(new)
        write(os.path.join(udir, f'fromGTF.{t}.txt'), header, out)
        print(f'{t}: {total} events in {len(fs)} contrasts -> {len(out)} distinct')

def split(umap, udir, cdir, cmap):
    os.makedirs(cmap, exist_ok=True)
    genes = set()   # the contrast's genes with events of any type (map_rmats_bipartitions.R lists parts for these)
    for t in TYPES:
        cf = os.path.join(cdir, f'fromGTF.{t}.txt')
        if os.path.exists(cf):
            ch, crows = read(cf); genes |= {r[ch.index('GeneID')].strip('"') for r in crows}
    for t in TYPES:
        cf = os.path.join(cdir, f'fromGTF.{t}.txt'); uf = os.path.join(udir, f'fromGTF.{t}.txt')
        mf = os.path.join(umap, f'combined.fromGTF.{t}.txt'); df = os.path.join(umap, f'combined.dexseq.{t}.mapped.txt')
        if not (os.path.exists(cf) and os.path.exists(mf)): continue
        uh, urows = read(uf)
        ukey = {r[uh.index('ID')]: key_of(uh, r) for r in urows}
        mh, mrows = read(mf)                       # union events + mapping columns
        extra = [c for c in mh if c not in uh]     # DexseqFragment, DexseqRefFrag, bipartID, inclSide
        by_key = {ukey[r[mh.index('ID')]]: [r[mh.index(c)] for c in extra] for r in mrows}
        ch, crows = read(cf)
        u2c = {}                                    # union ID -> this contrast's ID
        kid = {key_of(ch, r): r[ch.index('ID')] for r in crows}
        for uid, k in ukey.items():
            if k in kid: u2c[uid] = kid[k]
        write(os.path.join(cmap, f'combined.fromGTF.{t}.txt'), ch + extra,
              [r + by_key.get(key_of(ch, r), ['NA'] * len(extra)) for r in crows])
        # per-part event lists: keep this contrast's events, under its own IDs
        dh, drows = read(df)
        pre = f'{t}_'
        def tr(cell):
            if cell in ('NA', ''): return cell
            ids = [pre + u2c[x[len(pre):]] for x in cell.split(',') if x[len(pre):] in u2c]
            return ','.join(ids) if ids else 'NA'
        icols = [i for i, c in enumerate(dh) if c.startswith('rMATS_ID_')]
        gi = dh.index('GeneID')
        out = []
        for r in drows:
            if r[gi].strip('"') not in genes: continue
            r = list(r)
            for i in icols: r[i] = tr(r[i])
            out.append(r)
        write(os.path.join(cmap, f'combined.dexseq.{t}.mapped.txt'), dh, out)

if __name__ == '__main__':
    if len(sys.argv) >= 4 and sys.argv[1] == 'union': union(sys.argv[2], sys.argv[3:])
    elif len(sys.argv) == 6 and sys.argv[1] == 'split': split(*sys.argv[2:])
    else: sys.exit(__doc__)
