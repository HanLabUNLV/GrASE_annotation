"""Transcript-level read-through filter (step 03, stage D2).

usage: python readthrough_tx_filter.py <in.gtf> <out.gtf> <ref.gtf> <union.tracking> <gtf_list> <conditions.tsv>
                                       [max_ratio=0.5] [elev_frac=0.5] [elev_fold=3] [report.tsv] [compare=max]
The junction-level filter (stage D) only sees junctions whose two ends lie in different reference
genes. Read-through also arrives as an exon running across the gap between two genes, or as a
junction into an unannotated exon just outside a gene model; exon-overlap gene grouping then
fuses everything such a transcript touches. Here a transcript is BRIDGING when its exons overlap
exons of >= 2 reference genes on the same strand, however they are linked.
Support of a transcript = num_samples (samples whose assembly has exactly its intron chain).
A bridging transcript is dropped when every reference gene it touches also has non-bridging
transcripts and its support is < max_ratio x the best non-bridging support of the STRONGER of
those genes (compare=max, default since 2026-11-13; compare=min uses the weaker gene, i.e. the
bridge must be minor next to every gene). Kept otherwise: bridging transcripts that are the
dominant form (a gene the reference split in two, or a real fusion), and stress-elevated ones:
in some taxon, the stress group has >= elev_frac of its samples assembling the transcript and
>= elev_fold x the control group's fraction (taxon/condition from conditions.tsv; samples from
the union's tracking file, q1..qN in gtf_list order).
"""
import sys, re, bisect, collections

src, out, ref, tracking, gtflist, conds = sys.argv[1:7]
R = float(sys.argv[7]) if len(sys.argv) > 7 else 0.5
EF = float(sys.argv[8]) if len(sys.argv) > 8 else 0.5
EX = float(sys.argv[9]) if len(sys.argv) > 9 else 3.0
report = sys.argv[10] if len(sys.argv) > 10 and sys.argv[10] else None
cmp_ = (max if (sys.argv[11] if len(sys.argv) > 11 else 'max') == 'max' else min)
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]

ex, info, sup = collections.defaultdict(list), {}, {}
for l in open(src):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'transcript':
        m = re.search(r'num_samples "(\d+)"', f[8]); sup[tid(f[8])] = int(m[1]) if m else 1
    elif f[2] == 'exon':
        t = tid(f[8]); ex[t].append((int(f[3]), int(f[4]))); info[t] = (f[0], f[6])

rx = collections.defaultdict(list)
for l in open(ref):
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon':
        rx[(f[0], f[6])].append((int(f[3]), int(f[4]), re.search(r'gene_id "([^"]+)"', f[8])[1]))
idx = {k: ([a for a, _, _ in sorted(v)], sorted(v), max(b - a for a, b, _ in v)) for k, v in rx.items()}
def refgenes(t):
    h = set()
    if info[t] not in idx: return h
    st, v, mx = idx[info[t]]
    for a, b in ex[t]:
        i = bisect.bisect_left(st, a - mx)
        while i < len(v) and v[i][0] <= b:
            if v[i][1] >= a: h.add(v[i][2])
            i += 1
    return h
rg = {t: refgenes(t) for t in ex}
best = collections.defaultdict(int)
for t, g in rg.items():
    if len(g) == 1: best[next(iter(g))] = max(best[next(iter(g))], sup[t])
bridging = [t for t, g in rg.items() if len(g) >= 2]

# per-sample presence from the tracking file (only needed for bridging transcripts)
samples = [re.sub(r'\.scallop2\.gtf$', '', l.strip().split('/')[-1]) for l in open(gtflist) if l.strip()]
grp = {}
for l in open(conds):
    s, tx_, c = l.rstrip('\n').split('\t')[:3]; grp[s] = (tx_, c)
size = collections.Counter(grp[s] for s in samples)
need = set(bridging); present = {}
for l in open(tracking):
    f = l.rstrip('\n').split('\t'); t = f[0].split('|')[0]   # "TCONS_00000001|<exons>|<length>"
    if t in need:
        if len(f) - 4 != len(samples): sys.exit(f'tracking has {len(f) - 4} sample columns, gtf_list {len(samples)}')
        present[t] = [samples[i] for i, x in enumerate(f[4:]) if x != '-']

def elevated(t):
    n = collections.Counter(grp[s] for s in present.get(t, []))
    hits = []
    for (taxon, c), k in size.items():
        if c == 'control': continue
        fs = n[(taxon, c)] / k; fc = n[(taxon, 'control')] / size[(taxon, 'control')]
        if fs >= EF and fs >= EX * fc: hits.append(f'{taxon}:{c}')
    return hits

drop, kept_elev, rows = set(), 0, []
for t in bridging:
    g = rg[t]; b = [best[x] for x in g]
    minor = all(v > 0 for v in b) and sup[t] < R * cmp_(b)
    el = elevated(t) if minor else []
    decision = 'drop' if minor and not el else ('keep_stress' if el else ('keep_dominant' if all(v > 0 for v in b) else 'keep_no_alternative'))
    if decision == 'drop': drop.add(t)
    kept_elev += decision == 'keep_stress'
    rows.append((t, info[t][0], info[t][1], min(a for a, _ in ex[t]), max(b2 for _, b2 in ex[t]), sup[t],
                 ','.join(f'{x}:{best[x]}' for x in sorted(g)), ','.join(el), decision))
with open(out, 'w') as o:
    for l in open(src):
        if l[0] == '#' or tid(l.split('\t')[8]) not in drop: o.write(l)
if report:
    with open(report, 'w') as o:
        o.write('transcript_id\tchrom\tstrand\tstart\tend\tnum_samples\tref_genes:best_own_support\tstress_elevated_in\tdecision\n')
        for r in rows: o.write('\t'.join(map(str, r)) + '\n')
c = collections.Counter(r[-1] for r in rows)
print(f'bridging transcripts (exons on >= 2 reference genes): {len(bridging)} of {len(ex)}; '
      f'dropped {len(drop)} (support < {R} x the {cmp_.__name__}imum best own-transcript support of the touched genes); kept: '
      f'{c["keep_dominant"]} dominant/comparable, {c["keep_no_alternative"]} where a gene has no own transcript, '
      f'{kept_elev} stress-elevated')
