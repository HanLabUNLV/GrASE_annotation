"""Step 11: filter rMATS JC results (tag jc10_fdr0.05_absdpsi0.1), as documented in
taxon*/README.md for the original ../filter_rmats_jc_fdr_dpsi.py (not on this server):
  1. mean (IJC + SJC) per sample >= MIN_JC in at least one group
  2. FDR < MAX_FDR
  3. |IncLevelDifference| >= MIN_DPSI
Adds ref_gene_id (MAKER-style id from the GTF's ref_gene_id attribute, "novel_locus" if none).

Junction QC (when --junc_dir is given): every junction of the event (inclusion and
skipping forms) must pass the annotation's junction rules (step 01): canonical motif
(GT-AG/GC-AG/AT-AC on the event strand), length <= --max_intron, and >= --min_reads unique
reads in >= --min_samples samples (counted from the BAMs). Annotated junctions pass by
construction; this removes --novelSS events whose novel junctions fail the rules.
Overhang rule (default with --junc_dir since 2026-11-12; --no_anchor turns it off): those
samples must also each have >= 1 read with >= 12 bp aligned on both sides of the junction
(--anchor_dir, from lib/junction_anchor_counts.sh); SE/MXE inner exons < 12 bp are exempt.
Output then goes to filtered_<tag>_jqc_anc12/ (without it: filtered_<tag>_jqc/).

usage: python 11_filter_rmats.py <rmats_base_dir> <annotation.gtf> [min_jc=10] [max_fdr=0.05] [min_dpsi=0.1]
         [--junc_dir support/junction_counts --genome <genome.fa> --max_intron 20000 --min_samples 3 --min_reads 3
          --anchor_dir support/junction_counts_anchor --anchor_min 12 | --no_anchor]
  e.g.   python pipeline/11_filter_rmats.py rmats_scallop2_prefiltered assembly_unique/annotation.final.gtf \\
           --junc_dir support/junction_counts --genome <genome.fa>
  reads    <base>/<group>/<contrast>/out/{SE,A5SS,A3SS,MXE,RI}.MATS.JC.txt   (any <group>, any <a>_vs_<b>)
  writes   <base>/<taxon>/<contrast>/filtered_<tag>/<event>.MATS.JC.filtered.csv        (steps 1-3)
           <base>/<taxon>/<contrast>/filtered_<tag>_jqc[_anc12]/<event>.MATS.JC.filtered.csv (+ junction QC)
           <base>/summary_counts_<tag>[_anc12].tsv
"""
import sys, os, re, glob, argparse
from collections import defaultdict
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('base'); ap.add_argument('gtf')
ap.add_argument('min_jc', nargs='?', type=float, default=10)
ap.add_argument('max_fdr', nargs='?', type=float, default=0.05)
ap.add_argument('min_dpsi', nargs='?', type=float, default=0.1)
ap.add_argument('--junc_dir'); ap.add_argument('--genome')
ap.add_argument('--max_intron', type=int, default=20000)
ap.add_argument('--min_samples', type=int, default=3)
ap.add_argument('--min_reads', type=int, default=3)
ap.add_argument('--anchor_dir', default='support/junction_counts_anchor')
ap.add_argument('--anchor_min', type=int, default=12)
ap.add_argument('--no_anchor', action='store_true')
a = ap.parse_args()
base, gtf, min_jc, max_fdr, min_dpsi = a.base, a.gtf, a.min_jc, a.max_fdr, a.min_dpsi
tag = f"jc{min_jc:g}_fdr{max_fdr:g}_absdpsi{min_dpsi:g}"
EVENTS = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

refmap = defaultdict(set)
for l in open(gtf):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] != 'transcript': continue
    g = re.search(r'gene_id "([^"]+)"', f[8])[1]
    r = re.search(r'ref_gene_id "([^"]+)"', f[8])
    if r: refmap[g].add(r[1])

def mean_counts(col):
    def m(s):
        v = [float(x) for x in str(s).split(',') if x not in ('', 'NA')]
        return sum(v) / len(v) if v else 0.0
    return col.map(m)

# ---- junction QC (rules shared with step 12: pipeline/lib/rmats_junction_rules.py) ----
jqc = a.junc_dir is not None
if jqc:
    if not a.genome: sys.exit('--genome is required with --junc_dir')
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
    from rmats_junction_rules import JunctionRules
    anchor = None if a.no_anchor else a.anchor_dir
    rules = JunctionRules(a.junc_dir, a.genome, a.max_intron, a.min_samples, a.min_reads, anchor, a.anchor_min)
    passes = rules.passes
    qtag = f'_jqc_anc{a.anchor_min}' if anchor else '_jqc'
    print(f'junction QC: {a.junc_dir}' + (f', overhang rule from {anchor} (exempt inner exons < {a.anchor_min} bp)' if anchor else ''))

# ---- filter ------------------------------------------------------------------------
rows = []
for d in sorted(glob.glob(os.path.join(base, '*', '*_vs_*'))):
    if not os.path.isdir(os.path.join(d, 'out')): continue
    taxon, contrast = d.rstrip('/').split(os.sep)[-2:]
    out = os.path.join(d, f'filtered_{tag}'); os.makedirs(out, exist_ok=True)
    if jqc: outq = os.path.join(d, f'filtered_{tag}{qtag}'); os.makedirs(outq, exist_ok=True)
    for ev in EVENTS:
        p = os.path.join(d, 'out', f'{ev}.MATS.JC.txt')
        if not os.path.exists(p): continue
        df = pd.read_csv(p, sep='\t', dtype={'GeneID': str})
        jc1 = mean_counts(df.IJC_SAMPLE_1) + mean_counts(df.SJC_SAMPLE_1)
        jc2 = mean_counts(df.IJC_SAMPLE_2) + mean_counts(df.SJC_SAMPLE_2)
        keep = ((jc1 >= min_jc) | (jc2 >= min_jc)) & (df.FDR < max_fdr) & (df.IncLevelDifference.abs() >= min_dpsi)
        f = df[keep].copy()
        f['ref_gene_id'] = f.GeneID.map(lambda g: ','.join(sorted(refmap[g])) if refmap.get(g) else 'novel_locus')
        f.to_csv(os.path.join(out, f'{ev}.MATS.JC.filtered.csv'), index=False)
        row = [taxon, contrast, ev, len(df), int(((jc1 >= min_jc) | (jc2 >= min_jc)).sum()), len(f)]
        if jqc:
            fq = f[passes(ev, f)] if len(f) else f
            fq.to_csv(os.path.join(outq, f'{ev}.MATS.JC.filtered.csv'), index=False)
            row.append(len(fq))
        rows.append(row)

cols = ['taxon', 'contrast', 'event', 'tested', 'pass_jc', 'significant'] + (['significant_junction_qc'] if jqc else [])
s = pd.DataFrame(rows, columns=cols)
s.to_csv(os.path.join(base, f'summary_counts_{tag}{qtag[4:] if jqc else ""}.tsv'), sep='\t', index=False)
for v in cols[5:]:
    print(f'== {v}'); print(s.pivot_table(index=['taxon', 'contrast'], columns='event', values=v, aggfunc='sum').to_string())
