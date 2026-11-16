"""Select transcripts from a junction-filtered union of per-sample assemblies (step 03).

usage: python select_transcripts.py <in.gtf> <out.gtf> <unique_junc_dir> <multi_junc_dir>
                                    [min_samples=5] [max_multi_frac=0.5]
in.gtf:  gffcompare -i union after the junction rules (step 01), with num_samples "<n>" per
         transcript (samples whose assembly has exactly this intron chain).
1. Paralog filter: drop transcripts with a junction whose multi-mapped share of reads,
   multi / (unique + multi) summed over all samples, exceeds max_multi_frac. Unique reads
   (MAPQ 255) from <unique_junc_dir>/*.junc.tsv, multi-mapped primary reads (MAPQ < 255)
   from <multi_junc_dir>/*.junc_multi.tsv (lib/junction_multimap_counts.sh).
2. Keep transcripts assembled in >= min_samples samples, then rescue rarer ones that carry a
   junction no kept transcript has: in order of decreasing num_samples, then most new
   junctions, a transcript is added if it adds a junction. No junction that survived the
   rules is lost; rare recombinations of known junctions (which multiply GrASE's paths) are.
Single-exon transcripts are dropped (no junctions; strand often unknown).
Writes the kept transcripts' lines unchanged; prints counts.
"""
import sys, re, glob, collections

src, out, udir, mdir = sys.argv[1:5]
K = int(sys.argv[5]) if len(sys.argv) > 5 else 5
maxf = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]

ex, sup, chrom = collections.defaultdict(list), {}, {}
for l in open(src):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'transcript':
        t = tid(f[8]); chrom[t] = f[0]
        m = re.search(r'num_samples "(\d+)"', f[8]); sup[t] = int(m[1]) if m else 1
    elif f[2] == 'exon':
        ex[tid(f[8])].append((int(f[3]), int(f[4])))
intr = {}
for t, e in ex.items():
    e.sort(); intr[t] = {(chrom[t], b + 1, a - 1) for (_, b), (a, _) in zip(e, e[1:])}
multi = [t for t in ex if intr[t]]
alljunc = set().union(*(intr[t] for t in multi))

def counts(d, pat):
    c = collections.Counter()
    fns = glob.glob(f'{d}/{pat}')
    if not fns: sys.exit(f'no {pat} in {d}')
    for fn in fns:
        for l in open(fn):
            x = l.split('\t'); k = (x[0], int(x[1]), int(x[2]))
            if k in alljunc: c[k] += int(x[3])
    return c
u, m = counts(udir, '*.junc.tsv'), counts(mdir, '*.junc_multi.tsv')
bad = {j for j in alljunc if m[j] > maxf * (u[j] + m[j])}
cand = [t for t in multi if not (intr[t] & bad)]
print(f'input: {len(ex)} transcripts ({len(ex) - len(multi)} single-exon dropped), {len(alljunc)} junctions')
print(f'paralog filter: {len(bad)} junctions with multi-mapped share > {maxf}; '
      f'{len(multi) - len(cand)} transcripts dropped, {len(cand)} left')

keep = {t for t in cand if sup[t] >= K}
J = set().union(*(intr[t] for t in keep)) if keep else set()
rescued = 0
for t in sorted((t for t in cand if sup[t] < K), key=lambda t: (-sup[t], -len(intr[t] - J))):
    if intr[t] - J:
        keep.add(t); J |= intr[t]; rescued += 1
candJ = set().union(*(intr[t] for t in cand))
print(f'support >= {K} samples: {len(keep) - rescued}; rescued for a new junction: {rescued}; '
      f'kept {len(keep)} transcripts, {len(J)} of {len(candJ)} junctions')
with open(out, 'w') as o:
    for l in open(src):
        if l[0] == '#' or tid(l.split('\t')[8]) in keep: o.write(l)
