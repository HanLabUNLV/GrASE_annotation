"""Protein-homology check of bridged reference genes (step 03, stage D3).

usage: python homology_bridges.py <in.gtf> <out.gtf> <ref.gtf> <genome.fa> <blast_db> <cache_dir>
                                  <report.tsv> [blastp] [threads=16] [max_ovl=0.3]
A transcript that survived stage D2 and touches >= 2 reference genes either joins one gene that
the reference split in two (keep), or two genes for unrelated proteins / two copies of a
duplicated gene (a fusion or chimera, not a splicing unit: drop). Expression cannot tell these
apart; protein homology can. For each pair of reference genes a transcript joins, the
reference proteins (first CDS model of each gene, translated) are searched against Swiss-Prot
(blastp, E <= 1e-5, top 20 hits) and the pair classified:
  split       a shared subject protein hit in different regions (overlap <= max_ovl of the
              shorter hit): halves of one protein -> keep
  duplicate   a shared subject hit in the same region: paralog copies -> drop
  different   both have hits, none shared: unrelated proteins -> drop
  unknown     no CDS or no hit on one or both sides -> keep (no evidence)
A transcript is dropped if any pair it joins is duplicate or different.
BLAST results are cached per gene in <cache_dir>/ref_proteins.blastp.tsv (+ .done), so reruns
and other annotations of the same reference only search new genes.
"""
import sys, re, os, subprocess, bisect, collections

(src, out, ref, genome, db, cache, report) = sys.argv[1:8]
BLASTP = sys.argv[8] if len(sys.argv) > 8 and sys.argv[8] else 'blastp'
THREADS = sys.argv[9] if len(sys.argv) > 9 else '16'
MAXOVL = float(sys.argv[10]) if len(sys.argv) > 10 else 0.3
tid = lambda s: re.search(r'transcript_id "([^"]+)"', s)[1]

# transcripts and the reference genes their exons touch
ex, info = collections.defaultdict(list), {}
for l in open(src):
    if l[0] == '#': continue
    f = l.split('\t')
    if f[2] == 'exon':
        t = tid(f[8]); ex[t].append((int(f[3]), int(f[4]))); info[t] = (f[0], f[6])
rx = collections.defaultdict(list); cds = collections.defaultdict(list); first_tx = {}
for l in open(ref):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    g = re.search(r'gene_id "([^"]+)"', f[8])
    if not g: continue
    g = g[1]
    if f[2] == 'exon': rx[(f[0], f[6])].append((int(f[3]), int(f[4]), g))
    elif f[2] == 'CDS':
        t = re.search(r'transcript_id "([^"]+)"', f[8])[1]; first_tx.setdefault(g, t)
        if t == first_tx[g]: cds[g].append((f[0], int(f[3]), int(f[4]), f[6]))
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
bridges = {t: sorted(g) for t, g in rg.items() if len(g) >= 2}
need = sorted({g for gs in bridges.values() for g in gs})

# BLAST cache
os.makedirs(cache, exist_ok=True)
ctsv, cdone = f'{cache}/ref_proteins.blastp.tsv', f'{cache}/ref_proteins.blastp.done'
done = set(open(cdone).read().split()) if os.path.exists(cdone) else set()
todo = [g for g in need if g not in done]
if todo:
    chroms = {cds[g][0][0] for g in todo if cds.get(g)}
    seq, name, buf = {}, None, []
    for l in open(genome):
        if l[0] == '>':
            if name in chroms: seq[name] = ''.join(buf)
            name, buf = l[1:].split()[0], []
        elif name in chroms: buf.append(l.strip().upper())
    if name in chroms: seq[name] = ''.join(buf)
    comp = str.maketrans('ACGTN', 'TGCAN'); b = 'TCAG'
    aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    code = {x + y + z: aa[16 * i + 4 * j + k] for i, x in enumerate(b) for j, y in enumerate(b) for k, z in enumerate(b)}
    qf = f'{cache}/query.tmp.fa'
    with open(qf, 'w') as o:
        for g in todo:
            if not cds.get(g): continue
            v = sorted(cds[g], key=lambda x: x[1]); c, st = v[0][0], v[0][3]
            s = ''.join(seq[c][a - 1:e] for _, a, e, _ in v)
            if st == '-': s = s.translate(comp)[::-1]
            p = ''.join(code.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3)).rstrip('*')
            if len(p) >= 30: o.write(f'>{g}\n{p}\n')
    print(f'blastp: {len(todo)} new genes ({len(need) - len(todo)} cached) against {db}', flush=True)
    r = subprocess.run([BLASTP, '-query', qf, '-db', db, '-evalue', '1e-5', '-max_target_seqs', '20',
                        '-num_threads', THREADS, '-outfmt',
                        '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen stitle'],
                       capture_output=True, text=True)
    if r.returncode: sys.exit(r.stderr)
    with open(ctsv, 'a') as o: o.write(r.stdout)
    with open(cdone, 'a') as o: o.write(''.join(g + '\n' for g in todo))
    os.remove(qf)

hits = collections.defaultdict(dict)     # gene -> subject -> (sstart, send, title)  (best HSP)
for l in open(ctsv) if os.path.exists(ctsv) else []:
    x = l.rstrip('\n').split('\t')
    if x[0] in hits and x[1] in hits[x[0]]: continue
    hits[x[0]][x[1]] = (int(x[6]), int(x[7]), x[11].split(' OS=')[0])

def classify(a, b):
    A, B = hits.get(a, {}), hits.get(b, {})
    if not A or not B: return 'unknown', ''
    common = [s for s in A if s in B]
    if not common: return 'different', f'{next(iter(A.values()))[2]} | {next(iter(B.values()))[2]}'
    s = common[0]; (a1, a2, d), (b1, b2, _) = A[s], B[s]
    ovl = max(0, min(a2, b2) - max(a1, b1) + 1); short = min(a2 - a1, b2 - b1) + 1
    return ('split' if ovl <= MAXOVL * short else 'duplicate'), f'{d}: {a1}-{a2} / {b1}-{b2}'

pair_cls, drop = {}, set()
for t, gs in bridges.items():
    bad = False
    for i in range(len(gs)):
        for j in range(i + 1, len(gs)):
            k = (gs[i], gs[j])
            if k not in pair_cls: pair_cls[k] = classify(*k)
            bad |= pair_cls[k][0] in ('duplicate', 'different')
    if bad: drop.add(t)
with open(report, 'w') as o:
    o.write('gene_a\tgene_b\tclass\tevidence\n')
    for (a, b), (c, e) in sorted(pair_cls.items()): o.write(f'{a}\t{b}\t{c}\t{e}\n')
with open(out, 'w') as o:
    for l in open(src):
        if l[0] == '#' or tid(l.split('\t')[8]) not in drop: o.write(l)
c = collections.Counter(v[0] for v in pair_cls.values())
print(f'bridging transcripts: {len(bridges)}; reference-gene pairs: {len(pair_cls)} '
      f'(split {c["split"]}, different {c["different"]}, duplicate {c["duplicate"]}, unknown {c["unknown"]}); '
      f'dropped {len(drop)} transcripts joining different or duplicated genes')
