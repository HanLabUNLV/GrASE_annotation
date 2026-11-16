import sys
fa, tsv, out = sys.argv[1:4]
want = {}
for l in open(tsv):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'chrom': continue
    want.setdefault(f[0], []).append(f)
seq = {}; name = None; buf = []
def flush():
    if name in want: seq[name] = ''.join(buf)
for l in open(fa):
    if l[0] == '>':
        flush(); name = l[1:].split()[0]; buf = []
    elif name in want: buf.append(l.strip().upper())
flush()
comp = str.maketrans('ACGT', 'TGCA')
with open(out, 'w') as o:
    for c, rows in want.items():
        s0 = seq[c]
        for f in rows:
            s, e = int(f[2]), int(f[3])
            m = s0[s-1:s+1] + '-' + s0[e-2:e]
            if f[1] == '-': m = m[::-1].translate(comp)
            o.write('\t'.join(f + [m]) + '\n')
