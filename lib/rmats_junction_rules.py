"""Junction rules for rMATS events (shared by step 11 and step 12).

An event passes when every junction of every isoform (inclusion and skipping forms) is
canonical (GT-AG/GC-AG/AT-AC on the event strand), <= max_intron bp, and has >= min_reads
unique reads in >= min_samples samples (per-sample counts from the BAMs, step 01 cache).
These are the annotation's own junction rules (step 01), so events built only from
annotated junctions pass by construction.

Overhang (anchor) rule, when anchor_dir is given (support/junction_counts_anchor, from
lib/junction_anchor_counts.sh with min_anchor 12): of the >= min_samples samples with
>= min_reads reads, each must also have >= 1 read with >= 12 aligned bases on both sides of
the junction (STAR's SJ.out.tab overhang rule, outSJfilterOverhangMin 12). Junctions next to
an internal exon shorter than anchor_min (micro-exon; SE/MXE inner exon) are exempt, since
no read can have a 12 bp block inside such an exon. Removes junctions supported mainly by
untrimmed read ends forced into a splice (e.g. 90M65N3M7S).
"""
import os, glob
from collections import defaultdict

CANON = {'GT-AG', 'GC-AG', 'AT-AC'}
_COMP = str.maketrans('ACGT', 'TGCA')


class JunctionRules:
    def __init__(self, junc_dir, genome, max_intron=20000, min_samples=3, min_reads=3,
                 anchor_dir=None, anchor_min=12):
        self.max_intron, self.min_samples, self.anchor_min = max_intron, min_samples, anchor_min
        self.nsamp = defaultdict(int)
        self.nanch = None   # samples with >= min_reads reads and >= 1 anchored read
        if anchor_dir:
            self.nanch = defaultdict(int)
            fns = glob.glob(os.path.join(anchor_dir, '*.junc_anchor.tsv'))
            if not fns: raise FileNotFoundError(f'no *.junc_anchor.tsv in {anchor_dir}')
            for fn in fns:
                for l in open(fn):
                    c, s, e, n, a = l.split('\t')
                    if int(n) >= min_reads and int(a) >= 1: self.nanch[(c, int(s), int(e))] += 1
        for fn in glob.glob(os.path.join(junc_dir, '*.junc.tsv')):
            for l in open(fn):
                c, s, e, n = l.split('\t')
                if int(n) >= min_reads: self.nsamp[(c, int(s), int(e))] += 1
        self.seq, name, buf = {}, None, []
        for l in open(genome):
            if l[0] == '>':
                if name: self.seq[name] = ''.join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip().upper())
        self.seq[name] = ''.join(buf)

    def junction_ok(self, c, strand, s, e, anchor_exempt=False):
        """s, e: first and last intron base, 1-based."""
        if e < s or e - s + 1 > self.max_intron or self.nsamp[(c, s, e)] < self.min_samples: return False
        if self.nanch is not None and not anchor_exempt and self.nanch[(c, s, e)] < self.min_samples: return False
        m = self.seq[c][s-1:s+1] + '-' + self.seq[c][e-2:e]
        if strand == '-': m = m[::-1].translate(_COMP)
        return m in CANON

    @staticmethod
    def event_junctions(ev, r):
        """Junctions (first, last intron base; 1-based) of every isoform of an rMATS event row.
        rMATS: *ES / *Start_0base are 0-based exon starts, *EE / *End are 1-based exon ends."""
        g = r.__getitem__ if hasattr(r, '__getitem__') else lambda k: getattr(r, k)
        if ev == 'SE':
            return [(g('upstreamEE') + 1, g('exonStart_0base')), (g('exonEnd') + 1, g('downstreamES')),
                    (g('upstreamEE') + 1, g('downstreamES'))]
        if ev == 'MXE':
            return [(g('upstreamEE') + 1, g('1stExonStart_0base')), (g('1stExonEnd') + 1, g('downstreamES')),
                    (g('upstreamEE') + 1, g('2ndExonStart_0base')), (g('2ndExonEnd') + 1, g('downstreamES'))]
        if ev in ('A5SS', 'A3SS'):
            if g('flankingES') > g('longExonEnd'):   # flanking exon downstream
                return [(g('longExonEnd') + 1, g('flankingES')), (g('shortEE') + 1, g('flankingES'))]
            return [(g('flankingEE') + 1, g('longExonStart_0base')), (g('flankingEE') + 1, g('shortES'))]
        if ev == 'RI':
            return [(g('upstreamEE') + 1, g('downstreamES'))]
        raise ValueError(ev)

    def micro_exempt(self, ev, r):
        """Per junction of event_junctions(ev, r): True if it borders an inner exon < anchor_min bp."""
        g = r.__getitem__ if hasattr(r, '__getitem__') else lambda k: getattr(r, k)
        m = self.anchor_min
        if ev == 'SE':
            x = g('exonEnd') - g('exonStart_0base') < m
            return [x, x, False]
        if ev == 'MXE':
            x1 = g('1stExonEnd') - g('1stExonStart_0base') < m
            x2 = g('2ndExonEnd') - g('2ndExonStart_0base') < m
            return [x1, x1, x2, x2]
        return [False] * len(self.event_junctions(ev, r))

    def passes(self, ev, df):
        """Boolean list: does each row of an rMATS event table pass?"""
        return [all(self.junction_ok(r['chr'], r['strand'], s, e, x)
                    for (s, e), x in zip(self.event_junctions(ev, r), self.micro_exempt(ev, r)))
                for _, r in df.iterrows()]
