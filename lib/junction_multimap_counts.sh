#!/usr/bin/env bash
# Per-sample junction counts from multi-mapped reads, from one BAM (paralog filter, step 03).
#
# usage: junction_multimap_counts.sh <bam> <out.tsv> [samtools]
# out: chrom <tab> first intron base <tab> last intron base <tab> multi-mapped spliced reads
# Reads: PRIMARY alignments (-F 0x904) with MAPQ < 255 (STAR: 3 two loci, 1 three-four, 0 five+),
# so each multi-mapped read is counted once, at STAR's chosen primary locus. Together with the
# unique counts (support/junction_counts, MAPQ 255) this gives the multi-mapped fraction of a
# junction's reads: high values mark junctions in near-identical paralogs, where even "unique"
# placements are unreliable.
set -euo pipefail
bam=$1; out=$2; samtools=${3:-samtools}
"$samtools" view -F 0x904 "$bam" 2>/dev/null | awk -F'\t' -v OFS='\t' '
$5 < 255 && $6 ~ /N/ {
  p = $4; cig = $6
  while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
    l = substr(cig, RSTART, RLENGTH - 1) + 0; op = substr(cig, RSTART + RLENGTH - 1, 1); cig = substr(cig, RSTART + RLENGTH)
    if (op == "N") j[$3 OFS p OFS p + l - 1]++
    if (op ~ /[MDN=X]/) p += l
  }
}
END { for (k in j) print k, j[k] }' > "$out.tmp"
mv "$out.tmp" "$out"
