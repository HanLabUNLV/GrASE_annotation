#!/usr/bin/env bash
# Per-sample junction counts with anchor (overhang) information, from one BAM.
#
# usage: junction_anchor_counts.sh <bam> <out.tsv> [min_anchor=12] [samtools]
# out: chrom <tab> first intron base <tab> last intron base <tab> unique spliced reads
#      <tab> reads whose aligned bases on BOTH sides of this junction are >= min_anchor
# Reads: primary, not QC-fail/supplementary (-F 0x904), uniquely mapped (MAPQ 255), as in
# the step 01 cache. Anchor of a junction = aligned (M/=/X) bases in the alignment block
# immediately before / after that N operation (STAR's overhang notion).
set -euo pipefail
bam=$1; out=$2; amin=${3:-12}; samtools=${4:-samtools}
"$samtools" view -q 255 -F 0x904 "$bam" 2>/dev/null | awk -F'\t' -v OFS='\t' -v amin="$amin" '
$6 ~ /N/ {
  p = $4; cig = $6; nb = 0; blk = 0
  # split CIGAR into blocks separated by N; record matched bases per block and N positions
  while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
    l = substr(cig, RSTART, RLENGTH - 1) + 0; op = substr(cig, RSTART + RLENGTH - 1, 1); cig = substr(cig, RSTART + RLENGTH)
    if (op == "N") { nb++; ns[nb] = p; ne[nb] = p + l - 1; mb[nb] = blk; blk = 0 }
    else if (op ~ /[M=X]/) blk += l
    if (op ~ /[MDN=X]/) p += l
  }
  mb[nb + 1] = blk
  for (i = 1; i <= nb; i++) {
    k = $3 OFS ns[i] OFS ne[i]; t[k]++
    if (mb[i] >= amin && mb[i + 1] >= amin) a[k]++
  }
}
END { for (k in t) print k, t[k], a[k] + 0 }' > "$out.tmp"
mv "$out.tmp" "$out"
