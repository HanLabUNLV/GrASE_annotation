#!/usr/bin/env bash
# Step 10: GrASE Stage 1 counting. DEXSeq exonic-part counts per BAM against the
# concatenated per-gene GFF from step 09 (counts must use that GFF to join the graphs).
#   --paired yes --stranded reverse (fr-firststrand / dUTP) --order pos -a 10
#   (-a 10 skips multi-mapped reads: STAR gives them MAPQ <= 3)
# Count files: $WD/DEXSeq/count_files/<taxon>_<condition>/<S>_counts.txt, condition from
# $SAMPLES (pipeline/config.sh). One directory per taxon x condition, so the
# control directory serves both contrasts of a taxon (--cond1/--cond2 in later stages).
# Restartable: finished samples are skipped.
#
# Usage: WD=<grase workdir> pipeline/10_grase_count.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
WD=$(readlink -m "${WD:-grase}")
BAM_DIR=$(readlink -f "$BAM_DIR")
export PATH=$GRASE_ENV/bin:$PATH

GFF=$(ls "$WD"/ref/*.dexseq.bygene.gff)
DEXSEQ_COUNT=$(Rscript -e 'cat(system.file("python","dexseq_count.py",package="grase"))')
[[ -s $GFF && -s $DEXSEQ_COUNT ]] || { echo "missing $GFF or dexseq_count.py"; exit 1; }

count_one() {   # sample taxon condition
  local s=$1 d=$WD/DEXSeq/count_files/$2_$3 out
  out=$d/${s}_counts.txt; mkdir -p "$d"
  [[ -s $out ]] && { echo "skip $s"; return; }
  python "$DEXSEQ_COUNT" --format bam --order pos --paired "$PAIRED" --stranded "$DEXSEQ_STRANDED" -a 10 \
    "$GFF" "$BAM_DIR/$s$BAM_SUFFIX" "$out.tmp" 2> "$d/${s}_counts.log" \
    && mv "$out.tmp" "$out" && echo "ok   $s" || { echo "FAIL $s (see $d/${s}_counts.log)"; return 1; }
}
export -f count_one; export WD GFF DEXSEQ_COUNT BAM_DIR BAM_SUFFIX PAIRED DEXSEQ_STRANDED
echo "== DEXSeq counts ($(date)): $GFF"
xargs -a "$SAMPLES" -P "$JOBS" -n 3 bash -c 'count_one "$0" "$1" "$2"' \
  | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
echo "== done ($(date))"
