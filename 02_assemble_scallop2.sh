#!/usr/bin/env bash
# Step 02: per-sample Scallop2 assemblies on uniquely mapped reads (input to step 03).
#   Scallop2 1.1.2 (bioconda build h503566f_7; newer builds need AVX2 and crash on CPUs without it):
#     --library_type $LIBTYPE                     first = fr-firststrand (dUTP), second, unstranded
#     --min_mapping_quality $MIN_MAPQ             255 = unique for STAR (MAPQ 255 unique, 3 two loci,
#                                                 1 three-four, 0 five+); 60 for HISAT2. Applied while
#                                                 reading, so no filtered BAM copies are needed.
#     --min_transcript_coverage $MIN_TX_COV (3.0; Scallop2 default 0.5): drops faint
#                                                 intergenic/read-through transcripts
#   -> $OUT/scallop2/<group>/<sample>.scallop2.gtf, $OUT/gtf_list.scallop2.txt
# Why Scallop2 per sample and not a merged assembly: multi-mapped reads joined paralogs, and mergers
# (stringtie --merge, TACO) kept only 45-57% of well-supported junctions vs 73% in the per-sample
# assemblies; step 03 builds the annotation as a filtered union instead (pipeline/README.md).
# COMPARE_MERGERS=1 also writes those mergers' output for comparison (not used downstream):
#   stringtie --merge -g 0 -c 3 -f 0.1 -> merge/stringtie_merge.gtf; taco_run --isoform-frac 0.1 -> merge/taco/
# Settings: pipeline/config.sh. Restartable: finished samples are skipped.
#
# Usage: [OUT=assembly_unique MIN_TX_COV=3.0 COMPARE_MERGERS=0] pipeline/02_assemble_scallop2.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
OUT=$(readlink -m "${OUT:-assembly_unique}")
MIN_TX_COV=${MIN_TX_COV:-3.0}
COMPARE_MERGERS=${COMPARE_MERGERS:-0}
BAM_DIR=$(readlink -f "$BAM_DIR")

mkdir -p "$OUT"
awk -F'\t' '{print $2, $1}' "$SAMPLES" | sort -k1,1 -k2,2V > "$OUT/samples.txt"   # group sample
echo "== $(wc -l < "$OUT/samples.txt") samples -> $OUT ($(date))"
echo "   Scallop2: MIN_MAPQ=$MIN_MAPQ MIN_TX_COV=$MIN_TX_COV LIBTYPE=$LIBTYPE"

run_scallop2() {
  local grp=$1 s=$2 out=$OUT/scallop2/$1/$2.scallop2.gtf
  [[ -s $out ]] && { echo "skip $s"; return; }
  mkdir -p "$(dirname "$out")"
  "$SCALLOP2" -i "$BAM_DIR/$s$BAM_SUFFIX" -o "$out.tmp" --library_type "$LIBTYPE" \
    --min_mapping_quality "$MIN_MAPQ" --min_transcript_coverage "$MIN_TX_COV" --verbose 0 \
    > "$OUT/scallop2/$grp/$s.log" 2>&1 && mv "$out.tmp" "$out" && echo "ok   $s" || { echo "FAIL $s"; return 1; }
}
export -f run_scallop2; export OUT MIN_TX_COV
echo "== Scallop2 ($(date))"
xargs -a "$OUT/samples.txt" -P "$JOBS" -n 2 bash -c 'run_scallop2 "$0" "$1"' \
  | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
awk -v d="$OUT/scallop2" '{print d"/"$1"/"$2".scallop2.gtf"}' "$OUT/samples.txt" > "$OUT/gtf_list.scallop2.txt"
while read -r f; do [[ -s $f ]] || { echo "missing $f"; exit 1; }; done < "$OUT/gtf_list.scallop2.txt"
echo "   transcripts per sample: $(while read -r f; do awk '$3=="transcript"' "$f" | wc -l; done < "$OUT/gtf_list.scallop2.txt" | sort -n | awk '{a[NR]=$1} END{print "min "a[1]", median "a[int((NR+1)/2)]", max "a[NR]}')"

if [[ $COMPARE_MERGERS == 1 ]]; then
  mkdir -p "$OUT/merge"; M=$OUT/merge/stringtie_merge.gtf
  [[ -s $M ]] || "$STRINGTIE" --merge -g 0 -c 3 -f 0.1 -o "$M" "$OUT/gtf_list.scallop2.txt"
  [[ -s $OUT/merge/taco/assembly.gtf ]] || { rm -rf "$OUT/merge/taco"; "$TACO_RUN" -o "$OUT/merge/taco" -p "$THREADS" \
    --gtf-expr-attr cov --filter-min-expr 2.0 --isoform-frac 0.1 "$OUT/gtf_list.scallop2.txt" > "$OUT/merge/taco.log" 2>&1; }
  echo "   comparison merges in $OUT/merge/"
fi
echo "== done ($(date))"
