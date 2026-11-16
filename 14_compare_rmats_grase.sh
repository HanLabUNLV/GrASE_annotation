#!/usr/bin/env bash
# Step 14: GrASE vs rMATS, per group and contrast (GrASE/scripts/GSE188847 stage 80, without DEXSeq).
#   1. map_rmats_bipartitions.R once per class (internal, terminal), with all contrasts in one call: each
#      rMATS event is placed on the splicing graph and matched to the GrASE bipartition whose two sides
#      hold its two forms. The script maps the distinct events of all contrasts once (the contrasts' event
#      lists overlap heavily: 17% of event rows were distinct in the water study) and writes each
#      contrast's own mapping -> <contrast>/map.<class>/
#   2. compare_rmats_grase.R per contrast and class: events matched by both tools, called by either,
#      direction of change; exonic-part summary -> <contrast>/grase_vs_rmats.<class>.txt (+ .parts.txt,
#      discrepancies)
#   3. compare_parts_upset.R per contrast (GrASE internal + terminal vs rMATS on exonic parts, no DEXSeq)
#      -> <contrast>/upset/
# Calls: GrASE `significant` rule; rMATS FDR < 0.05 and |IncLevelDifference| >= 0.1 (both at --pval 0.05).
# rMATS was run with b1 = treatment, so --rmats_sample1 is the treatment group.
# Restartable: finished parts are skipped.
#
# Usage: [RMATS=rmats_prefiltered WD=grase OUT=comparison_rmats_grase THREADS=16] pipeline/14_compare_rmats_grase.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
RMATS=$(readlink -f "${RMATS:-rmats_prefiltered}")
WD=$(readlink -f "${WD:-grase}")
OUT=$(readlink -m "${OUT:-comparison_rmats_grase}")
CLASSES=${CLASSES:-internal terminal}
export PATH=$GRASE_ENV/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
S=$(readlink -f "$GRASE_DIR")/Rpkg/scripts
echo "== GrASE ($WD) vs rMATS ($RMATS) -> $OUT ($(date))"

contrast_dirs=()
for g in $(pipeline_groups); do for c in $(pipeline_treatments "$g"); do
  contrast_dirs+=("$g/${c}_vs_$CONTROL")
done; done

# 1. map all contrasts once per class
mkdir -p "$OUT"
for a in $CLASSES; do
  todo=0; for d in "${contrast_dirs[@]}"; do [[ -s $OUT/$d/map.$a/combined.fromGTF.SE.txt ]] || todo=1; done
  if (( todo )); then
    echo "-- mapping $a ($(date))"
    rdirs=$(printf "$RMATS/%s/out," "${contrast_dirs[@]}"); odirs=$(printf "$OUT/%s/map.$a," "${contrast_dirs[@]}")
    (cd "$OUT" && Rscript "$S/map_rmats_bipartitions.R" --graphdir "$WD/graphml" --splitdir "$WD/bipartition.filtered.anc12" \
      --class "$a" --rmatsdir "${rdirs%,}" --outdir "${odirs%,}" --mc_cores "$THREADS" > "$OUT/map.$a.log" 2>&1) \
      || { echo "mapping $a FAILED (see $OUT/map.$a.log)"; exit 1; }
  fi
  grep -E 'distinct|^==|events .*matched' "$OUT/map.$a.log" | sed "s/^/   $a /"
done

for d in "${contrast_dirs[@]}"; do
  g=${d%%/*}; c=${d#*/}; c=${c%_vs_*}; o=$OUT/$d; mkdir -p "$o"
  echo "-- $g $c ($(date))"
  # 2. event-level comparison per class
  for a in $CLASSES; do
    (cd "$o" && Rscript "$S/compare_rmats_grase.R" \
      --grase "$WD/test_anc12/$g/grase_exsj/test_bipartition.exsj.${a}_betabinom_EBmap.annotated.txt" \
      --rmats_dir "$RMATS/$d/out" --map_dir "map.$a/" --contrast "${g}_${c}_vs_${g}_$CONTROL" \
      --rmats_sample1 "${g}_$c" --pval 0.05 --output "grase_vs_rmats.$a.txt" > "compare.$a.log" 2>&1) \
      || { echo "   compare $a FAILED (see $o/compare.$a.log)"; continue; }
  done
  # 3. UpSet plots on exonic parts, GrASE (both classes) vs rMATS
  grase_files=$(for a in $CLASSES; do printf '%s,' "$WD/test_anc12/$g/grase_exsj/test_bipartition.exsj.${a}_betabinom_EBmap.annotated.txt"; done)
  map_dirs=$(for a in $CLASSES; do printf '%s,' "$o/map.$a"; done)
  (cd "$o" && Rscript "$S/compare_parts_upset.R" --grase "${grase_files%,}" --map_dirs "${map_dirs%,}" \
    --class "$(echo $CLASSES | tr ' ' ',')" --rmats_dir "$RMATS/$d/out" --contrast "${g}_${c}_vs_${g}_$CONTROL" \
    --pval 0.05 --outdir upset > upset.log 2>&1) || echo "   UpSet FAILED (see $o/upset.log)"
done
echo "== done ($(date))"
