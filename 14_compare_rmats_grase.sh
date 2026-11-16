#!/usr/bin/env bash
# Step 14: GrASE vs rMATS, per group and contrast (GrASE/scripts/GSE188847 stage 80, without DEXSeq).
#   1. lib/rmats_event_union.py union: the distinct rMATS events of all contrasts (rMATS detects events
#      per contrast, but the lists overlap: 17% of event rows were distinct in the water study)
#   2. map_rmats_bipartitions.R once per class (internal, terminal) on that union: each rMATS event is
#      placed on the splicing graph and matched to the GrASE bipartition whose two sides hold its two
#      forms. The match depends only on the event and its gene, so mapping the union once is the same
#      as mapping every contrast (and 1/8 of the work here).
#   3. lib/rmats_event_union.py split: each contrast's own events and IDs -> <contrast>/map.<class>/
#   4. compare_rmats_grase.R per contrast and class: events matched by both tools, called by either,
#      direction of change; exonic-part summary -> <contrast>/grase_vs_rmats.<class>.txt (+ .parts.txt,
#      discrepancies)
#   5. compare_parts_upset.R per contrast (GrASE internal + terminal vs rMATS on exonic parts, no DEXSeq)
#      -> <contrast>/upset/
# Calls: GrASE `significant` rule; rMATS FDR < 0.05 and |IncLevelDifference| >= 0.1 (both at --pval 0.05).
# rMATS was run with b1 = treatment, so --rmats_sample1 is the treatment group.
# Restartable: finished parts are skipped.
#
# Usage: [RMATS=rmats_prefiltered WD=grase OUT=comparison_rmats_grase JOBS=16] pipeline/14_compare_rmats_grase.sh
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
U=$OUT/union
echo "== GrASE ($WD) vs rMATS ($RMATS) -> $OUT ($(date))"

contrast_dirs=()
for g in $(pipeline_groups); do for c in $(pipeline_treatments "$g"); do
  contrast_dirs+=("$g/${c}_vs_$CONTROL")
done; done

# 1. union of events
if [[ ! -s $U/events/fromGTF.SE.txt ]]; then
  "$PYTHON" -I pipeline/lib/rmats_event_union.py union "$U/events" \
    $(for d in "${contrast_dirs[@]}"; do echo "$RMATS/$d/out"; done) | sed 's/^/   /'
fi

# 2. map the union once per class
for a in $CLASSES; do
  if [[ ! -s $U/map.$a/combined.fromGTF.SE.txt ]]; then
    echo "-- mapping $a ($(date))"
    (cd "$U" && Rscript "$S/map_rmats_bipartitions.R" --graphdir "$WD/graphml" --splitdir "$WD/bipartition.filtered.anc12" \
      --class "$a" --rmatsdir "$U/events" --outdir "$U/map.$a" --mc_cores "$THREADS" > "$U/map.$a.log" 2>&1) \
      || { echo "mapping $a FAILED (see $U/map.$a.log)"; exit 1; }
  fi
  grep -E 'events .*matched' "$U/map.$a.log" | sed "s/^/   $a /"
done

for d in "${contrast_dirs[@]}"; do
  g=${d%%/*}; c=${d#*/}; c=${c%_vs_*}; o=$OUT/$d; mkdir -p "$o"
  echo "-- $g $c ($(date))"
  # 3. this contrast's events and IDs
  for a in $CLASSES; do
    [[ -s $o/map.$a/combined.fromGTF.SE.txt ]] || \
      "$PYTHON" -I pipeline/lib/rmats_event_union.py split "$U/map.$a" "$U/events" "$RMATS/$d/out" "$o/map.$a"
  done
  # 4. event-level comparison per class
  for a in $CLASSES; do
    (cd "$o" && Rscript "$S/compare_rmats_grase.R" \
      --grase "$WD/test_anc12/$g/grase_exsj/test_bipartition.exsj.${a}_betabinom_EBmap.annotated.txt" \
      --rmats_dir "$RMATS/$d/out" --map_dir "map.$a/" --contrast "${g}_${c}_vs_${g}_$CONTROL" \
      --rmats_sample1 "${g}_$c" --pval 0.05 --output "grase_vs_rmats.$a.txt" > "compare.$a.log" 2>&1) \
      || { echo "   compare $a FAILED (see $o/compare.$a.log)"; continue; }
  done
  # 5. UpSet plots on exonic parts, GrASE (both classes) vs rMATS
  grase_files=$(for a in $CLASSES; do printf '%s,' "$WD/test_anc12/$g/grase_exsj/test_bipartition.exsj.${a}_betabinom_EBmap.annotated.txt"; done)
  map_dirs=$(for a in $CLASSES; do printf '%s,' "$o/map.$a"; done)
  (cd "$o" && Rscript "$S/compare_parts_upset.R" --grase "${grase_files%,}" --map_dirs "${map_dirs%,}" \
    --class "$(echo $CLASSES | tr ' ' ',')" --rmats_dir "$RMATS/$d/out" --contrast "${g}_${c}_vs_${g}_$CONTROL" \
    --pval 0.05 --outdir upset > upset.log 2>&1) || echo "   UpSet FAILED (see $o/upset.log)"
done
echo "== done ($(date))"
