#!/usr/bin/env bash
# Downstream analyses on one annotation: rMATS (08 -> 12 -> 11) and GrASE (09 -> 10 -> 13) in parallel.
# Outputs: rmats[_TAG]/, rmats[_TAG]_prefiltered/, grase[_TAG]/ (this project's 2026-11 build used TAG=scallop2).
# Usage: [GTF=assembly_unique/annotation.final.gtf TAG= RMATS_THREADS=4 GRASE_JOBS=12] pipeline/run_downstream.sh
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
GTF=$(readlink -f "${GTF:-assembly_unique/annotation.final.gtf}")
TAG=${TAG:-}
R=$PROJ/rmats${TAG:+_$TAG}; RP=${R}_prefiltered; W=$PROJ/grase${TAG:+_$TAG}
[[ -s $GTF ]] || { echo "missing annotation $GTF (run steps 02 and 03)"; exit 1; }
echo "== downstream on $GTF ($(date)): rMATS -> $R, $RP; GrASE -> $W"
(
  mkdir -p "$R"
  GTF=$GTF BASE=$R NTHREAD=${RMATS_THREADS:-4} pipeline/08_run_rmats.sh > "$R/run_rmats.log" 2>&1 || { echo "rMATS 08 FAILED (see $R/run_rmats.log)"; exit 1; }
  echo "rMATS 08 done ($(date))"
  mkdir -p "$RP"
  SRC=$R BASE=$RP NTHREAD=${RMATS_THREADS:-4} pipeline/12_rmats_prefiltered.sh > "$RP/run.log" 2>&1 || { echo "rMATS 12 FAILED (see $RP/run.log)"; exit 1; }
  echo "rMATS 12 done ($(date))"
  "$PYTHON" -I pipeline/11_filter_rmats.py "$RP" "$GTF" > "$RP/filter.log" 2>&1 || { echo "rMATS 11 FAILED (see $RP/filter.log)"; exit 1; }
  echo "rMATS 11 done ($(date))"; cat "$RP/filter.log"
) &
pr=$!
(
  mkdir -p "$W"
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 GTF=$GTF WD=$W JOBS=${GRASE_JOBS:-$JOBS} \
    pipeline/09_grase_annotation.sh > "$W/09_annotation.log" 2>&1 || { echo "GrASE 09 FAILED (see $W/09_annotation.log)"; exit 1; }
  echo "GrASE 09 done ($(date))"; grep -aE 'genes:|bipartition files|internal:|events' "$W/09_annotation.log" | tail -4
  WD=$W JOBS=${GRASE_JOBS:-$JOBS} pipeline/10_grase_count.sh > "$W/10_count.log" 2>&1 || { echo "GrASE 10 FAILED (see $W/10_count.log)"; exit 1; }
  echo "GrASE 10 done ($(date))"
  WD=$W JOBS=${GRASE_JOBS:-$JOBS} pipeline/13_grase_test.sh > "$W/13_grase_test.log" 2>&1 || { echo "GrASE 13 FAILED (see $W/13_grase_test.log)"; exit 1; }
  echo "GrASE 13 done ($(date))"; grep '   3 ' "$W/13_grase_test.log"
) &
pg=$!
st=0; wait $pr || st=1; wait $pg || st=1
echo "== all done ($(date)), status $st"
exit $st
