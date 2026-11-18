#!/usr/bin/env bash
# Whole pipeline: setup checks -> 02 (Scallop2) -> 03 (annotation) -> rMATS (08, 09, 10) and GrASE
# (11, 12, 13) in parallel. Settings: pipeline/config.sh. Each step skips work already done, so a
# stopped run can be restarted with the same command.
# Usage (from anywhere): nohup pipeline/run_all.sh > run_all.log 2>&1 &
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
. "${CONFIG:-pipeline/config.sh}"
echo "== run_all ($(date)) in $PWD"
pipeline/setup.sh
mkdir -p assembly_unique
echo "== step 02 ($(date))"; pipeline/02_assemble_scallop2.sh > assembly_unique/step02.log 2>&1 \
  || { tail -5 assembly_unique/step02.log; echo "step 02 FAILED (assembly_unique/step02.log)"; exit 1; }
grep -E 'samples ->|transcripts per sample|FAIL' assembly_unique/step02.log
echo "== step 03 ($(date))"; pipeline/03_union_filter.sh > assembly_unique/step03.log 2>&1 \
  || { tail -5 assembly_unique/step03.log; echo "step 03 FAILED (assembly_unique/step03.log)"; exit 1; }
grep -E '^E |done' assembly_unique/step03.log
echo "== downstream ($(date))"; pipeline/run_downstream.sh
echo "== run_all done ($(date))"
