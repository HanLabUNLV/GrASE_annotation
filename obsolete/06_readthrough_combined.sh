#!/usr/bin/env bash
# Step 06: read-through filter on the combined annotation (steps 7-9 of step 01):
# drop transcripts with a minor gene-bridging junction (bridge reads < 0.2 x the bridged
# genes' internal junction reads) unless elevated under stress, then one gene_id per
# contiguous gene piece.
# Outputs: merged.junction.final.gtf, support/junction_qc/junction.bridging_ratio.tsv,
#          support/junction_qc/junction.readthrough_by_condition.tsv
#
# Usage: pipeline/06_readthrough_combined.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory

ONLY_READTHROUGH=1 FINAL_GTF=merged.junction.gtf NORT_GTF=merged.junction.final.gtf \
RT_PREFIX=support/junction_qc/junction \
pipeline/01_filter_support.sh
