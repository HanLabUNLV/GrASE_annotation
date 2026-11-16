#!/usr/bin/env bash
# Step 07: side-by-side QC of every annotation version (pipeline/lib/annotation_qc.py):
# transcripts, loci, isoforms per locus, intron lengths, non-canonical and weakly supported
# introns, gene-bridging introns, transcripts spanning >= 2 MAKER genes.
# Output: support/junction_qc/annotation_qc.tsv
#
# Usage: pipeline/07_annotation_qc.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PYTHON=${PYTHON:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/grase/bin/python}
OUT=${OUT:-support/junction_qc/annotation_qc.tsv}
A=assembly_scallop2_aletsch

"$PYTHON" -I pipeline/lib/annotation_qc.py \
  Metrosideros_polymrpha_incana.scaffolded.fa Metrosideros_polymrpha_incana.scaffolded.gtf support/junction_counts \
  "MAKER=Metrosideros_polymrpha_incana.scaffolded.gtf" \
  "ST.f0.01=merged.f0.01.gtf" \
  "ST.junction=merged.f0.01.canonical.sj3x3.max20kb.gtf" \
  "AL.raw=$A/aletsch/aletsch.gtf" \
  "AL.junction=$A/aletsch.canonical.sj3x3.max20kb.gtf" \
  "junction=merged.junction.gtf" \
  "junction.final=merged.junction.final.gtf" \
  "archived.combined.final=archive/intron_chain_filter/merged.combined.final.gtf" \
  > "$OUT"
column -t -s$'\t' "$OUT"
