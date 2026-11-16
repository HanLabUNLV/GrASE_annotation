#!/usr/bin/env bash
# Step 03: junction-level filter of the raw Aletsch meta-assembly (step 01 on aletsch.gtf):
# every intron canonical, <= 20 kb, >= 3 unique reads in >= 3 samples, each with >= 1 read
# anchored >= 12 bp on both sides (overhang rule, since 2026-11-12).
# Also writes the per-transcript sample-support table (exact intron-chain matches of the
# Scallop2 per-sample assemblies, gffcompare "="), used by the evidence rule in step 05.
# Outputs: assembly_scallop2_aletsch/aletsch.canonical.sj3x3.anc12.max20kb.gtf,
#          assembly_scallop2_aletsch/aletsch.anc12.intron_qc.tsv,
#          assembly_scallop2_aletsch/aletsch.transcript_support.tsv
#          (before 2026-11-12, no overhang rule: aletsch.canonical.sj3x3.max20kb.gtf, aletsch.intron_qc.tsv)
#
# Usage: pipeline/03_filter_aletsch.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
A=assembly_scallop2_aletsch

REF=$A/aletsch/aletsch.gtf LIST=$A/gtf_list.scallop2.txt GTF_SUFFIX=.scallop2.gtf TAXON_FIELD=3 JOBS=4 \
OUTDIR=$A/support_aletsch TABLE=$A/aletsch.transcript_support.tsv \
FINAL_GTF=$A/aletsch.canonical.sj3x3.anc12.max20kb.gtf INTRON_TABLE=$A/aletsch.anc12.intron_qc.tsv \
pipeline/01_filter_support.sh
