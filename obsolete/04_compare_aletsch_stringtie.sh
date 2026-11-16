#!/usr/bin/env bash
# Step 04: compare the filtered Aletsch (step 03) and StringTie (step 01) annotations
# with gffcompare. StringTie is the reference; the .tmap class codes are
# used by step 05 to assign added Aletsch transcripts to StringTie genes.
# Output: assembly_scallop2_aletsch/compare_vs_stringtie/AL_vs_ST.* (.stats, .tmap, .refmap)
#
# Usage: [ST=<stringtie.gtf> AL=<aletsch.gtf> D=<outdir>] pipeline/04_compare_aletsch_stringtie.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
GFFCOMPARE=${GFFCOMPARE:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/stringtie22/bin/gffcompare}
ST=$(readlink -f "${ST:-merged.f0.01.canonical.sj3x3.anc12.max20kb.gtf}")
AL=$(readlink -f "${AL:-assembly_scallop2_aletsch/aletsch.canonical.sj3x3.anc12.max20kb.gtf}")
D=${D:-assembly_scallop2_aletsch/compare_vs_stringtie}
mkdir -p "$D"
ln -sf "$AL" "$D/AL.final.gtf"
ln -sf "$ST" "$D/ST.final.gtf"
cd "$D"
"$GFFCOMPARE" -r ST.final.gtf -o AL_vs_ST AL.final.gtf
