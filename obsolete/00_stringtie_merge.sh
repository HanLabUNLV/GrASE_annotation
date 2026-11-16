#!/usr/bin/env bash
# Step 00: StringTie per-sample assemblies (done on the original cluster) and the merge.
#
# Per-sample assembly, StringTie 2.2.3, on the original PBS cluster
# (taxon*/jobs/00_stringtie_build_gtf.pbs; command copied from each GTF header, e.g. S15):
#   stringtie <S>.Aligned.sortedByCoord.out.bam \
#     -G .../commongarden/splicing_analysis_20260308/inputs/refs/annotation.rmats.gtf \
#     --rf -o taxon<X>/stringtie/per_sample_gtf/<S>.stringtie.gtf -p 22
# The guide annotation.rmats.gtf is not on this server, so that step is not rerun here.
#
# Guide GTF for the merge: Metrosideros_polymrpha_incana.scaffolded.gtf, converted from the
# .gff with AGAT v1.7.0 (GTFrelax output), from the user's shell history:
#   conda create -n agat_env -c conda-forge -c bioconda agat
#   agat_convert_sp_gff2gtf.pl --gff Metrosideros_polymrpha_incana.scaffolded.gff -o Metrosideros_polymrpha_incana.scaffolded.gtf
#
# Merge (StringTie 2.2.3; command copied from the merged.f0.01.gtf header, where the output
# was written as merged3.gtf and renamed). Output: merged.f0.01.gtf.
#
# Usage: pipeline/00_stringtie_merge.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
STRINGTIE=${STRINGTIE:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/stringtie22/bin/stringtie}
OUT=${OUT:-merged.f0.01.gtf}

"$STRINGTIE" --merge -G Metrosideros_polymrpha_incana.scaffolded.gtf -c 2 -T 1 -f 0.01 -m 200 \
  -o "$OUT" gtf_list.txt
