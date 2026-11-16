#!/usr/bin/env bash
# Step 05: combined annotation = StringTie (step 01) + Aletsch (step 03) transcripts that
# add a junction missing from StringTie (intron-retention variants, class "m", excluded).
# Added transcripts take the overlapping StringTie gene_id, or a new ALG.<n> locus id.
# Then the evidence rule (EVIDENCE=1, default): transcripts whose exact intron chain no
# sample assembled are removed when all their junctions remain in other transcripts
# (pipeline/lib/evidence_filter.py; support tables from steps 01 and 03).
# Details: pipeline/lib/build_combined_gtf.py. Output: merged.junction.gtf
#
# Usage: [ST=<stringtie.gtf> AL=<aletsch.gtf> TMAP=<step 04 .tmap> OUT=<out.gtf>] pipeline/05_build_combined.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PYTHON=${PYTHON:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/grase/bin/python}
OUT=${OUT:-merged.junction.gtf}

ST=${ST:-merged.f0.01.canonical.sj3x3.anc12.max20kb.gtf}
AL=${AL:-assembly_scallop2_aletsch/aletsch.canonical.sj3x3.anc12.max20kb.gtf}
TMAP=${TMAP:-assembly_scallop2_aletsch/compare_vs_stringtie/AL_vs_ST.AL.final.gtf.tmap}   # from step 04

ST_SUPPORT=${ST_SUPPORT:-support/merged.f0.01.transcript_support.tsv}              # from step 01
AL_SUPPORT=${AL_SUPPORT:-assembly_scallop2_aletsch/aletsch.transcript_support.tsv}   # from step 03
EVIDENCE=${EVIDENCE:-1}   # 1 = apply the evidence rule (pipeline/lib/evidence_filter.py)

if [[ $EVIDENCE == 1 ]]; then
  tmp=$(mktemp -p "$(dirname "$OUT")" .combined.XXXX.gtf); trap 'rm -f "$tmp"' EXIT
  "$PYTHON" -I pipeline/lib/build_combined_gtf.py "$ST" "$AL" "$TMAP" "$tmp"
  "$PYTHON" -I pipeline/lib/evidence_filter.py "$tmp" "$OUT" "$ST_SUPPORT" "$AL_SUPPORT"
else
  "$PYTHON" -I pipeline/lib/build_combined_gtf.py "$ST" "$AL" "$TMAP" "$OUT"
fi
