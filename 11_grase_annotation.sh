#!/usr/bin/env bash
# Step 11: GrASE annotation files from one GTF (GrASE README, "Building the annotation
# files from scratch", sections A-C):
#   A. per-gene GTFs, per-gene DEXSeq GFFs (+ concatenated GFF for counting), splicing
#      graphs (GrASE/scripts/createfiles.sh)
#   B. bubble enumeration -> bipartitions (Rpkg/scripts/bubble_path_split.R)
#   C. split internal vs terminal events  (Rpkg/scripts/filterTerminal.R)
#   D. overhang filter (since 2026-11-12; ANCHOR_DIR="" skips it): drop events with a path
#      through a junction that fails the step 01 overhang rule (lib/grase_anchor_filter.py).
#      A no-op for annotations built with that rule; a post-hoc filter for older ones.
#      Use bipartition.filtered.anc12/ downstream.
# Output under $WD: ref/, gtf/, dexseq.gff/, graphml/, bipartition/, bipartition.filtered/,
#                   bipartition.filtered.anc12/ (+ .removed.tsv)
#
# Usage: GTF=<annotation.gtf> WD=<grase workdir> pipeline/11_grase_annotation.sh
#   default: assembly_unique/annotation.final.gtf -> grase/
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
GTF=$(readlink -f "${GTF:-assembly_unique/annotation.final.gtf}")
WD=$(readlink -m "${WD:-grase}")
ENV=$GRASE_ENV   # R 4.4 + grase, Python + HTSeq, GNU parallel
export PATH=$ENV/bin:$PATH
GRASE=$(readlink -f "$GRASE_DIR")
ANCHOR_DIR=${ANCHOR_DIR-$PROJ/support/junction_counts_anchor}   # per-sample overhang counts (step 01)

echo "== GrASE annotation: $GTF -> $WD ($(date))"
mkdir -p "$WD"

# A. per-gene files and splicing graphs (skipped if graphs already exist)
if [[ -z $(ls -A "$WD/graphml" 2>/dev/null) ]]; then
  rm -rf "$WD/gtf"   # createfiles.sh appends to per-gene GTFs, so start clean
  PY=$ENV/bin/python bash "$GRASE/scripts/createfiles.sh" "$WD" "$GTF" "$JOBS"
fi
echo "   genes: $(wc -l < "$WD/ref/genelist"), graphs: $(ls "$WD/graphml" | wc -l)  ($(date))"

# B. bipartitions (restartable: genes with output or a .running file are skipped)
rm -f "$WD"/bipartition/*.running
Rscript "$GRASE/Rpkg/scripts/bubble_path_split.R" \
  --graphdir="$WD/graphml" --outdir="$WD/bipartition" --split=bipartition --mc_cores="$JOBS"
echo "   bipartition files: $(ls "$WD/bipartition" | grep -c '\.bipartition\.txt$')  ($(date))"

# C. internal vs terminal
Rscript "$GRASE/Rpkg/scripts/filterTerminal.R" \
  --split_dir="$WD/bipartition" --outdir="$WD/bipartition.filtered" --split_type=bipartition
echo "   internal: $(ls "$WD/bipartition.filtered" | grep -c 'internal'), terminal: $(ls "$WD/bipartition.filtered" | grep -c 'terminal')  ($(date))"

# D. overhang filter
if [[ -n $ANCHOR_DIR ]]; then
  python -I "$PROJ/pipeline/lib/grase_anchor_filter.py" "$WD" "$GTF" "$ANCHOR_DIR" 3 3 12
fi
echo "== done ($(date))"
