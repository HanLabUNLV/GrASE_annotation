#!/usr/bin/env bash
# Step 09: rMATS statistics on the events that pass the junction rules only.
# Step 08 (--novelSS --statoff) reads the BAMs once and detects and counts all events,
# including events at novel splice sites that never passed the annotation's junction
# rules. Here, per contrast:
#   1. keep only events whose every junction is canonical, <= 20 kb and has >= 3 unique
#      reads in >= 3 samples (pipeline/lib/prefilter_rmats_events.py; same rules as step 01),
#      and (since 2026-11-12; ANCHOR_DIR="" turns it off) those samples each have >= 1 read
#      with >= 12 bp aligned on both sides (overhang rule; SE/MXE inner exons < 12 bp exempt).
#   2. copy step 08's event and count tables restricted to those events
#      (fromGTF.*, JC/JCEC.raw.input.*; pipeline/lib/subset_rmats_counts.py)
#   3. rmats.py --task stat on them, so FDR is computed over the kept events only.
# Nothing is recounted: counts, lengths and p-values of kept events are those of step 08.
# (rmats.py --task post --fixed-event-set was tried first: it swaps IncFormLen/SkipFormLen
#  for ~40% of MXE events when reloading them, changing their p-values and dPSI.)
# Output: $BASE/<group>/<contrast>/{events,out}; then filter with step 10.
#
# Usage: [SRC=rmats BASE=rmats_prefiltered ANCHOR_DIR=support/junction_counts_anchor] \
#          pipeline/09_rmats_prefiltered.sh
# (contrasts with $BASE/<group>/<contrast>/out/summary.txt are skipped: use a new BASE to redo them)
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
. "${CONFIG:-pipeline/config.sh}"
SRC=$(readlink -f "${SRC:-rmats}")
BASE=$(readlink -m "${BASE:-rmats_prefiltered}")
ENV=$RMATS_ENV
NTHREAD=${NTHREAD:-$THREADS}
ANCHOR_DIR=${ANCHOR_DIR-support/junction_counts_anchor}   # "" = no overhang rule
RMATS=("$ENV/bin/python" "$ENV/rMATS/rmats.py")
export PATH=$ENV/bin:$PATH

echo "rMATS $("${RMATS[@]}" --version)  stats on prefiltered events from $SRC  start $(date)"
for d in "$SRC"/*/*_vs_"$CONTROL"; do
  [[ -d $d/out ]] || continue
  taxon=$(basename "$(dirname "$d")"); contrast=$(basename "$d")
  o=$BASE/$taxon/$contrast
  if [[ -s $o/out/summary.txt ]]; then echo "skip $taxon $contrast (done)"; continue; fi
  rm -rf "$o"; mkdir -p "$o/events"
  cp "$d/b1.txt" "$d/b2.txt" "$o/"
  echo "== $taxon $contrast  $(date)"
  "$PYTHON" -I pipeline/lib/prefilter_rmats_events.py "$d/out" "$o/events" support/junction_counts \
    "$GENOME" 20000 3 3 "$ANCHOR_DIR" > "$o/events_kept.tsv" 2> "$o/prefilter.log" || { cat "$o/prefilter.log"; exit 1; }
  awk '{printf "   %s %d -> %d\n", $1, $2, $3}' "$o/events_kept.tsv"
  "$PYTHON" -I pipeline/lib/subset_rmats_counts.py "$d/out" "$o/events" "$o/out" > "$o/subset.log"
  "${RMATS[@]}" --task stat --od "$o/out" --tmp "$o/tmp" --nthread "$NTHREAD" > "$o/rmats.log" 2>&1 \
    && echo "   done $taxon $contrast  $(date)" \
    || { echo "   FAILED $taxon $contrast (see $o/rmats.log)"; exit 1; }
  rm -rf "$o/tmp"
done
echo "all contrasts finished $(date)"
