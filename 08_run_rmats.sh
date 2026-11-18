#!/usr/bin/env bash
# Step 08: rMATS per group and contrast (each condition vs $CONTROL within a group; groups and
# conditions from $SAMPLES) with the step 03 annotation and --novelSS: reads the BAMs, detects and
# counts events. With the default --statoff, statistics are run in step 09 on the events passing the
# junction rules. b1 = treatment, b2 = control. Library settings from pipeline/config.sh.
# Contrasts with a finished summary.txt are skipped, so the script can be restarted.
# (This project's groups equal the original taxon*/metadata/<contrast>/b{1,2}.txt lists.)
#
# Usage: [GTF=<gtf> BASE=<outdir> NTHREAD=<n> RMATS_EXTRA=<opts>] nohup pipeline/08_run_rmats.sh > <outdir>/run_rmats.log 2>&1 &
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
GTF=$(readlink -f "${GTF:-assembly_unique/annotation.final.gtf}")
BASE=$(readlink -m "${BASE:-rmats}")       # output root; use a new one per GTF
NTHREAD=${NTHREAD:-$THREADS}
RMATS_EXTRA=${RMATS_EXTRA---novelSS --statoff}   # default: detect unannotated splice sites, count only

# rmats.py's shebang resolves to the system python; call the env's python directly
RMATS=("$RMATS_ENV/bin/python" "$RMATS_ENV/rMATS/rmats.py")
export PATH=$RMATS_ENV/bin:$PATH
bam() { readlink -f "$BAM_DIR/$1$BAM_SUFFIX"; }

echo "rMATS $("${RMATS[@]}" --version)  GTF=$GTF  $RMATS_READTYPE $RMATS_LIBTYPE read length $READ_LEN  extra: ${RMATS_EXTRA:-none}  start $(date)"
for g in $(pipeline_groups); do
  for c in $(pipeline_treatments "$g"); do
    contrast=${c}_vs_$CONTROL; d=$BASE/$g/$contrast
    if [[ -s $d/out/summary.txt ]]; then echo "skip $g $contrast (done)"; continue; fi
    mkdir -p "$d"
    pipeline_samples "$g" "$c" | while read -r s; do bam "$s"; done | paste -sd, > "$d/b1.txt"
    pipeline_samples "$g" "$CONTROL" | while read -r s; do bam "$s"; done | paste -sd, > "$d/b2.txt"
    for b in b1 b2; do
      [[ -s $d/$b.txt ]] || { echo "no samples for $b of $g $contrast"; exit 1; }
      tr ',' '\n' < "$d/$b.txt" | while read -r f; do [[ -s $f ]] || { echo "missing BAM $f"; exit 1; }; done
    done
    rm -rf "$d/out" "$d/tmp"; mkdir -p "$d/out" "$d/tmp"
    echo "== $g $contrast  $(date)"
    "${RMATS[@]}" --b1 "$d/b1.txt" --b2 "$d/b2.txt" --gtf "$GTF" --od "$d/out" --tmp "$d/tmp" \
      --readLength "$READ_LEN" --variable-read-length --nthread "$NTHREAD" -t "$RMATS_READTYPE" \
      --libType "$RMATS_LIBTYPE" $RMATS_EXTRA > "$d/rmats.log" 2>&1 \
      && echo "   done $g $contrast  $(date)" \
      || { echo "   FAILED $g $contrast (see $d/rmats.log)"; exit 1; }
  done
done
echo "all contrasts finished $(date)"
