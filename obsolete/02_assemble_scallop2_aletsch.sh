#!/usr/bin/env bash
# Alternative assemblies from the same 63 STAR BAMs:
#   1. Scallop2 v1.1.2 per sample (bioconda build h503566f_7: the newest build h5ca1c30_8
#      uses AVX2/BMI2 and crashes with "Illegal instruction" on this Xeon E5-2643)
#        -> scallop2/<taxon>/<S>.scallop2.gtf (+ gtf_list.scallop2.txt)
#   2. Aletsch (bioconda) multi-sample meta-assembly of all BAMs -> aletsch/aletsch.gtf
# Library: 101 bp paired-end, fr-firststrand (StringTie --rf / rMATS fr-firststrand).
# Restartable: finished samples / steps are skipped.
#
# Usage (from the project directory):
#   nohup pipeline/02_assemble_scallop2_aletsch.sh > assembly_scallop2_aletsch/run_assembly.log 2>&1 &
set -euo pipefail

PROJ=/srv/jupyter-homes/jupyter-eynizade/waterstudy_metrosideros
OUT=$PROJ/assembly_scallop2_aletsch
BAM_DIR=$PROJ/BAMs/star_alignmentuntrimmed_results
BAM_SUFFIX=.Aligned.sortedByCoord.out.bam
SCALLOP2=/srv/jupyter-homes/jupyter-eynizade/.conda/envs/scallop2/bin/scallop2
ALETSCH=/srv/jupyter-homes/jupyter-eynizade/.conda/envs/aletsch/bin/aletsch
JOBS=${JOBS:-12}
THREADS=${THREADS:-16}

cd "$PROJ"
# sample -> taxon, from the StringTie sample list
sed 's#^\(taxon[A-Z]\)/.*/\(S[0-9]*\)\.stringtie\.gtf$#\1 \2#' gtf_list.txt > "$OUT/samples.txt"

# --- 1. Scallop2 per sample ----------------------------------------------------
run_scallop2() {
  local taxon=$1 s=$2 out=$OUT/scallop2/$1/$2.scallop2.gtf
  [[ -s $out ]] && { echo "skip $s"; return; }
  mkdir -p "$(dirname "$out")"
  "$SCALLOP2" -i "$BAM_DIR/$s$BAM_SUFFIX" -o "$out.tmp" --library_type first --verbose 0 \
    > "$OUT/scallop2/$taxon/$s.log" 2>&1 && mv "$out.tmp" "$out" && echo "ok   $s" || { echo "FAIL $s"; return 1; }
}
export -f run_scallop2; export OUT BAM_DIR BAM_SUFFIX SCALLOP2
echo "== Scallop2 ($(date))"
xargs -a "$OUT/samples.txt" -P "$JOBS" -n 2 bash -c 'run_scallop2 "$0" "$1"' | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
awk -v d="assembly_scallop2_aletsch/scallop2" '{print d"/"$1"/"$2".scallop2.gtf"}' "$OUT/samples.txt" > "$OUT/gtf_list.scallop2.txt"

# --- 2. Aletsch ---------------------------------------------------------------
A=$OUT/aletsch; mkdir -p "$A/profile" "$A/gtf"
awk -v d="$BAM_DIR" -v x="$BAM_SUFFIX" '{print d"/"$2 x" "d"/"$2 x".bai paired_end"}' "$OUT/samples.txt" > "$A/bam_list.txt"
if [[ ! -s $A/aletsch.gtf ]]; then
  nsamp=$(wc -l < "$A/bam_list.txt")
  if [[ $(ls "$A/profile"/*.profile 2>/dev/null | wc -l) -ne $nsamp ]]; then
    echo "== Aletsch profile ($(date))"
    # aletsch prints one debug line per read ("Hit ..."); keep those out of the logs
    set +e
    "$ALETSCH" --profile -i "$A/bam_list.txt" -p "$A/profile" -t "$THREADS" 2>&1 | grep -v '^Hit ' > "$A/profile.log"
    st=${PIPESTATUS[0]}; set -e
    [[ $st -eq 0 ]] || { echo "Aletsch profile FAILED (exit $st)"; exit 1; }
  fi
  # Aletsch infers strandedness only from XS tags; these STAR BAMs have none
  # (outSAMstrandField None), so every sample profiles as unstranded (library_type 0)
  # and no splice graphs get built. Set fr-firststrand (FR_FIRST = 1 in
  # util/constants.h); strand is then assigned per read from the SAM flags.
  mkdir -p "$A/profile.inferred"
  for f in "$A/profile"/*.profile; do
    [[ -s $A/profile.inferred/$(basename "$f") ]] || cp "$f" "$A/profile.inferred/"
    sed -i 's/^library_type .*/library_type 1/' "$f"
  done
  rm -rf "$A/gtf"; mkdir -p "$A/gtf"
  echo "== Aletsch assemble ($(date))"
  set +e
  "$ALETSCH" -i "$A/bam_list.txt" -o "$A/aletsch.gtf.tmp" -p "$A/profile" -d "$A/gtf" \
    -t "$THREADS" -c $((2 * nsamp)) 2>&1 | grep -v '^Hit ' > "$A/aletsch.log"
  st=${PIPESTATUS[0]}; set -e
  [[ $st -eq 0 ]] || { echo "Aletsch assemble FAILED (exit $st)"; exit 1; }
  mv "$A/aletsch.gtf.tmp" "$A/aletsch.gtf"
fi
echo "== done ($(date))"
