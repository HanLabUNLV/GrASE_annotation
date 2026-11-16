#!/usr/bin/env bash
# Step 02b (comparison, 2026-11-13): StringTie 3.0.0-only assembly with strict settings, on uniquely
# mapped reads, as an alternative to step 02 (Scallop2 + merge).
#   1. per sample: samtools view -h -q MIN_MAPQ -F 0x904 <bam> | stringtie - --rf -g 0 -c 5 -f 0.15 -j 3
#        -q 255 -F 0x904   unique, primary, not QC-fail/supplementary (StringTie has no MAPQ
#                          option and counts multi-mapped reads as 1/NH; streamed, so no
#                          filtered BAM copies; output identical to a filtered BAM, tested)
#        --rf              fr-firststrand (dUTP)
#        -g 0              no bridging of gaps between read mappings (default 50 bp)
#        -c 5              min reads per bp for multi-exon transcripts (default 1)
#        -f 0.15           min isoform fraction (default 0.01)
#        -j 3              min reads per junction (default 1)
#      no -G: the guide annotation would force existing gene-model fusions back in
#      -> stringtie/<taxon>/<S>.stringtie.gtf
#   2. stringtie --merge -g 0 -c 3 -T 2 -f 0.1 (no -G) -> merge/stringtie_merge.gtf
#      (StringTie GTFs carry TPM, so -T applies here, unlike on Scallop2 output)
# Restartable: finished samples are skipped.
# Usage: [OUT=assembly_stringtie_strict JOBS=4] nohup pipeline/02b_assemble_stringtie_strict.sh > <OUT>/run_assembly.log 2>&1 &
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
OUT=$(readlink -m "${OUT:-assembly_stringtie_strict}")
MIN_MAPQ=${MIN_MAPQ:-255}
CONDITIONS=${CONDITIONS:-$PROJ/support/junction_qc/sample_conditions.tsv}
BAM_DIR=${BAM_DIR:-$PROJ/BAMs/star_alignmentuntrimmed_results}
BAM_SUFFIX=${BAM_SUFFIX:-.Aligned.sortedByCoord.out.bam}
ENVS=/srv/jupyter-homes/jupyter-eynizade/.conda/envs
SAMTOOLS=$ENVS/rmats420/bin/samtools
STRINGTIE=${STRINGTIE:-$ENVS/stringtie3/bin/stringtie}   # 3.0.0, bioconda build h29c0135_0 (3.0.1-3.0.3 builds need AVX2; crash on this CPU)
JOBS=${JOBS:-4}
ST_OPTS=${ST_OPTS:---rf -g 0 -c 5 -f 0.15 -j 3}
MERGE_OPTS=${MERGE_OPTS:--g 0 -c 3 -T 2 -f 0.1}

mkdir -p "$OUT/merge"
awk -F'\t' '{print $2, $1}' "$CONDITIONS" | sort -k1,1 -k2,2V > "$OUT/samples.txt"
echo "== $(wc -l < "$OUT/samples.txt") samples -> $OUT ($(date)); MAPQ>=$MIN_MAPQ; stringtie $ST_OPTS; merge $MERGE_OPTS"

run_st() {
  local taxon=$1 s=$2 out=$OUT/stringtie/$1/$2.stringtie.gtf
  [[ -s $out ]] && { echo "skip $s"; return; }
  mkdir -p "$(dirname "$out")"
  set -o pipefail
  "$SAMTOOLS" view -h -q "$MIN_MAPQ" -F 0x904 "$BAM_DIR/$s$BAM_SUFFIX" 2> "$out.samtools.log" \
    | "$STRINGTIE" - $ST_OPTS -l "$s" -o "$out.tmp" > "$OUT/stringtie/$taxon/$s.log" 2>&1 \
    && mv "$out.tmp" "$out" && rm -f "$out.samtools.log" && echo "ok   $s" || { echo "FAIL $s"; return 1; }
}
export -f run_st; export OUT BAM_DIR BAM_SUFFIX SAMTOOLS STRINGTIE MIN_MAPQ ST_OPTS
echo "== StringTie per sample ($(date))"
xargs -a "$OUT/samples.txt" -P "$JOBS" -n 2 bash -c 'run_st "$0" "$1"' \
  | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
awk -v d="$OUT/stringtie" '{print d"/"$1"/"$2".stringtie.gtf"}' "$OUT/samples.txt" > "$OUT/gtf_list.stringtie.txt"
while read -r f; do [[ -s $f ]] || { echo "missing $f"; exit 1; }; done < "$OUT/gtf_list.stringtie.txt"

M=$OUT/merge/stringtie_merge.gtf
if [[ ! -s $M ]]; then
  echo "== stringtie --merge ($(date))"
  "$STRINGTIE" --merge $MERGE_OPTS -o "$M.tmp" "$OUT/gtf_list.stringtie.txt" && mv "$M.tmp" "$M"
fi
echo "   $M: $(awk '$3=="transcript"' "$M" | wc -l) transcripts, $(awk -F'\t' '$3=="transcript"{match($9,/gene_id "[^"]+"/); print substr($9,RSTART,RLENGTH)}' "$M" | sort -u | wc -l) genes"
echo "== done ($(date))"
