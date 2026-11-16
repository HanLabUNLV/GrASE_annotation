#!/usr/bin/env bash
# Step 13: GrASE Stages 2, 2b and 3 (GrASE/README.md "Quick Start"): differential path usage
# per group ($SAMPLES), each condition vs $CONTROL, for internal and terminal events.
# Follows the repo's worked example GrASE/scripts/GSE188847/ (stages 30, 40, 50, 60), with all
# conditions per group (--conditions) and all contrasts in one run instead of one.
#   2.  exoncnt.R: sum DEXSeq exonic-part counts (step 10) onto each bipartition's shared
#       reference and two distinct sets -> <group>/bipartition.ex.<class>.counts/
#   2b. bipartition_sjcnt.R + merge_exon_sj_counts.R: where a side has no exclusive exonic part,
#       use its exclusive junctions' split reads instead -> <group>/bipartition.exsj.<class>.counts/
#       Junction reads come from $SJ_DIR/<group>_<condition>/<S>.SJ.out.tab, written by this
#       script from the step 01 BAM cache (support/junction_counts: unique, MAPQ 255), in
#       SJ.out.tab layout (unique reads in column 7; other columns 0). STAR's own SJ.out.tab is
#       not used, for the reasons in pipeline/README.md ("Why junction support is counted
#       from the BAMs").
#   3.  exontest.R with the GSE188847 stage 60 settings: model betabinom_EBmap, --dpi_gate=both,
#       --padj_threshold=0.05, --independent_filtering, per-base pi (--gff_dir, --read_length 101),
#       phi file phi.merged.<class>.txt; one run per group and class with all its contrasts, so phi
#       is estimated over all the group's conditions -> <group>/grase_exsj/ (both classes)
#       Other thresholds are exontest.R defaults (min_dpi 0.05, min_reads 10 in one group,
#       min_effect |log OR| 0.5); adjustment PADJ_METHOD=BH (plain Benjamini-Hochberg per
#       contrast; nested_BH = GrASE's gene screen, its default before 2026-11-14); MIN_DPI=<x> overrides min_dpi.
# Events: $BIP, by default bipartition.filtered.anc12/ (step 09 D, overhang filter).
# Restartable: a stage whose output exists is skipped (delete the directory to redo it).
#
# Groups and conditions come from $SAMPLES (pipeline/config.sh): per group, every condition is tested
# against $CONTROL in one run (so the dispersion is estimated over all of the group's conditions).
# GrASE's R scripts write error logs to their working directory: they are run from $OUT/<group>/.
# Usage: [WD=grase RUN_GROUPS="<g> ..." CLASSES="internal terminal" JOBS=14] pipeline/13_grase_test.sh
# Results: $OUT/<group>/grase_exsj/test_bipartition.exsj.<class>_betabinom_EBmap.annotated.txt
#   (column contrast: <group>_<condition>_vs_<group>_<CONTROL>)
# A rerun reuses existing fits and only redoes the annotation (exontest.R); to re-test, delete
#   $OUT/<group>/grase_exsj/test_bipartition.exsj.<class>_* and phi.merged.<class>* first.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
PROJ=$PWD
. "${CONFIG:-pipeline/config.sh}"
WD=$(readlink -m "${WD:-grase}")
BIP=$(readlink -m "${BIP:-$WD/bipartition.filtered.anc12}")
OUT=$(readlink -m "${OUT:-$WD/test_anc12}")
SJ_DIR=$(readlink -m "${SJ_DIR:-$WD/SJ}")
RUN_GROUPS=${RUN_GROUPS:-$(pipeline_groups | paste -sd' ')}
CLASSES=${CLASSES:-internal terminal}
MODEL=${MODEL:-betabinom_EBmap}
PADJ_METHOD=${PADJ_METHOD:-BH}   # plain BH over all events of a contrast (since 2026-11-14); nested_BH = gene screen
MIN_DPI=${MIN_DPI:-}          # empty = exontest.R default (0.05)
READ_LENGTH=${READ_LENGTH:-$READ_LEN}
export PATH=$GRASE_ENV/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
S=$(readlink -f "$GRASE_DIR")/Rpkg/scripts
GFF=$(ls "$WD"/ref/*.dexseq.bygene.gff)

echo "== GrASE stages 2-3: $BIP -> $OUT ($(date))"

# junction files for 2b (from the step 01 cache; rebuilt if missing)
while read -r s t c; do
  f=$SJ_DIR/${t}_$c/$s.SJ.out.tab
  [[ -s $f ]] && continue
  mkdir -p "$(dirname "$f")"
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,0,0,0,$4,0,0}' "$PROJ/support/junction_counts/$s.junc.tsv" \
    | sort -k1,1 -k2,2n -k3,3n > "$f.tmp" && mv "$f.tmp" "$f"
done < "$CONDITIONS"

for t in $RUN_GROUPS; do
  conds=$(pipeline_conditions "$t" | sed "s/^/${t}_/" | paste -sd,)
  ctr=$(pipeline_treatments "$t" | sed "s/.*/${t}_&:${t}_$CONTROL/" | paste -sd,)
  for a in $CLASSES; do
    o=$OUT/$t; ex=$o/bipartition.ex.$a.counts; sj=$o/sjcnt.$a; exsj=$o/bipartition.exsj.$a.counts
    test=$o/grase_exsj; mkdir -p "$o" "$test"
    echo "-- $t $a ($(date))"
    # 2. exonic-part counts onto bipartitions
    if [[ ! -s $ex/bipartition.ex.$a.cnt.combined.txt ]]; then
      rm -rf "$ex"
      (cd "$o" && Rscript "$S/exoncnt.R" -c "$WD/DEXSeq/count_files" -t bipartition --conditions="$conds" -a "$a" \
        -i "$BIP" -o "$ex" --mc_cores "$JOBS" > "$o/exoncnt.$a.log" 2>&1)
    fi
    echo "   2  exon counts: $(ls "$ex" | grep -c '\.cnt\.txt$') genes"
    # 2b. junction counts where a distinct set is empty (start each output dir empty)
    if [[ ! -s $exsj/bipartition.exsj.$a.cnt.combined.txt ]]; then
      rm -rf "$sj" "$exsj"
      (cd "$o" && Rscript "$S/bipartition_sjcnt.R" --inputdir "$BIP" --graphmldir "$WD/graphml" --gff "$GFF" \
        --sjdir "$SJ_DIR" --conditions "$conds" --output "$sj" --type "$a" --cores "$JOBS" > "$o/sjcnt.$a.log" 2>&1)
      Rscript "$S/merge_exon_sj_counts.R" --exon_counts "$ex" --sj_counts "$sj" --output "$exsj" > "$o/merge.$a.log" 2>&1
      for f in "$ex"/*.bipartition.txt; do ln -sf "$(realpath "$f")" "$exsj/"; done
      python3 "$S/combine_exoncnt.py" "$exsj" "$exsj/bipartition.exsj.$a.cnt.combined.txt.tmp" > "$o/combine.$a.log" 2>&1
      mv "$exsj/bipartition.exsj.$a.cnt.combined.txt.tmp" "$exsj/bipartition.exsj.$a.cnt.combined.txt"
    fi
    [[ -s $sj/bipartition_sjcnt.errors.log ]] && echo "   2b WARNING: $(wc -l < "$sj/bipartition_sjcnt.errors.log") gene errors, see $sj/bipartition_sjcnt.errors.log"
    echo "   2b $(grep -h 'from SJ' "$o/merge.$a.log" | tr '\n' ' ')"
    # 3. test
    (cd "$o" && Rscript "$S/exontest.R" --file="bipartition.exsj.$a.cnt.combined.txt" --outdir="$test" --countdir="$exsj/" \
      --splittype=bipartition --gff_dir="$WD/dexseq.gff" --read_length="$READ_LENGTH" --dpi_gate=both \
      --padj_threshold=0.05 --padj_method="$PADJ_METHOD" --independent_filtering ${MIN_DPI:+--min_dpi="$MIN_DPI"} --phi="phi.merged.$a.txt" \
      --model="$MODEL" --mc_cores="$JOBS" --contrasts="$ctr" \
      > "$test/test_$a.log" 2>&1)
    r=$test/test_bipartition.exsj.${a}_$MODEL.annotated.txt
    [[ -s $r ]] || { echo "   3  FAILED: no $r (see $test/test_$a.log)"; exit 1; }
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
      {n[$h["contrast"]]++; if($h["significant"]=="TRUE") s[$h["contrast"]]++}
      END{for(c in n) printf "   3  %s: %d tests, %d significant\n", c, n[c], s[c]}' "$r"
  done
done
echo "== done ($(date))"
