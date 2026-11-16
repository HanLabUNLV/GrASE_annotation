#!/usr/bin/env bash
# Step 01: filter a merged/assembled GTF by junction evidence.
#
# Default (junction level): keep a transcript if every intron is canonical
# (GT-AG/GC-AG/AT-AC), <= MAX_INTRON bp, and has >= INTRON_MIN_READS unique reads in
# >= INTRON_MIN_SAMPLES samples (counted from the BAMs), and (overhang rule, since
# 2026-11-12; ANCHOR_MIN=0 turns it off) each of those samples has >= 1 read with
# >= ANCHOR_MIN (12) bp aligned on both sides of the intron. Introns next to an internal
# exon shorter than ANCHOR_MIN are exempt from the overhang rule (no read can anchor there).
#
# Optional transcript-level filter (MIN > 0; off by default since 2026-11-11): keep a
# transcript only if its exact intron chain (gffcompare class "=") is assembled in >= MIN
# samples within one taxon (taxon = path component TAXON_FIELD of each $LIST entry).
# It removed real, well-supported junctions (212 of 349 original significant SE events
# not testable in rMATS lacked such a junction); outputs made with it are in
# archive/intron_chain_filter/.
#
# Steps:
#   1-2 run when MIN > 0 or SUPPORT_TABLE=1 (default; the table feeds the evidence rule
#   in step 05); 3 runs only when MIN > 0 (otherwise $OUT_GTF = $REF):
#   1. gffcompare each per-sample GTF against the merged GTF -> support/<sample>.*
#      (.tmap/.refmap land next to the input GTF; samples with a complete
#      .stats and .tmap are skipped unless FORCE=1)
#   2. per-transcript support table   -> $TABLE
#   3. filtered GTF                   -> $OUT_GTF
#   4-6. intron filters on $OUT_GTF (drop transcripts with any failing intron)
#      - non-canonical splice motif (not GT-AG/GC-AG/AT-AC; genome = $GENOME)
#      - intron longer than MAX_INTRON (default 20 kb; MAKER has 3 of 129,883 introns
#        > 20 kb; longer junctions were weakly supported and chained genes together)
#      - < INTRON_MIN_READS unique spliced reads in >= INTRON_MIN_SAMPLES samples,
#        counted from the BAMs (per-sample counts cached in $JUNC_DIR)
#      - overhang (ANCHOR_MIN > 0): of those samples, < INTRON_MIN_SAMPLES have a read
#        with >= ANCHOR_MIN aligned bases on both sides (STAR's SJ.out.tab overhang
#        rule; per-sample counts cached in $ANCHOR_DIR by lib/junction_anchor_counts.sh).
#        Catches untrimmed read ends forced into a splice (e.g. 90M65N3M7S). Annotations
#        built before 2026-11-12 (merged.junction.final.gtf) did not apply it; their
#        rMATS/GrASE results are filtered after the fact (steps 11 and 09 D).
#      intron table -> $INTRON_TABLE, final GTF -> $FINAL_GTF
#      (default: merged.f0.01.canonical.sj3x3.anc12.max20kb.gtf;
#       before 2026-11-12, without the overhang rule: merged.f0.01.canonical.sj3x3.max20kb.gtf)
#   7-9. read-through filter (READTHROUGH_FILTER=1 or ONLY_READTHROUGH=1; see step 7)
#      and one gene_id per contiguous gene piece -> $NORT_GTF (default <FINAL_GTF>.nort.gtf)
#
# Usage (paths are relative to the project directory; the script cd's there):
#   pipeline/01_filter_support.sh            # junction-level filter of merged.f0.01.gtf (pipeline default)
#   REF=<gtf> FINAL_GTF=<out.gtf> INTRON_TABLE=<tsv> pipeline/01_filter_support.sh   # another assembly (step 03)
#   ONLY_READTHROUGH=1 FINAL_GTF=<in.gtf> NORT_GTF=<out.gtf> pipeline/01_filter_support.sh   # read-through only (step 06)
#   MIN=3 MAX_INTRON=0 READTHROUGH_FILTER=1 pipeline/01_filter_support.sh
#                                           # the archived transcript-level version (merged.support3pertaxon*)
#
# ---------------------------------------------------------------------------
# Provenance of the current files (copied verbatim from file headers / .stats)
# ---------------------------------------------------------------------------
# Per-sample assemblies (StringTie 2.2.3, original cluster, header of each
# taxon*/stringtie/per_sample_gtf/<S>.stringtie.gtf; shown for S15):
#   /storage/stacylab/WaterStudyRNAseq/splicing_waterstudy/conda_envs/stringtie/bin/stringtie \
#     /storage/stacylab/WaterStudyRNAseq/star_alignmentuntrimmed_results/S15.Aligned.sortedByCoord.out.bam \
#     -G /storage/stacylab/WaterStudyRNAseq/commongarden/splicing_analysis_20260308/inputs/refs/annotation.rmats.gtf \
#     --rf -o /storage/stacylab/WaterStudyRNAseq/splicing_waterstudy/taxonG/stringtie/per_sample_gtf/S15.stringtie.gtf -p 22
#
# Merged GTF (StringTie 2.2.3, header of merged.f0.01.gtf; output was written as
# merged3.gtf and later renamed to merged.f0.01.gtf):
#   stringtie --merge -G Metrosideros_polymrpha_incana.scaffolded.gtf -c 2 -T 1 -f 0.01 -m 200 -o merged3.gtf gtf_list.txt
# (guide GTF converted from the .gff with AGAT v1.7.0, GTFrelax output, env agat_env:
#   agat_convert_sp_gff2gtf.pl --gff Metrosideros_polymrpha_incana.scaffolded.gff -o Metrosideros_polymrpha_incana.scaffolded.gtf)
#
# Per-sample gffcompare v0.12.10, with merged.gtf -> merged.f0.01.gtf at the time
# (command line in line 2 of support/<S>.stats):
#   32 samples (taxonG, taxonI, S11), run as jupyter-eynizade:
#     gffcompare -r merged.gtf -o support/S1 taxonI/stringtie/per_sample_gtf/S1.stringtie.gtf
#   31 samples (rest of taxonL, taxonT), run as mirahan on 2026-10-10:
#     /srv/jupyter-homes/jupyter-eynizade/.conda/envs/stringtie22/bin/gffcompare -r merged.gtf \
#       -o support/S12 taxonT/stringtie/per_sample_gtf/S12.stringtie.gtf
#
# Support table + filter: originally done with ad-hoc awk one-liners (MIN=3,
# class code "=" only, taxon from path). Steps 2-3 of this script were written
# afterwards and reproduce support/merged.transcript_support.tsv and
# merged.support3pertaxon.gtf byte-for-byte from the existing .tmap files.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
. "${CONFIG:-pipeline/config.sh}"

REF=${REF:-merged.f0.01.gtf}                   # unfiltered merged GTF (merged.gtf now -> filtered output)
LIST=${LIST:-gtf_list.txt}                     # one per-sample GTF path per line
MIN=${MIN:-0}                                  # transcript-level filter: min samples within one taxon (0 = off)
JOBS=${JOBS:-6}                                # parallel gffcompare runs
FORCE=${FORCE:-0}                              # 1 = rerun gffcompare even if outputs exist
OUTDIR=${OUTDIR:-support}                      # gffcompare -o prefix directory
SUPPORT_TABLE=${SUPPORT_TABLE:-1}              # 1 = write the per-transcript sample-support table even when MIN=0
if [[ $MIN -gt 0 ]]; then TABLE=${TABLE:-support/merged.transcript_support.tsv}; else TABLE=${TABLE:-support/$(basename "${REF%.gtf}").transcript_support.tsv}; fi
if [[ $MIN -gt 0 ]]; then OUT_GTF=${OUT_GTF:-merged.support${MIN}pertaxon.gtf}; else OUT_GTF=$REF; fi
GTF_SUFFIX=${GTF_SUFFIX:-.stringtie.gtf}           # per-sample GTF name = <sample><GTF_SUFFIX>
TAXON_FIELD=${TAXON_FIELD:-1}                  # "/"-separated path component in $LIST holding the taxon
GFFCOMPARE=${GFFCOMPARE:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/stringtie22/bin/gffcompare}
# intron filters (steps 4-6); MOTIF_FILTER=0 and INTRON_MIN_SAMPLES=0 skip them
MOTIF_FILTER=${MOTIF_FILTER:-1}                # 1 = require GT-AG / GC-AG / AT-AC introns
INTRON_MIN_SAMPLES=${INTRON_MIN_SAMPLES:-3}    # intron needs >= this many samples ...
INTRON_MIN_READS=${INTRON_MIN_READS:-3}        # ... each with >= this many unique spliced reads
MAX_INTRON=${MAX_INTRON:-20000}                # max intron length in bp (0 = no limit)
GENOME=${GENOME:-Metrosideros_polymrpha_incana.scaffolded.fa}
FAI=${FAI:-support/$(basename "$GENOME").fai}  # samtools faidx index (kept out of the FASTA's dir)
BAM_DIR=${BAM_DIR:-BAMs/star_alignmentuntrimmed_results}
BAM_SUFFIX=${BAM_SUFFIX:-.Aligned.sortedByCoord.out.bam}
JUNC_DIR=${JUNC_DIR:-support/junction_counts}  # per-sample junction counts (cache)
ANCHOR_MIN=${ANCHOR_MIN:-12}                   # overhang rule: min aligned bp each side (0 = off)
[[ $ANCHOR_MIN == 12 ]] && _ad=support/junction_counts_anchor || _ad=support/junction_counts_anchor$ANCHOR_MIN
ANCHOR_DIR=${ANCHOR_DIR:-$_ad}                 # per-sample overhang counts (cache, per ANCHOR_MIN)
SAMTOOLS=${SAMTOOLS:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/rmats420/bin/samtools}
sfx=""
[[ $MOTIF_FILTER == 1 ]] && sfx+=.canonical
[[ $INTRON_MIN_SAMPLES -gt 0 ]] && sfx+=".sj${INTRON_MIN_SAMPLES}x${INTRON_MIN_READS}"
[[ $INTRON_MIN_SAMPLES -gt 0 && $ANCHOR_MIN -gt 0 ]] && sfx+=".anc${ANCHOR_MIN}"
[[ $MAX_INTRON -gt 0 ]] && sfx+=".max$((MAX_INTRON / 1000))kb"
true
FINAL_GTF=${FINAL_GTF:-${OUT_GTF%.gtf}$sfx.gtf}
INTRON_TABLE=${INTRON_TABLE:-support/$(basename "${FINAL_GTF%.gtf}").intron_qc.tsv}
# read-through filter (steps 7-9): off unless READTHROUGH_FILTER=1 or ONLY_READTHROUGH=1
READTHROUGH_FILTER=${READTHROUGH_FILTER:-0}
[[ ${ONLY_READTHROUGH:-0} == 1 ]] && READTHROUGH_FILTER=1
MAKER_GTF=${MAKER_GTF:-Metrosideros_polymrpha_incana.scaffolded.gtf}   # gene models defining "genes"
CONDITIONS=${CONDITIONS:-support/samples.tsv}   # sample <tab> group <tab> condition (= SAMPLES in config.sh)
RT_MAX_RATIO=${RT_MAX_RATIO:-0.2}     # read-through = bridge reads < this x gene-internal junction reads
RT_FOLD=${RT_FOLD:-3}                 # stress-elevated: stress ratio >= RT_MAX_RATIO and >= RT_FOLD x control
RT_MIN_READS=${RT_MIN_READS:-10}      # ... with >= this many bridging reads in the stress group
QC_DIR=${QC_DIR:-pipeline/lib}        # bridging_ratio.py, readthrough_*.py, split_disjoint_genes.py
PYTHON=${PYTHON:-/srv/jupyter-homes/jupyter-eynizade/.conda/envs/grase/bin/python}

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
mkdir -p "$OUTDIR" "$(dirname "$TABLE")"

# tmap path gffcompare writes for a given input GTF
tmap_of() { local d b; d=$(dirname "$1"); b=$(basename "$1"); echo "$d/${b%"$GTF_SUFFIX"}.$b.tmap"; }
gtf_summary() {
  awk -F'\t' '$3=="transcript"{t++; match($9,/gene_id "[^"]+"/); g[substr($9,RSTART,RLENGTH)]} $3=="exon"{e++}
    END{printf "kept %d transcripts, %d exons, %d genes", t, e, length(g)}' "$1"
}
# sample names / taxa from the GTF list
samples_of_list() { awk -F/ -v x="$GTF_SUFFIX" '{s=$NF; if(substr(s,length(s)-length(x)+1)==x) s=substr(s,1,length(s)-length(x)); print s}' "$LIST"; }

if [[ ${ONLY_READTHROUGH:-0} != 1 ]]; then   # ONLY_READTHROUGH=1: steps 7-9 on an existing $FINAL_GTF
if [[ $MIN -gt 0 || $SUPPORT_TABLE == 1 ]]; then   # steps 1-2: support table; step 3 only with MIN > 0
# --- 1. per-sample gffcompare -------------------------------------------------
run_one() {
  local f=$1 s
  s=$(basename "$f" "$GTF_SUFFIX")
  if [[ $FORCE != 1 && -s $(tmap_of "$f") ]] && grep -q 'consensus transcripts written' "$OUTDIR/$s.stats" 2>/dev/null; then
    echo "skip $s"; return
  fi
  "$GFFCOMPARE" -r "$REF" -o "$OUTDIR/$s" "$f" 2>>"$OUTDIR/$s.gffcompare.log" \
    && [[ -s $(tmap_of "$f") ]] && echo "ok   $s" || { echo "FAIL $s (see $OUTDIR/$s.gffcompare.log)"; return 1; }
}
export -f run_one tmap_of
export REF OUTDIR FORCE GFFCOMPARE GTF_SUFFIX
echo "== gffcompare vs $REF"
xargs -a "$LIST" -P "$JOBS" -n 1 bash -c 'run_one "$0"' | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'

# --- 2. support table ---------------------------------------------------------
echo "== support table -> $TABLE"
# exact matches: ref_transcript <tab> taxon <tab> sample (deduplicated)
while read -r f; do
  awk -v t="$(cut -d/ -f"$TAXON_FIELD" <<<"$f")" -v s="$(basename "$f" "$GTF_SUFFIX")" -F'\t' \
    'FNR>1 && $3=="=" {print $2"\t"t"\t"s}' "$(tmap_of "$f")"
done < "$LIST" | sort -u > "$tmpdir/eq.tsv"

# transcript -> gene from the merged GTF
awk -F'\t' '$3=="transcript"{match($9,/gene_id "[^"]+"/); g=substr($9,RSTART+9,RLENGTH-10)
  match($9,/transcript_id "[^"]+"/); print substr($9,RSTART+15,RLENGTH-16)"\t"g}' "$REF" > "$tmpdir/t2g.tsv"

taxa=$(cut -d/ -f"$TAXON_FIELD" "$LIST" | sort -u | paste -sd' ')
awk -F'\t' -v OFS='\t' -v taxa="$taxa" -v min="$MIN" '
  NR==FNR {c[$1 SUBSEP $2]++; next}
  FNR==1  {nt=split(taxa,T," "); h="transcript_id\tgene_id\ttotal"; for(i=1;i<=nt;i++) h=h"\t"T[i]
           print h, "max_taxon_support", "keep"}
  { tot=0; m=0; cols=""
    for(i=1;i<=nt;i++){v=c[$1 SUBSEP T[i]]+0; tot+=v; if(v>m)m=v; cols=cols"\t"v}
    print $1"\t"$2"\t"tot cols, m, (m>=min?1:0) }' "$tmpdir/eq.tsv" "$tmpdir/t2g.tsv" > "$TABLE"

# --- 3. filtered GTF ----------------------------------------------------------
if [[ $MIN -gt 0 ]]; then
echo "== filtered GTF -> $OUT_GTF"
awk -F'\t' 'NR==FNR {if(FNR>1 && $NF==1) k[$1]=1; next}
  /^#/ {print; next}
  {match($9,/transcript_id "[^"]+"/); if(substr($9,RSTART+15,RLENGTH-16) in k) print}' "$TABLE" "$REF" > "$OUT_GTF"

echo "$(gtf_summary "$OUT_GTF") (MIN=$MIN per taxon)"
else
  echo "== transcript-level filter off (MIN=0): table only; intron filters on $REF"
fi
else
  echo "== transcript-level filter off (MIN=0), no support table: intron filters on $REF"
fi

# --- 4-6. intron filters on $OUT_GTF ------------------------------------------
# 4. splice motif   : every intron must be GT-AG, GC-AG or AT-AC on the transcript strand
# 5. read support   : every intron must have >= INTRON_MIN_READS uniquely mapped spliced reads
#                     (MAPQ 255, primary) in >= INTRON_MIN_SAMPLES samples, counted from the
#                     BAMs (NOT STAR SJ.out.tab, which omits junctions that fail STAR's
#                     outSJfilter* rules even though the reads are in the BAM)
# 6. write $FINAL_GTF (single-exon transcripts have no introns and always pass)
if [[ $MOTIF_FILTER == 1 || $INTRON_MIN_SAMPLES -gt 0 || $MAX_INTRON -gt 0 ]]; then
mkdir -p "$(dirname "$INTRON_TABLE")"

# unique introns of $OUT_GTF: chrom strand start end n_tx
awk -F'\t' -v OFS='\t' '$3=="exon"{match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
    if(t==pt) n[$1 OFS $7 OFS pe+1 OFS $4-1]++; pt=t; pe=$5}
  END{for(k in n) print k, n[k]}' "$OUT_GTF" | sort -k1,1 -k3,3n -k4,4n -k2,2 > "$tmpdir/introns.tsv"
echo "== intron QC: $(wc -l < "$tmpdir/introns.tsv") unique introns"

# 4. motifs: donor = first 2 intron bases, acceptor = last 2 (+ strand), reverse-complemented for -
echo "   splice motifs from $GENOME"
awk -F'\t' '{print $1":"$3"-"$3+1; print $1":"$4-1"-"$4}' "$tmpdir/introns.tsv" > "$tmpdir/regions.txt"
"$SAMTOOLS" faidx --fai-idx "$FAI" "$GENOME" -r "$tmpdir/regions.txt" -n 1000 \
  | awk '/^>/{next} {print toupper($0)}' | paste - - > "$tmpdir/dinuc.tsv"
paste "$tmpdir/introns.tsv" "$tmpdir/dinuc.tsv" | awk -F'\t' -v OFS='\t' '
  function rc(s,  i,o,c){o=""; for(i=length(s);i>0;i--){c=substr(s,i,1); o=o (c=="A"?"T":c=="C"?"G":c=="G"?"C":c=="T"?"A":"N")}; return o}
  {m=($2=="-") ? rc($7)"-"rc($6) : $6"-"$7
   print $1,$2,$3,$4,$5,m,(m=="GT-AG"||m=="GC-AG"||m=="AT-AC")?1:0}' > "$tmpdir/introns.motif.tsv"

# 5. per-sample junction counts from the BAMs (cached in $JUNC_DIR; independent of the GTF)
if [[ $INTRON_MIN_SAMPLES -gt 0 ]]; then
  mkdir -p "$JUNC_DIR"
  count_junc() {
    local s=$1 bam=$BAM_DIR/$1$BAM_SUFFIX out=$JUNC_DIR/$1.junc.tsv
    [[ -s $out ]] && { echo "skip $s"; return; }
    [[ -s $bam ]] || { echo "FAIL $s (no BAM $bam)"; return 1; }
    # chrom <tab> first intron base <tab> last intron base <tab> unique spliced reads
    "$SAMTOOLS" view -@2 -q 255 -F 0x904 "$bam" 2>/dev/null | awk -F'\t' -v OFS='\t' '$6~/N/{
        p=$4; cig=$6
        while(match(cig,/[0-9]+[MIDNSHP=X]/)){l=substr(cig,RSTART,RLENGTH-1)+0; op=substr(cig,RSTART+RLENGTH-1,1)
          cig=substr(cig,RSTART+RLENGTH); if(op=="N") j[$3 OFS p OFS p+l-1]++; if(op~/[MDN=X]/) p+=l}}
      END{for(k in j) print k, j[k]}' > "$out.tmp" && mv "$out.tmp" "$out" && echo "ok   $s"
  }
  export -f count_junc
  export SAMTOOLS BAM_DIR BAM_SUFFIX JUNC_DIR
  echo "   junction reads from BAMs in $BAM_DIR (cache: $JUNC_DIR)"
  samples_of_list | xargs -P "$JOBS" -n 1 bash -c 'count_junc "$0"' \
    | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
  samples_of_list | while read -r s; do cat "$JUNC_DIR/$s.junc.tsv"; done \
    | awk -F'\t' -v OFS='\t' -v r="$INTRON_MIN_READS" '{k=$1 OFS $2 OFS $3; t[k]+=$4; if($4>=r) n[k]++}
        END{for(k in t) print k, t[k], n[k]+0}' > "$tmpdir/junc_all.tsv"
  # overhang counts: chrom, first, last intron base, unique reads, reads with >= ANCHOR_MIN bp each side
  if [[ $ANCHOR_MIN -gt 0 ]]; then
    mkdir -p "$ANCHOR_DIR"
    count_anchor() {
      local s=$1 out=$ANCHOR_DIR/$1.junc_anchor.tsv
      [[ -s $out ]] && { echo "skip $s"; return; }
      bash pipeline/lib/junction_anchor_counts.sh "$BAM_DIR/$s$BAM_SUFFIX" "$out" "$ANCHOR_MIN" "$SAMTOOLS" \
        && echo "ok   $s" || { echo "FAIL $s (overhang counts)"; return 1; }
    }
    export -f count_anchor
    export ANCHOR_DIR ANCHOR_MIN
    echo "   junction overhang counts (>= $ANCHOR_MIN bp each side; cache: $ANCHOR_DIR)"
    samples_of_list | xargs -P "$JOBS" -n 1 bash -c 'count_anchor "$0"' \
      | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
    # samples with >= INTRON_MIN_READS reads AND >= 1 anchored read
    samples_of_list | while read -r s; do cat "$ANCHOR_DIR/$s.junc_anchor.tsv"; done \
      | awk -F'\t' -v OFS='\t' -v r="$INTRON_MIN_READS" '$4>=r && $5>=1 {n[$1 OFS $2 OFS $3]++}
          END{for(k in n) print k, n[k]}' > "$tmpdir/junc_anchor.tsv"
  fi
fi
[[ -e $tmpdir/junc_anchor.tsv ]] || : > "$tmpdir/junc_anchor.tsv"   # overhang rule off
# introns next to an internal exon < ANCHOR_MIN bp (micro-exon): exempt from the overhang rule
awk -F'\t' -v OFS='\t' -v amin="$ANCHOR_MIN" '
  function flush(  i) { for (i = 2; i < n; i++) if (e[i] - s[i] + 1 < amin) {
      print c[i], st[i], e[i-1] + 1, s[i] - 1; print c[i], st[i], e[i] + 1, s[i+1] - 1 } }
  $3=="exon" {match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
              if (t != pt) { flush(); n = 0 }; n++; c[n]=$1; st[n]=$7; s[n]=$4; e[n]=$5; pt=t}
  END {flush()}' "$OUT_GTF" | sort -u > "$tmpdir/micro.tsv"

awk -F'\t' -v OFS='\t' -v motif="$MOTIF_FILTER" -v ms="$INTRON_MIN_SAMPLES" -v r="$INTRON_MIN_READS" -v maxi="$MAX_INTRON" \
    -v amin="$ANCHOR_MIN" '
  FILENAME==ARGV[1] {t[$1 OFS $2 OFS $3]=$4; n[$1 OFS $2 OFS $3]=$5; next}
  FILENAME==ARGV[2] {na[$1 OFS $2 OFS $3]=$4; next}
  FILENAME==ARGV[3] {mx[$1 OFS $2 OFS $3 OFS $4]=1; next}
  FNR==1  {h="chrom\tstrand\tstart\tend\tlength\tn_tx\tmotif\tcanonical\tunique_reads_all_samples\tn_samples_ge" r "reads"
           if (amin > 0 && ms > 0) h=h "\tn_samples_ge" r "reads_anchored" amin "\tmicroexon_exempt"
           print h "\tpass"}
  {k=$1 OFS $3 OFS $4; rd=(k in t)?t[k]:0; ns=(k in n)?n[k]:0; a=(k in na)?na[k]:0; x=(($1 OFS $2 OFS $3 OFS $4) in mx)
   pass=(motif!=1 || $7==1) && (ms<=0 || ns>=ms) && (maxi<=0 || $4-$3+1<=maxi) && (amin<=0 || ms<=0 || x || a>=ms)
   if (amin > 0 && ms > 0) print $1,$2,$3,$4,$4-$3+1,$5,$6,$7,rd,ns,a,x,pass?1:0
   else print $1,$2,$3,$4,$4-$3+1,$5,$6,$7,rd,ns,pass?1:0}' \
  <(cat "$tmpdir/junc_all.tsv" 2>/dev/null || true) "$tmpdir/junc_anchor.tsv" "$tmpdir/micro.tsv" \
  "$tmpdir/introns.motif.tsv" > "$INTRON_TABLE"

# 6. drop transcripts with any failing intron, and transcripts without a strand
#    (single-exon, strand "."; HTSeq/dexseq_prepare_annotation.py in GrASE fails on them)
echo "== intron-filtered GTF -> $FINAL_GTF"
awk -F'\t' -v OFS='\t' 'NR==FNR {if(FNR>1 && $NF==0) bad[$1 OFS $2 OFS $3 OFS $4]=1; next}
  $3=="exon" {match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
              if($7!="+" && $7!="-") {drop[t]=1; nostrand[t]=1}
              if(t==pt && (($1 OFS $7 OFS pe+1 OFS $4-1) in bad)) drop[t]=1; pt=t; pe=$5}
  END{for(t in drop) print t; printf "%d\n", length(nostrand) > "/dev/stderr"}' "$INTRON_TABLE" "$OUT_GTF" \
  > "$tmpdir/drop.txt" 2> "$tmpdir/nostrand.txt"
awk -F'\t' 'NR==FNR {d[$1]=1; next} /^#/ {print; next}
  {match($9,/transcript_id "[^"]+"/); if(!(substr($9,RSTART+15,RLENGTH-16) in d)) print}' "$tmpdir/drop.txt" "$OUT_GTF" > "$FINAL_GTF"

awk -F'\t' -v maxi="$MAX_INTRON" 'NR>1{n++; if(!$8)nc++; if(maxi>0 && $5>maxi)lg++; if($NF==0)f++} END{printf "   introns: %d, non-canonical: %d, longer than max: %d, failing filters: %d\n", n, nc, lg, f}' "$INTRON_TABLE"
echo "   dropped $(wc -l < "$tmpdir/drop.txt") transcripts with a failing intron or no strand ($(cat "$tmpdir/nostrand.txt") without strand)"
echo "$(gtf_summary "$FINAL_GTF") (motif filter=$MOTIF_FILTER, intron support >=$INTRON_MIN_READS reads in >=$INTRON_MIN_SAMPLES samples, overhang >=${ANCHOR_MIN} bp, max intron ${MAX_INTRON} bp)"
else
  FINAL_GTF=$OUT_GTF   # no intron filters: read-through filter works on $OUT_GTF
fi

fi

# --- 7-9. read-through filter on $FINAL_GTF -----------------------------------
# Read-through: transcription continues past the end of one gene into the next gene on
# the same strand, and splicing joins an exon of gene A to an exon of gene B. Seen as a
# junction with donor in one $MAKER_GTF gene and acceptor in another ("bridging").
# 7. ratio = bridge reads / min(median reads of junctions inside gene A, inside gene B)
#    (unique reads summed over all samples, from $JUNC_DIR). Low ratio = minor
#    read-through; ratio ~1 = MAKER probably split one gene into two models.
# 8. bridging junctions with ratio < RT_MAX_RATIO are dropped unless elevated under
#    stress in some taxon (ratio >= RT_MAX_RATIO, >= RT_FOLD x control, >= RT_MIN_READS
#    bridging reads in that group); edit the decision column of the table to override.
# 9. drop transcripts using a dropped junction, then give each disjoint piece of a
#    gene_id its own id (<gene_id>.p2, ...) -> $NORT_GTF
[[ $READTHROUGH_FILTER == 1 ]] || exit 0
[[ $INTRON_MIN_SAMPLES -gt 0 ]] || { echo "read-through filter needs junction counts (INTRON_MIN_SAMPLES > 0)"; exit 1; }
NORT_GTF=${NORT_GTF:-${FINAL_GTF%.gtf}.nort.gtf}
RT_PREFIX=${RT_PREFIX:-support/$(basename "${FINAL_GTF%.gtf}")}
[[ -s $CONDITIONS ]] || { echo "missing sample sheet $CONDITIONS (SAMPLES in pipeline/config.sh)"; exit 1; }
echo "== read-through filter -> $NORT_GTF"
"$PYTHON" -I "$QC_DIR/bridging_ratio.py" "$FINAL_GTF" "$MAKER_GTF" "$JUNC_DIR" "$RT_PREFIX.bridging_ratio.tsv" | sed 's/^/   /'
"$PYTHON" -I "$QC_DIR/readthrough_by_condition.py" "$FINAL_GTF" "$MAKER_GTF" "$RT_PREFIX.bridging_ratio.tsv" \
  "$CONDITIONS" "$JUNC_DIR" "$RT_PREFIX.readthrough_by_condition.tsv" "$RT_MAX_RATIO" "$RT_FOLD" "$RT_MIN_READS" | sed 's/^/   /'
"$PYTHON" -I "$QC_DIR/readthrough_filter.py" "$FINAL_GTF" "$RT_PREFIX.readthrough_by_condition.tsv" "$tmpdir/nort.gtf" | sed 's/^/   /'
"$PYTHON" -I "$QC_DIR/split_disjoint_genes.py" "$tmpdir/nort.gtf" "$NORT_GTF" | sed 's/^/   /'
echo "$(gtf_summary "$NORT_GTF") (read-through filter: ratio < $RT_MAX_RATIO dropped unless stress-elevated)"
