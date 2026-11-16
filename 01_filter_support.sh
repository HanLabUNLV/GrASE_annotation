#!/usr/bin/env bash
# Step 01: junction-evidence filters for a GTF. Not run on its own: step 03 calls it for the junction
# rules (stage B) and for the junction-level read-through filter (stage D, ONLY_READTHROUGH=1).
#
# Intron filters (steps 1-3 below): a transcript is kept if every intron
#   - has a canonical splice motif (GT-AG/GC-AG/AT-AC on the transcript strand; MOTIF_FILTER=1),
#   - is <= MAX_INTRON bp (20 kb),
#   - has >= INTRON_MIN_READS (3) uniquely mapped spliced reads in >= INTRON_MIN_SAMPLES (3) samples,
#     counted from the BAMs (MAPQ 255, primary; not STAR's SJ.out.tab, which omits junctions failing
#     STAR's outSJfilter* rules although the reads are in the BAM),
#   - overhang (ANCHOR_MIN=12; 0 = off): each of those samples has >= 1 read with >= ANCHOR_MIN bp
#     aligned on both sides of the intron (STAR's SJ.out.tab overhang rule; catches untrimmed read
#     ends forced into a splice, e.g. 90M65N3M7S). Introns next to an internal exon shorter than
#     ANCHOR_MIN are exempt (no read can anchor there).
#   Transcripts without a strand are dropped too (GrASE's DEXSeq preparation fails on them).
#   Per-sample counts are cached: $JUNC_DIR (unique reads), $ANCHOR_DIR (overhang).
#   -> $FINAL_GTF, intron table $INTRON_TABLE (motif, reads, samples, anchored samples, pass)
#
# Read-through filter (steps 4-6; READTHROUGH_FILTER=1 or ONLY_READTHROUGH=1), with reference gene
# models $MAKER_GTF: a junction with donor in one reference gene and acceptor in another ("bridging")
#   4. ratio = bridge reads / min(median reads of the junctions inside each gene), all samples;
#      low = minor read-through, ~1 or more = the reference probably split one gene in two
#   5. bridging junctions with ratio < RT_MAX_RATIO (0.2) are dropped unless elevated under a
#      non-control condition in some group (ratio >= RT_MAX_RATIO, >= RT_FOLD x control, >=
#      RT_MIN_READS bridging reads; step 03 disables this by default)
#   6. drop transcripts using a dropped junction; each disjoint piece of a gene_id gets its own id
#   -> $NORT_GTF, tables $RT_PREFIX.{bridging_ratio,readthrough_by_condition}.tsv
#
# Samples: $SAMPLES (pipeline/config.sh), BAMs <BAM_DIR>/<sample><BAM_SUFFIX>.
# Usage (as step 03 calls it):
#   REF=<in.gtf> FINAL_GTF=<out.gtf> INTRON_TABLE=<tsv> pipeline/01_filter_support.sh
#   ONLY_READTHROUGH=1 FINAL_GTF=<in.gtf> NORT_GTF=<out.gtf> RT_PREFIX=<prefix> pipeline/01_filter_support.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
. "${CONFIG:-pipeline/config.sh}"

MOTIF_FILTER=${MOTIF_FILTER:-1}                # 1 = require GT-AG / GC-AG / AT-AC introns
INTRON_MIN_SAMPLES=${INTRON_MIN_SAMPLES:-3}    # intron needs >= this many samples ...
INTRON_MIN_READS=${INTRON_MIN_READS:-3}        # ... each with >= this many unique spliced reads
MAX_INTRON=${MAX_INTRON:-20000}                # max intron length in bp (0 = no limit)
ANCHOR_MIN=${ANCHOR_MIN:-12}                   # overhang rule: min aligned bp each side (0 = off)
FAI=${FAI:-support/$(basename "$GENOME").fai}  # samtools faidx index (kept out of the FASTA's dir)
JUNC_DIR=${JUNC_DIR:-support/junction_counts}  # per-sample junction counts (cache)
[[ $ANCHOR_MIN == 12 ]] && _ad=support/junction_counts_anchor || _ad=support/junction_counts_anchor$ANCHOR_MIN
ANCHOR_DIR=${ANCHOR_DIR:-$_ad}                 # per-sample overhang counts (cache, per ANCHOR_MIN)
READTHROUGH_FILTER=${READTHROUGH_FILTER:-0}
[[ ${ONLY_READTHROUGH:-0} == 1 ]] && READTHROUGH_FILTER=1
RT_MAX_RATIO=${RT_MAX_RATIO:-0.2}     # read-through = bridge reads < this x gene-internal junction reads
RT_FOLD=${RT_FOLD:-3}                 # condition-elevated: ratio >= RT_MAX_RATIO and >= RT_FOLD x control
RT_MIN_READS=${RT_MIN_READS:-10}      # ... with >= this many bridging reads in that condition
QC_DIR=pipeline/lib
if [[ ${ONLY_READTHROUGH:-0} != 1 ]]; then
  [[ -n ${REF:-} && -s $REF ]] || { echo "01: REF=<gtf> required"; exit 1; }
  FINAL_GTF=${FINAL_GTF:-${REF%.gtf}.filtered.gtf}
  INTRON_TABLE=${INTRON_TABLE:-${FINAL_GTF%.gtf}.intron_qc.tsv}
fi
[[ -n ${FINAL_GTF:-} ]] || { echo "01: FINAL_GTF=<gtf> required"; exit 1; }

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
gtf_summary() {
  awk -F'\t' '$3=="transcript"{t++; match($9,/gene_id "[^"]+"/); g[substr($9,RSTART,RLENGTH)]} $3=="exon"{e++}
    END{printf "kept %d transcripts, %d exons, %d genes", t, e, length(g)}' "$1"
}
samples() { cut -f1 "$SAMPLES"; }

if [[ ${ONLY_READTHROUGH:-0} != 1 ]]; then
mkdir -p "$(dirname "$INTRON_TABLE")"

# unique introns of $REF: chrom strand start end n_tx
awk -F'\t' -v OFS='\t' '$3=="exon"{match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
    if(t==pt) n[$1 OFS $7 OFS pe+1 OFS $4-1]++; pt=t; pe=$5}
  END{for(k in n) print k, n[k]}' "$REF" | sort -k1,1 -k3,3n -k4,4n -k2,2 > "$tmpdir/introns.tsv"
echo "== intron QC: $(wc -l < "$tmpdir/introns.tsv") unique introns"

# 1. motifs: donor = first 2 intron bases, acceptor = last 2 (+ strand), reverse-complemented for -
echo "   splice motifs from $GENOME"
awk -F'\t' '{print $1":"$3"-"$3+1; print $1":"$4-1"-"$4}' "$tmpdir/introns.tsv" > "$tmpdir/regions.txt"
"$SAMTOOLS" faidx --fai-idx "$FAI" "$GENOME" -r "$tmpdir/regions.txt" -n 1000 \
  | awk '/^>/{next} {print toupper($0)}' | paste - - > "$tmpdir/dinuc.tsv"
paste "$tmpdir/introns.tsv" "$tmpdir/dinuc.tsv" | awk -F'\t' -v OFS='\t' '
  function rc(s,  i,o,c){o=""; for(i=length(s);i>0;i--){c=substr(s,i,1); o=o (c=="A"?"T":c=="C"?"G":c=="G"?"C":c=="T"?"A":"N")}; return o}
  {m=($2=="-") ? rc($7)"-"rc($6) : $6"-"$7
   print $1,$2,$3,$4,$5,m,(m=="GT-AG"||m=="GC-AG"||m=="AT-AC")?1:0}' > "$tmpdir/introns.motif.tsv"

# 2. per-sample junction counts from the BAMs (cached; independent of the GTF)
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
  export JUNC_DIR
  echo "   junction reads from BAMs in $BAM_DIR (cache: $JUNC_DIR)"
  samples | xargs -P "$JOBS" -n 1 bash -c 'count_junc "$0"' \
    | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
  samples | while read -r s; do cat "$JUNC_DIR/$s.junc.tsv"; done \
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
    samples | xargs -P "$JOBS" -n 1 bash -c 'count_anchor "$0"' \
      | sort | awk '{n[$1]++} $1=="FAIL"{print} END{for(k in n) printf "  %s: %d\n",k,n[k]}'
    # samples with >= INTRON_MIN_READS reads AND >= 1 anchored read
    samples | while read -r s; do cat "$ANCHOR_DIR/$s.junc_anchor.tsv"; done \
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
  END {flush()}' "$REF" | sort -u > "$tmpdir/micro.tsv"

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

# 3. drop transcripts with any failing intron, and transcripts without a strand
echo "== intron-filtered GTF -> $FINAL_GTF"
awk -F'\t' -v OFS='\t' 'NR==FNR {if(FNR>1 && $NF==0) bad[$1 OFS $2 OFS $3 OFS $4]=1; next}
  $3=="exon" {match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
              if($7!="+" && $7!="-") {drop[t]=1; nostrand[t]=1}
              if(t==pt && (($1 OFS $7 OFS pe+1 OFS $4-1) in bad)) drop[t]=1; pt=t; pe=$5}
  END{for(t in drop) print t; printf "%d\n", length(nostrand) > "/dev/stderr"}' "$INTRON_TABLE" "$REF" \
  > "$tmpdir/drop.txt" 2> "$tmpdir/nostrand.txt"
awk -F'\t' 'NR==FNR {d[$1]=1; next} /^#/ {print; next}
  {match($9,/transcript_id "[^"]+"/); if(!(substr($9,RSTART+15,RLENGTH-16) in d)) print}' "$tmpdir/drop.txt" "$REF" > "$FINAL_GTF"

awk -F'\t' -v maxi="$MAX_INTRON" 'NR>1{n++; if(!$8)nc++; if(maxi>0 && $5>maxi)lg++; if($NF==0)f++} END{printf "   introns: %d, non-canonical: %d, longer than max: %d, failing filters: %d\n", n, nc, lg, f}' "$INTRON_TABLE"
echo "   dropped $(wc -l < "$tmpdir/drop.txt") transcripts with a failing intron or no strand ($(cat "$tmpdir/nostrand.txt") without strand)"
echo "$(gtf_summary "$FINAL_GTF") (motif filter=$MOTIF_FILTER, intron support >=$INTRON_MIN_READS reads in >=$INTRON_MIN_SAMPLES samples, overhang >=${ANCHOR_MIN} bp, max intron ${MAX_INTRON} bp)"
fi

# --- 4-6. read-through filter on $FINAL_GTF -----------------------------------
[[ $READTHROUGH_FILTER == 1 ]] || exit 0
[[ -n $MAKER_GTF && -s $MAKER_GTF ]] || { echo "read-through filter needs reference gene models (REF_GTF)"; exit 1; }
NORT_GTF=${NORT_GTF:-${FINAL_GTF%.gtf}.nort.gtf}
RT_PREFIX=${RT_PREFIX:-${FINAL_GTF%.gtf}}
[[ -s $SAMPLES ]] || { echo "missing sample sheet $SAMPLES (pipeline/config.sh)"; exit 1; }
echo "== read-through filter -> $NORT_GTF"
"$PYTHON" -I "$QC_DIR/bridging_ratio.py" "$FINAL_GTF" "$MAKER_GTF" "$JUNC_DIR" "$RT_PREFIX.bridging_ratio.tsv" | sed 's/^/   /'
"$PYTHON" -I "$QC_DIR/readthrough_by_condition.py" "$FINAL_GTF" "$MAKER_GTF" "$RT_PREFIX.bridging_ratio.tsv" \
  "$SAMPLES" "$JUNC_DIR" "$RT_PREFIX.readthrough_by_condition.tsv" "$RT_MAX_RATIO" "$RT_FOLD" "$RT_MIN_READS" | sed 's/^/   /'
"$PYTHON" -I "$QC_DIR/readthrough_filter.py" "$FINAL_GTF" "$RT_PREFIX.readthrough_by_condition.tsv" "$tmpdir/nort.gtf" | sed 's/^/   /'
"$PYTHON" -I "$QC_DIR/split_disjoint_genes.py" "$tmpdir/nort.gtf" "$NORT_GTF" | sed 's/^/   /'
echo "$(gtf_summary "$NORT_GTF") (read-through filter: ratio < $RT_MAX_RATIO dropped unless condition-elevated)"
