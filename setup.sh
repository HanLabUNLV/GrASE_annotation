#!/usr/bin/env bash
# Setup and checks before a run (pipeline/run_all.sh calls it first). Reads pipeline/config.sh.
#   1. sample sheet: 3 tab-separated columns, unique samples, every group has $CONTROL and >= 1 other
#      condition; warns for conditions with < 3 samples (the junction rules need >= 3 samples)
#   2. BAMs: <BAM_DIR>/<sample><BAM_SUFFIX> exist; missing .bai indexes are built
#   3. genome, its samtools index (support/<genome>.fai), reference GTF (optional; needs CDS for the
#      homology check)
#   4. tools: each runs on this CPU (bioconda builds for newer CPUs fail with "Illegal instruction"),
#      R package grase, python with pandas/scipy, the GrASE clone
#   5. Swiss-Prot (HOMOLOGY_DB set): downloads UniProt Swiss-Prot and builds the BLAST database if
#      missing; release recorded in <db dir>/RELEASE.txt
# Usage: pipeline/setup.sh   (exit status 1 if a check fails)
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
. "${CONFIG:-pipeline/config.sh}"
fail=0
ok()   { echo "  ok    $*"; }
bad()  { echo "  FAIL  $*"; fail=1; }
warn() { echo "  warn  $*"; }

echo "== 1. sample sheet $SAMPLES"
if [[ ! -s $SAMPLES ]]; then bad "missing $SAMPLES (sample <tab> group <tab> condition)"
else
  awk -F'\t' 'NF!=3{print "line "NR": "NF" columns"; e=1} END{exit e}' "$SAMPLES" >/dev/null || bad "not 3 tab-separated columns per line"
  d=$(cut -f1 "$SAMPLES" | sort | uniq -d | head -3 | paste -sd,); [[ -z $d ]] || bad "duplicate samples: $d"
  for g in $(pipeline_groups); do
    pipeline_conditions "$g" | grep -qx "$CONTROL" || bad "group $g has no '$CONTROL' samples"
    [[ -n $(pipeline_treatments "$g") ]] || bad "group $g has no condition other than '$CONTROL'"
    for c in $(pipeline_conditions "$g"); do
      n=$(pipeline_samples "$g" "$c" | wc -l); [[ $n -ge 3 ]] || warn "$g/$c has $n samples (< 3)"
    done
  done
  ok "$(wc -l < "$SAMPLES") samples, groups: $(pipeline_groups | paste -sd' '); contrasts: $(for g in $(pipeline_groups); do for c in $(pipeline_treatments "$g"); do printf '%s/%s_vs_%s ' "$g" "$c" "$CONTROL"; done; done)"
fi

echo "== 2. BAMs in $BAM_DIR"
nb=0
for s in $(cut -f1 "$SAMPLES" 2>/dev/null); do
  b=$BAM_DIR/$s$BAM_SUFFIX
  if [[ ! -s $b ]]; then bad "missing $b"; continue; fi
  if [[ ! -s $b.bai && ! -s ${b%.bam}.bai ]]; then "$SAMTOOLS" index "$b" && echo "  index $b" || bad "cannot index $b"; fi
  nb=$((nb + 1))
done
ok "$nb BAMs"

echo "== 3. genome and reference"
if [[ -z $GENOME ]]; then bad "GENOME not set (pipeline.config.sh)"; elif [[ -s $GENOME ]]; then ok "genome $GENOME"; else bad "missing genome $GENOME"; fi
FAI=support/$(basename "$GENOME").fai; mkdir -p support
if [[ -s $GENOME && ! -s $FAI ]]; then "$SAMTOOLS" faidx --fai-idx "$FAI" "$GENOME" && ok "indexed -> $FAI" || bad "samtools faidx failed"; fi
if [[ -n $REF_GTF ]]; then
  if [[ -s $REF_GTF ]]; then
    ok "reference $REF_GTF: $(awk -F'\t' '$3=="gene"' "$REF_GTF" | wc -l) genes, $(awk -F'\t' '$3=="CDS"' "$REF_GTF" | wc -l) CDS lines"
  else bad "missing reference $REF_GTF (set REF_GTF= for reference-free)"; fi
else warn "no reference gene models: read-through and homology checks (step 03 D-D3) are skipped"; fi

echo "== 4. tools"
chk() { local name=$1; shift; local out; out=$("$@" 2>&1); local st=$?
  if [[ $st -eq 0 ]]; then ok "$name: $(echo "$out" | grep -m1 . | cut -c1-60)"
  elif [[ $st -eq 132 ]] || grep -q 'Illegal instruction' <<<"$out"; then bad "$name: Illegal instruction (build needs a newer CPU; install an older build)"
  else bad "$name: $(echo "$out" | tail -1 | cut -c1-100)"; fi; }
chk scallop2 "$SCALLOP2" --version
chk gffcompare "$GFFCOMPARE" --version
chk samtools "$SAMTOOLS" --version
chk rMATS "$RMATS_ENV/bin/python" "$RMATS_ENV/rMATS/rmats.py" --version
chk "R grase" "$GRASE_ENV/bin/Rscript" -e 'suppressMessages(library(grase)); cat("grase", as.character(packageVersion("grase")), "padj default", formals(grase::adjust_pvalues)$method, "\n")'
chk "python pandas/scipy" "$PYTHON" -c 'import pandas, scipy; print("pandas", pandas.__version__, "scipy", scipy.__version__)'
[[ -s $GRASE_DIR/Rpkg/scripts/exontest.R ]] && ok "GrASE clone $GRASE_DIR" || bad "missing GrASE clone $GRASE_DIR (github.com/HanLabUNLV/GrASE)"
[[ -n $HOMOLOGY_DB ]] && { chk blastp "$BLASTP" -version; chk makeblastdb "$MAKEBLASTDB" -version; }

echo "== 5. Swiss-Prot"
if [[ -z $HOMOLOGY_DB ]]; then warn "HOMOLOGY_DB empty: homology check (step 03 D3) skipped"
elif [[ -s $HOMOLOGY_DB.pin || -s $HOMOLOGY_DB.00.pin || -s $HOMOLOGY_DB.pal ]]; then ok "BLAST db $HOMOLOGY_DB ($(head -2 "$(dirname "$HOMOLOGY_DB")/RELEASE.txt" 2>/dev/null | tail -1))"
else
  D=$(dirname "$HOMOLOGY_DB"); mkdir -p "$D"; U=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete
  echo "  downloading Swiss-Prot -> $D"
  if curl -sfL --max-time 1800 -o "$D/uniprot_sprot.fasta.gz" "$U/uniprot_sprot.fasta.gz" \
     && gunzip -f "$D/uniprot_sprot.fasta.gz" \
     && "$MAKEBLASTDB" -in "$D/uniprot_sprot.fasta" -dbtype prot -out "$HOMOLOGY_DB" -parse_seqids > "$D/makeblastdb.log" 2>&1; then
    curl -sfL --max-time 60 "$U/reldate.txt" | head -3 > "$D/RELEASE.txt"
    ok "built $HOMOLOGY_DB ($(sed -n 2p "$D/RELEASE.txt"))"
  else bad "Swiss-Prot download or makeblastdb failed (set HOMOLOGY_DB= to skip the homology check)"; fi
fi

echo "== disk: $(df -h . | awk 'NR==2{print $4" free"}')"
[[ $fail -eq 0 ]] && echo "== setup OK" || echo "== setup FAILED"
exit $fail
