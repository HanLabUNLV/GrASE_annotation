#!/usr/bin/env bash
# Step 03 (since 2026-11-13): annotation from the per-sample Scallop2 assemblies of step 02, by an
# evidence-based union instead of a merger. stringtie --merge and TACO keep only the dominant
# isoforms of a locus (here 45-57% of well-supported junctions vs 73% in the per-sample
# assemblies), which is what gene-level expression wants and the opposite of what splicing needs.
#   A. gffcompare -i: non-redundant union of the 63 assemblies (identical intron chains collapsed);
#      num_samples = samples whose assembly has exactly that intron chain  -> union/union.combined.gtf
#   B. junction rules (step 01 engine, 01_filter_support.sh): every intron canonical (GT-AG/GC-AG/
#      AT-AC), <= 20 kb, >= 3 unique reads in >= 3 samples, each of those with >= 1 read anchored
#      >= 12 bp on both sides (micro-exon introns exempt); single-exon and unstranded transcripts
#      dropped                                                      -> union/union.junction.gtf
#   C. lib/select_transcripts.py:
#      - paralog filter: drop transcripts using a junction whose reads are mostly multi-mapped
#        (multi / (unique + multi) > MAX_MULTI_FRAC, summed over samples; multi-mapped primary
#        reads from support/junction_counts_multi, lib/junction_multimap_counts.sh)
#      - keep transcripts assembled in >= MIN_SAMPLES samples, plus rarer ones that carry a
#        junction no kept transcript has (no junction that passed B is lost) -> union/union.selected.gtf
#   D. read-through filter (step 01 stages 7-9, as step 06): junctions joining two reference (MAKER)
#      genes with < 0.2 x the genes' own junction reads are dropped (stress-elevated ones are kept
#      only with RT_STRESS=1; off by default)
#                                                                   -> union/union.nort.gtf
#   D2. lib/readthrough_tx_filter.py: transcripts whose exons touch >= 2 reference genes (however
#      linked: exon across the intergenic gap, junction into an unannotated exon) are dropped when
#      assembled in < RT_TX_RATIO (0.5) x the samples of the stronger touched gene's best own
#      transcript (RT_TX_COMPARE=max; min = weaker gene)
#      (stress exemption only with RT_STRESS=1); dominant bridges (reference split one gene) are kept
#                                                                   -> union/union.nort_tx.gtf
#                                                                      (+ union.readthrough_tx.tsv)
#   D3. lib/homology_bridges.py: for each pair of reference genes a kept bridge joins, the two
#      reference proteins are searched against Swiss-Prot (local blastp; HOMOLOGY_DB, release in
#      support/swissprot/RELEASE.txt). Halves of one protein (MAKER split one gene) -> keep;
#      unrelated proteins (fusion/read-through) or the same region of one protein (paralog
#      chimera) -> drop the bridging transcripts; no hit -> keep  -> union/union.nort_hom.gtf
#                                                                    (+ union.bridge_homology.tsv)
#   E. lib/assign_genes.py: gene_id by exon overlap on the same strand (not span overlap), but two
#      groups touching different reference genes are never joined by exon overlap alone (UTR
#      extensions of neighbours); only a transcript touching both (kept by D2) joins them;
#      ref_gene_id = overlapping reference genes ("novel" if none)  -> $FINAL
# Restartable: stages whose output exists are skipped (delete it to redo).
#
# Usage: [ASM=assembly_unique MIN_SAMPLES=5 MAX_MULTI_FRAC=0.5 RT_TX_RATIO=0.5 RT_TX_COMPARE=max RT_STRESS=0 REF_GTF=<maker.gtf>] pipeline/03_union_filter.sh
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."   # project directory
. "${CONFIG:-pipeline/config.sh}"
ASM=${ASM:-assembly_unique}
U=$ASM/union
FINAL=${FINAL:-$ASM/annotation.final.gtf}
MIN_SAMPLES=${MIN_SAMPLES:-5}
MAX_MULTI_FRAC=${MAX_MULTI_FRAC:-0.5}
RT_TX_RATIO=${RT_TX_RATIO:-0.5}
RT_TX_COMPARE=${RT_TX_COMPARE:-max}   # compare a bridge with the stronger (max) or weaker (min) touched gene
RT_STRESS=${RT_STRESS:-0}   # 1 = keep stress-elevated read-through (stages D and D2); off since 2026-11-13
if [[ $RT_STRESS == 1 ]]; then RT_MIN_READS_D=10; ELEV_FRAC=0.5; else RT_MIN_READS_D=1000000000; ELEV_FRAC=2; fi
# REF_GTF, GENOME, HOMOLOGY_DB, SAMPLES, BAM_*, tools: pipeline/config.sh
MULTI_DIR=${MULTI_DIR:-support/junction_counts_multi}
HOMOLOGY_CACHE=${HOMOLOGY_CACHE:-$(dirname "${HOMOLOGY_DB:-support/swissprot/x}")}
mkdir -p "$U"
n_tx() { awk '$3=="transcript"' "$1" | wc -l; }
echo "== step 03: $ASM -> $FINAL ($(date)); MIN_SAMPLES=$MIN_SAMPLES MAX_MULTI_FRAC=$MAX_MULTI_FRAC REF_GTF=${REF_GTF:-none}"

# A. union
if [[ ! -s $U/union.combined.gtf ]]; then
  "$GFFCOMPARE" -i "$ASM/gtf_list.scallop2.txt" -o "$U/union" > "$U/gffcompare.log" 2>&1
fi
echo "A  union: $(n_tx "$U/union.combined.gtf") transcripts"

# B. junction rules (step 01 engine; per-sample junction and overhang counts are cached)
if [[ ! -s $U/union.junction.gtf ]]; then
  REF=$U/union.combined.gtf LIST=$ASM/gtf_list.scallop2.txt GTF_SUFFIX=.scallop2.gtf MIN=0 SUPPORT_TABLE=0 \
  FINAL_GTF=$U/union.junction.gtf INTRON_TABLE=$U/union.intron_qc.tsv JOBS=$JOBS \
    pipeline/01_filter_support.sh > "$U/junction_rules.log" 2>&1 || { cat "$U/junction_rules.log"; exit 1; }
fi
grep -E 'introns:|dropped|kept' "$U/junction_rules.log" | sed 's/^ */B  /'

# C. multi-mapped junction counts (cached), paralog filter, sample support with junction rescue
mkdir -p "$MULTI_DIR"
cut -f1 "$CONDITIONS" | xargs -P "$JOBS" -I{} bash -c \
  '[[ -s '"$MULTI_DIR"'/{}.junc_multi.tsv ]] || bash pipeline/lib/junction_multimap_counts.sh '"$BAM_DIR"'/{}'"$BAM_SUFFIX"' '"$MULTI_DIR"'/{}.junc_multi.tsv '"$SAMTOOLS"
for s in $(cut -f1 "$CONDITIONS"); do [[ -s $MULTI_DIR/$s.junc_multi.tsv ]] || { echo "C  multi-mapped counts failed for $s"; exit 1; }; done
if [[ ! -s $U/union.selected.gtf ]]; then
  "$PYTHON" -I pipeline/lib/select_transcripts.py "$U/union.junction.gtf" "$U/union.selected.gtf" \
    support/junction_counts "$MULTI_DIR" "$MIN_SAMPLES" "$MAX_MULTI_FRAC" > "$U/select.log"
fi
sed 's/^/C  /' "$U/select.log"

# D. read-through filter (needs reference gene models)
if [[ -n $REF_GTF ]]; then
  if [[ ! -s $U/union.nort.gtf ]]; then
    ONLY_READTHROUGH=1 FINAL_GTF=$U/union.selected.gtf NORT_GTF=$U/union.nort.gtf RT_PREFIX=$U/union \
    MAKER_GTF=$REF_GTF CONDITIONS=$CONDITIONS RT_MIN_READS=$RT_MIN_READS_D pipeline/01_filter_support.sh > "$U/readthrough.log" 2>&1 \
      || { cat "$U/readthrough.log"; exit 1; }
  fi
  sed 's/^ */D  /' "$U/readthrough.log" | grep -v '^D  $'
  # D2. transcript-level read-through: exons on >= 2 reference genes, however linked
  if [[ ! -s $U/union.nort_tx.gtf ]]; then
    "$PYTHON" -I pipeline/lib/readthrough_tx_filter.py "$U/union.nort.gtf" "$U/union.nort_tx.gtf" "$REF_GTF" \
      "$U/union.tracking" "$ASM/gtf_list.scallop2.txt" "$CONDITIONS" "$RT_TX_RATIO" "$ELEV_FRAC" 3 \
      "$U/union.readthrough_tx.tsv" "$RT_TX_COMPARE" > "$U/readthrough_tx.log"
  fi
  sed 's/^/D2 /' "$U/readthrough_tx.log"
  # D3. protein homology of bridged reference genes (needs reference CDS and a protein database)
  if [[ -n $HOMOLOGY_DB ]]; then
    if [[ ! -s $U/union.nort_hom.gtf ]]; then
      "$PYTHON" -I pipeline/lib/homology_bridges.py "$U/union.nort_tx.gtf" "$U/union.nort_hom.gtf" "$REF_GTF" \
        "$GENOME" "$HOMOLOGY_DB" "$HOMOLOGY_CACHE" "$U/union.bridge_homology.tsv" "$BLASTP" "$THREADS" \
        > "$U/homology.log" 2>&1 || { cat "$U/homology.log"; exit 1; }
    fi
    grep -v '^blastp:' "$U/homology.log" | sed 's/^/D3 /'
    IN_E=$U/union.nort_hom.gtf
  else
    echo "D3 no protein database: homology check skipped"; IN_E=$U/union.nort_tx.gtf
  fi
else
  echo "D  no reference gene models: read-through filter skipped"; IN_E=$U/union.selected.gtf
fi

# E. genes by exon overlap; reference gene ids
"$PYTHON" -I pipeline/lib/assign_genes.py "$IN_E" "$FINAL" "$REF_GTF" SG | sed 's/^/E  /'
echo "== done: $FINAL ($(date))"
