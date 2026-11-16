# Project settings, read by every step (source'd after the step cd's to the project directory, the
# parent of pipeline/). Each value can be overridden from the environment: VAR=... pipeline/<step>
# Paths are relative to the project directory unless absolute.
# A different config file: CONFIG=/path/to/config.sh pipeline/<step>

# --- data ---------------------------------------------------------------------------
: "${SAMPLES:=support/samples.tsv}"          # sample <tab> group <tab> condition (no header)
: "${CONTROL:=control}"                      # condition every other condition is compared against, per group
: "${BAM_DIR:=BAMs/star_alignmentuntrimmed_results}"   # <BAM_DIR>/<sample><BAM_SUFFIX> (+ .bai)
: "${BAM_SUFFIX:=.Aligned.sortedByCoord.out.bam}"
: "${GENOME:=Metrosideros_polymrpha_incana.scaffolded.fa}"
: "${REF_GTF:=Metrosideros_polymrpha_incana.scaffolded.gtf}"   # reference gene models; "" = none (reference-free)

# --- library --------------------------------------------------------------------------
: "${LIBTYPE:=first}"        # first = fr-firststrand (dUTP), second = fr-secondstrand, unstranded
: "${PAIRED:=yes}"           # yes / no
: "${READ_LEN:=101}"         # read (mate) length
: "${MIN_MAPQ:=255}"         # unique alignments: 255 for STAR, 60 for HISAT2

# --- protein homology check (step 03 D3) ------------------------------------------------
: "${HOMOLOGY_DB=support/swissprot/sprot}"   # BLAST db of UniProt Swiss-Prot, built by setup.sh; "" = skip

# --- software ---------------------------------------------------------------------------
: "${ENVS:=/srv/jupyter-homes/jupyter-eynizade/.conda/envs}"
: "${GRASE_DIR:=GrASE}"      # clone of github.com/HanLabUNLV/GrASE (scripts are run from it)
: "${SCALLOP2:=$ENVS/scallop2/bin/scallop2}"
: "${GFFCOMPARE:=$ENVS/stringtie22/bin/gffcompare}"
: "${SAMTOOLS:=$ENVS/rmats420/bin/samtools}"
: "${BLASTP:=$ENVS/blast/bin/blastp}"
: "${MAKEBLASTDB:=$ENVS/blast/bin/makeblastdb}"
: "${RMATS_ENV:=$ENVS/rmats420}"             # python + rMATS/rmats.py
: "${GRASE_ENV:=$ENVS/grase}"                # R + grase package, python + HTSeq
: "${PYTHON:=$GRASE_ENV/bin/python}"         # pandas, scipy
: "${STRINGTIE:=$ENVS/stringtie3/bin/stringtie}"   # only for COMPARE_MERGERS=1 in step 02
: "${TACO_RUN:=$ENVS/taco/bin/taco_run}"           # only for COMPARE_MERGERS=1 in step 02

# --- resources --------------------------------------------------------------------------
: "${THREADS:=16}"
: "${JOBS:=12}"

# --- derived (do not edit) ------------------------------------------------------------
case $LIBTYPE in
  first)      RMATS_LIBTYPE=fr-firststrand; DEXSEQ_STRANDED=reverse ;;
  second)     RMATS_LIBTYPE=fr-secondstrand; DEXSEQ_STRANDED=yes ;;
  unstranded) RMATS_LIBTYPE=fr-unstranded;  DEXSEQ_STRANDED=no ;;
  *) echo "config: LIBTYPE must be first, second or unstranded (got '$LIBTYPE')" >&2; exit 1 ;;
esac
[[ $PAIRED == yes ]] && RMATS_READTYPE=paired || RMATS_READTYPE=single
CONDITIONS=$SAMPLES          # name used by the older step code
MAKER_GTF=$REF_GTF
# groups (column 2) and, per group, the conditions other than CONTROL
pipeline_groups() { cut -f2 "$SAMPLES" | sort -u; }
pipeline_treatments() { awk -F'\t' -v g="$1" -v c="$CONTROL" '$2==g && $3!=c {print $3}' "$SAMPLES" | sort -u; }
pipeline_conditions() { awk -F'\t' -v g="$1" '$2==g {print $3}' "$SAMPLES" | sort -u; }
pipeline_samples() { awk -F'\t' -v g="$1" -v c="$2" '$2==g && $3==c {print $1}' "$SAMPLES"; }
export SAMPLES CONTROL BAM_DIR BAM_SUFFIX GENOME REF_GTF MAKER_GTF CONDITIONS LIBTYPE PAIRED READ_LEN MIN_MAPQ \
  HOMOLOGY_DB ENVS GRASE_DIR SCALLOP2 GFFCOMPARE SAMTOOLS BLASTP MAKEBLASTDB RMATS_ENV GRASE_ENV PYTHON \
  STRINGTIE TACO_RUN THREADS JOBS RMATS_LIBTYPE DEXSEQ_STRANDED RMATS_READTYPE
