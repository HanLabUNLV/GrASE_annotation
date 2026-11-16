# Splicing annotation and analysis pipeline for non-model organisms

Builds a transcript annotation for splicing analysis from RNA-seq BAMs, then tests differential splicing with
rMATS and GrASE within each sample group (every condition vs a control). The annotation is built from the
reads; a reference gene-model set and Swiss-Prot are used only to catch read-through/fusions and to name genes,
and both are optional. Python helpers are in `lib/`.

Developed on a water-stress study of *Metrosideros* (63 STAR BAMs, 4 taxa x control/drought/waterlogging,
101 bp paired-end, dUTP; sample sheet in `examples/samples.tsv`). Numbers quoted below ("this build") are from
that study's run of 2026-11-13 to 15.

## Running it

1. Copy `pipeline/` into a project directory (scripts work relative to `pipeline/..`), with a clone of
   GrASE (`GrASE/`, or set `GRASE_DIR`).
2. Write the sample sheet `support/samples.tsv`: `sample <tab> group <tab> condition`, no header (example:
   `examples/samples.tsv`, the water-study design: 4 taxa x control/drought/waterlogging). Within each
   group, every condition is compared with `CONTROL` (default `control`). BAMs are
   `<BAM_DIR>/<sample><BAM_SUFFIX>`.
3. Write the project's settings in `<project>/pipeline.config.sh` (plain `VAR=value` lines; outside the
   repository): at least `GENOME`, `REF_GTF` (empty for reference-free), `BAM_DIR`, `BAM_SUFFIX`, `LIBTYPE`,
   `READ_LEN` and `ENVS` (conda envs). `pipeline/config.sh` holds the generic defaults and loads this file
   first; don't edit it per project. Environment variables override both; `PROJECT_CONFIG=<file>` selects
   another project file.
4. `nohup pipeline/run_all.sh > run_all.log 2>&1 &`. It runs `setup.sh` (checks the sample sheet, BAMs and
   indexes, genome index, that every tool runs on this CPU, and downloads/builds Swiss-Prot if needed), then
   steps 02 and 03, then `run_downstream.sh` (rMATS and GrASE in parallel). Every step skips finished work,
   so a stopped run restarts with the same command.

`run_all.sh` writes `rmats/`, `rmats_prefiltered/` and `grase/`; `TAG=<x>` gives `rmats_<x>/` etc. (the
water-study build used `TAG=scallop2`).

| Step | Script | What it does | Main output |
|---|---|---|---|
| | `setup.sh` | Checks and one-time setup (BAM indexes, genome index, Swiss-Prot BLAST db) | |
| 02 | `02_assemble_scallop2.sh` | Scallop2 per sample on uniquely mapped reads | `assembly_unique/scallop2/<group>/<S>.scallop2.gtf` |
| 03 | `03_union_filter.sh` | Union of the assemblies; junction rules; paralog filter; selection; read-through and homology checks; genes | `assembly_unique/annotation.final.gtf` |
| 08 | `08_run_rmats.sh` | rMATS: read BAMs, detect and count events (`--novelSS --statoff`) | `rmats/` |
| 12 | `12_rmats_prefiltered.sh` | rMATS statistics on events passing the junction rules | `rmats_prefiltered/` |
| 11 | `11_filter_rmats.py` | rMATS significance filter | `rmats_prefiltered/<group>/<contrast>/filtered_jc10_fdr0.05_absdpsi0.1/` |
| 09 | `09_grase_annotation.sh` | GrASE per-gene files, splicing graphs, bipartitions, overhang filter | `grase/` |
| 10 | `10_grase_count.sh` | GrASE Stage 1: DEXSeq exonic-part counts | `grase/DEXSeq/count_files/` |
| 13 | `13_grase_test.sh` | GrASE Stages 2, 2b, 3: counts onto bipartitions, junction substitution, tests | `grase/test_anc12/<group>/grase_exsj/` |
| | `run_downstream.sh` | Runs 08 -> 12 -> 11 and 09 -> 10 -> 13 in parallel on one annotation | `rmats[_TAG]/`, `grase[_TAG]/` |
| | `run_all.sh` | setup -> 02 -> 03 -> run_downstream | `run_all.log` |

`01_filter_support.sh` is not run on its own any more: step 03 calls it as the engine for the junction rules
(stage B) and the junction-level read-through filter (stage D). Its per-sample caches are shared by all steps:
`support/junction_counts/` (unique junction reads), `support/junction_counts_anchor/` (overhang),
`support/junction_counts_multi/` (multi-mapped reads). Sample groups: `support/samples.tsv`.

**History.** The water-study's first build (2026-10 to 2026-11-12) assembled with StringTie 2.2.3 per sample (cluster run,
multi-mapped reads included, MAKER guide) plus an Aletsch meta-assembly. It fused neighbouring genes and
near-identical paralogs into single loci (e.g. a receptor-like-protein cluster, a MATE transporter pair, a
read-through across four unrelated genes), and most GrASE events that "replicated" across taxa sat in such
loci. Its scripts are not in this repository, and its outputs were deleted. Lessons that shaped the current
design:
- Multi-mapped reads joined paralog copies; assembly now uses unique reads only (MAPQ 255).
- Strict assembler settings (StringTie 3 `-g 0 -c 5 -f 0.15 -j 3`) reduced fusions per sample but every
  merger re-fused them, and mergers (`stringtie --merge`, TACO) kept only 45-57% of well-supported junctions
  against 73% in the per-sample assemblies. The annotation is therefore a filtered union, not a merge.
- Fusions are handled explicitly after the union (read-through and protein-homology checks), not through
  isoform-fraction cut-offs.

---

## 02 `02_assemble_scallop2.sh` - Scallop2 on unique reads

- **Method:** Scallop2 1.1.2 per sample (bioconda build `h503566f_7`; newer builds need AVX2 and crash on CPUs
  without it): `--library_type $LIBTYPE --min_mapping_quality $MIN_MAPQ --min_transcript_coverage 3.0`.
  MAPQ 255 = unique for STAR (use 60 for HISAT2); applied while reading, so no filtered BAM copies are needed.
- Scallop2 prints its splice graphs even with `--verbose 0` (~3-4 GB per sample); only the last 100 lines of
  each log (`<sample>.log`) are kept.
- **This build:** 40,384-57,576 transcripts per sample (median 49,546); `assembly_unique/gtf_list.scallop2.txt`.
  ~3 h for 63 samples with 12 workers.

## 03 `03_union_filter.sh` - annotation by evidence-based union

**Goal.** Turn the per-sample assemblies of step 02 into one annotation that (1) keeps every splice junction
the reads support well, so every splicing event can be tested, and (2) does not contain transcripts that
the assembler built from artifacts. Downstream, rMATS and GrASE define their events from this annotation,
and GrASE tests events *within a gene*. Two kinds of artifact matter most:
- **false junctions** (misalignments, untrimmed read ends, reads from a near-identical paralog), which
  create events that do not exist;
- **fused genes**: one transcript that runs from one gene into the next (read-through) or joins two
  copies of a duplicated gene makes the two genes one locus, so GrASE tests "path through gene A vs path
  through gene B" as if it were splicing.

**Why a union and not a merge.** Tools that merge assemblies (`stringtie --merge`, TACO) keep only the
dominant isoforms of each locus: here they kept 45-57% of the well-supported junctions, against 73% present
in the per-sample assemblies. Step 03 therefore keeps all assembled transcripts and removes what fails
explicit, data-based rules, stage by stage. Each stage writes its own output in `assembly_unique/union/`, and
a rerun skips stages whose output exists (delete a stage's output to redo it and everything after it).

**Terms.** *Junction* (= intron): a splice from one exon end to the next exon start, identified by
chromosome and first/last intron base. *Well-supported junction*: >= 3 unique reads in >= 3 samples.
*Reference genes*: the gene models in `REF_GTF` (optional). *Bridging transcript*: a transcript whose exons
overlap exons of >= 2 different reference genes.

### Stage A - union of the assemblies (`gffcompare -i`)
All per-sample GTFs are combined; transcripts with an identical intron chain (same junctions in the same
order) become one. Each transcript gets `num_samples`, the number of samples whose own assembly contains
exactly that intron chain. This is the evidence stage C uses: a transcript assembled independently in many
samples is more credible than one assembled once.
-> `union.combined.gtf` (this build: 454,391 transcripts from 63 assemblies)

### Stage B - junction rules (`01_filter_support.sh`)
Every junction of the union is checked against read evidence counted directly from the BAMs (uniquely
mapped, primary alignments; see "Why junction support is counted from the BAMs"). A transcript is kept only
if all its junctions pass:

| Rule | Why |
|---|---|
| Canonical splice motif: GT-AG, GC-AG or AT-AC on the transcript strand (from the genome sequence) | >99% of real plant introns are canonical; non-canonical "introns" are mostly alignment artifacts |
| Intron <= 20 kb (`MAX_INTRON`) | very long introns were weakly supported and chained neighbouring genes together |
| >= 3 unique reads in >= 3 samples (`INTRON_MIN_READS`, `INTRON_MIN_SAMPLES`) | replicated evidence, not one sample's noise |
| Overhang: in each of those samples, >= 1 read with >= 12 aligned bases on both sides (`ANCHOR_MIN`) | reads that reach only a few bases past a "junction" are often untrimmed read ends forced into a splice (e.g. CIGAR `90M65N3M7S`); STAR applies the same 12 bp rule to its own junction report. Introns next to an exon shorter than 12 bp are exempt, because no read can anchor 12 bp in it |

Single-exon transcripts are dropped (they carry no junction and no testable splicing, and often no strand).
Per-sample counts are cached in `support/junction_counts/` and `support/junction_counts_anchor/`.
-> `union.junction.gtf`; `union.intron_qc.tsv` lists every junction with motif, reads, samples, anchored
samples and pass/fail (this build: 160,585 of 242,071 junctions pass, 19% were non-canonical; 332,472
transcripts left)

### Stage C - paralog filter and transcript selection (`lib/select_transcripts.py`)
1. **Paralog filter.** In families of near-identical gene copies, reads from one copy also align to the
   others; even reads that the aligner calls unique may be placed on the wrong copy, and the assembler then
   builds chimeric transcripts across copies. For every junction the reads from multi-mapped alignments are
   counted (`support/junction_counts_multi/`, `lib/junction_multimap_counts.sh`); a junction where more than
   half of the reads are multi-mapped (`MAX_MULTI_FRAC=0.5`) marks such a region, and transcripts using it
   are dropped.
2. **Selection by sample support.** Assemblers combine known junctions into many rare intron chains
   (about half of the union's transcripts were assembled in one sample only); most add no new junction, but each one adds paths that GrASE must
   enumerate. A transcript is kept if its intron chain was assembled in >= 5 samples (`MIN_SAMPLES`).
3. **Junction rescue.** Rarer transcripts are then added back, most supported first, whenever they carry a
   junction that no kept transcript has. So selection never loses a junction that passed stage B; it only
   removes redundant combinations.
-> `union.selected.gtf`, `select.log` (this build: 2,925 paralog junctions removed 9,328 transcripts;
101,160 transcripts with >= 5 samples + 31,168 rescued = 132,328, with all 153,489 junctions retained)

### Stage D - read-through, junction level (`01_filter_support.sh`, `lib/bridging_ratio.py`)
Read-through transcription continues past the end of one gene into the next and splices an exon of gene A
to an exon of gene B. It is usually a minor product. For each junction whose donor lies in one reference
gene and acceptor in another, the ratio of its reads to the typical (median) junction reads inside the two
genes is computed. A ratio < 0.2 means a minor bridge next to two active genes: the junction is dropped,
with the transcripts using it. A ratio around 1 or higher means the bridge is as strong as the genes'
own splicing, which usually means the reference split one real gene into two models; these are kept.
-> `union.nort.gtf`, `union.bridging_ratio.tsv`, `union.readthrough_by_condition.tsv` (this build: 216
transcripts dropped)

### Stage D2 - read-through, transcript level (`lib/readthrough_tx_filter.py`)
Stage D only sees junctions with both ends inside reference genes. Read-through also appears as an exon
running across the gap between two genes, or as a junction into an unannotated exon just outside a gene
model. Stage D2 therefore looks at whole transcripts: every bridging transcript is compared with the
transcripts that each touched gene has on its own (support = `num_samples`):
- **dropped** if it was assembled in fewer than 0.5 x the samples of the stronger gene's best own
  transcript (`RT_TX_RATIO=0.5`, `RT_TX_COMPARE=max`): a minor product next to established genes;
- **kept** if it is as common as the genes' own transcripts (the joined form is what is transcribed), or
  if a touched gene has no transcript of its own (the gene exists only as part of the joined form, as when
  the reference split one gene).
-> `union.nort_tx.gtf`, `union.readthrough_tx.tsv` with the decision per bridging transcript (this build:
1,344 dropped)

### Stage D3 - protein homology of bridged genes (`lib/homology_bridges.py`)
After D2 the remaining bridges are frequent, and frequency cannot distinguish a gene that the reference
split in two from two unrelated genes that are always co-transcribed. Protein homology can. For every pair
of reference genes that a remaining bridging transcript joins, both reference proteins (first coding
model of each gene, translated from the genome) are searched against UniProt Swiss-Prot with local
`blastp` (E <= 1e-5, top 20 hits):

| Result | Interpretation | Action |
|---|---|---|
| both hit the same protein, in different regions (e.g. residues 1-200 and 220-600) | halves of one protein: the reference split one gene | keep the joined transcripts |
| both have hits, to different proteins | two different genes: fusion / read-through | drop the bridging transcripts |
| both hit the same protein in the same region | two copies of a duplicated gene: paralog chimera | drop |
| no hit for one or both | no evidence either way | keep |

BLAST results are cached in `support/swissprot/ref_proteins.blastp.tsv`, so reruns only search new genes.
-> `union.nort_hom.gtf`, `union.bridge_homology.tsv` with the class of every gene pair (this build: 2,213
pairs: 567 split genes, 291 different proteins, 259 duplicates, 1,096 without hits; 1,706 transcripts
dropped)

### Stage E - genes (`lib/assign_genes.py`)
Transcripts are grouped into genes: transcripts on the same strand whose exons overlap belong to one gene
(overlapping *spans* are not enough, so a gene inside another gene's intron, or two neighbours, stay
separate). With a reference, two groups that touch *different* reference genes are not joined by exon
overlap alone (neighbouring genes whose UTRs extend into each other); only a transcript that itself touches
both genes (a bridge kept by D2/D3) joins them. Genes are numbered `SG000001`... by position, and every
transcript gets `ref_gene_id`: the reference genes its exons overlap, or `novel`.
-> `annotation.final.gtf` (this build: 6,819 exon-overlap joins between different reference genes blocked)

### Output and options
- **`assembly_unique/annotation.final.gtf`** (this build): **129,062 transcripts in 20,478 genes** (6.3 per
  gene, max 44); 150,813 well-supported junctions (63% of all well-supported junctions in the data; the rest
  are absent from every assembly or fail a rule above: motif, length, overhang, paralog or read-through); 1,330 genes (6.5%) still join >= 2 reference genes,
  mostly genes the reference split in two. Run time ~20 min.
- **Without a reference** (`REF_GTF=""`): stages D-D3 are skipped and E groups by exon overlap only.
  **Without Swiss-Prot** (`HOMOLOGY_DB=""`): D3 is skipped.
- **Settings:** `MIN_SAMPLES` (5), `MAX_MULTI_FRAC` (0.5), `RT_TX_RATIO` (0.5), `RT_TX_COMPARE` (max),
  `RT_STRESS` (0; 1 keeps read-through that is elevated under a treatment in some group, in D and D2).
- **Water-study reference:** all 27,953 MAKER genes (21,669 `maker-*`, 3,606 `augustus_masked-*`, 2,678
  `snap_masked-*`); Swiss-Prot release 2026_04. A web-BLAST check of 40 bridged gene pairs agreed with the
  local classification for 39.

## 08 `08_run_rmats.sh` - rMATS detection and counting

- rMATS 4.2.0 per group and contrast (`<condition>_vs_<CONTROL>`), b1 = condition, b2 = control, both built
  from `support/samples.tsv`; `--readLength $READ_LEN --variable-read-length -t paired|single` and `--libType`
  from `LIBTYPE`; `--novelSS --statoff`.
- Output `rmats/<group>/<contrast>/`; ~1 h per contrast with 4 threads (this build).

## 12 `12_rmats_prefiltered.sh` - rMATS statistics on junction-rule events

- `lib/prefilter_rmats_events.py` keeps events whose every junction passes the junction rules (incl.
  overhang; `lib/rmats_junction_rules.py`), `lib/subset_rmats_counts.py` restricts step 08's tables,
  `rmats.py --task stat`; FDR is over kept events. Output `rmats_prefiltered/`.

## 11 `11_filter_rmats.py` - rMATS significance filter

- The study's filter (tag `jc10_fdr0.05_absdpsi0.1`): mean (IJC + SJC) per sample >= 10 in one group,
  FDR < 0.05, |IncLevelDifference| >= 0.1; adds `ref_gene_id`. Optional `--junc_dir` re-applies the junction
  rules (`_jqc_anc12` output); not needed after step 12.
- **This build:** 16,344 significant events (A3SS 2,704, A5SS 1,919, MXE 197, RI 9,345, SE 2,179).
- **Caution:** rMATS p-values are inflated here (11-17% of events at p < 0.05 vs ~5% for GrASE and for a
  Welch t-test on per-sample PSI); in the first build only 32-41% of rMATS calls reached even unadjusted
  p < 0.05 on that t-test, because replicate PSI SD (0.12-0.15) exceeds the 0.1 dPSI cut-off. Use GrASE
  calls, or rMATS calls confirmed by a replicate-aware test.

## 09 `09_grase_annotation.sh` - GrASE annotation

- A. `GrASE/scripts/createfiles.sh` (per-gene GTFs and DEXSeq GFFs, concatenated GFF, splicing graphs);
  B. `bubble_path_split.R --split=bipartition`; C. `filterTerminal.R` (internal vs terminal);
  D. `lib/grase_anchor_filter.py` (overhang rule; a no-op for annotations built by step 03).
- Use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` and `JOBS` = cores for stage B.
- **This build:** 16,231 graphs; 12,992 internal and 15,225 terminal bipartition files, 268,400 events (14
  removed by D). ~12 h (stage B ~10 h with 12 workers).

## 10 `10_grase_count.sh` - GrASE counting

- GrASE's `dexseq_count.py` per BAM: `--format bam --order pos -a 10`, `--paired` from `PAIRED` and
  `--stranded` from `LIBTYPE` (reverse for dUTP) (multi-mapped reads skipped); `count_files/<group>_<condition>/`.
  ~12 h for 63 samples with 12 workers (this build).

## 13 `13_grase_test.sh` - GrASE tests

- Follows `GrASE/scripts/GSE188847/` (stages 30-60), with all conditions of a group and all its contrasts in one
  run, so the dispersion is estimated over the whole group:
  2. `exoncnt.R`; 2b. `bipartition_sjcnt.R` + `merge_exon_sj_counts.R` + `combine_exoncnt.py` (a side with no
  exclusive exonic part uses its exclusive junctions' reads, from the step 01 BAM counts written as
  `grase/SJ/<group>_<condition>/<S>.SJ.out.tab`); 3. `exontest.R`: `betabinom_EBmap`,
  `--dpi_gate=both --padj_threshold=0.05 --padj_method=BH --independent_filtering`, per-base pi.
- GrASE's R scripts write their error/warning logs to the directory they are run from; step 13 runs them from
  `grase/test_anc12/<group>/` (`exoncnt.bipartition.errors.log`; `phi.glmmtmb.errors.*` are convergence
  warnings from the dispersion fits). In this build stage 2 failed for one gene (`SG010750` internal, "`x`
  and `y` must share the same src") in all four taxa; every other gene was counted.
- **This build** (water study; significant / tested events; ~13 h total):

| Taxon | Internal: drought | Internal: waterlogging | Terminal: drought | Terminal: waterlogging |
|---|---|---|---|---|
| G | 14 / 124,273 | 14 / 123,330 | 34 / 250,305 | 31 / 243,363 |
| I | 0 / 123,467 | 1 / 123,686 | 3 / 245,105 | 22 / 244,146 |
| L | 20 / 124,265 | 14 / 124,551 | 30 / 245,704 | 26 / 249,671 |
| T | 15 / 124,407 | 6 / 124,107 | 58 / 254,184 | 53 / 251,740 |

- Signal is weak: p-values are near-uniform and replicate variation exceeds most effects (4-6 samples per
  group). Replication of significant events (p < 0.05, same direction) in an independent taxon under the same
  stress: 12% (chance 2.5%); different taxon and stress: 5.2%. Genes significant in >= 3 of the 8 contrasts:
  `SG017962` (DTX42-like MATE transporter; G and L, both stresses, path proportion -0.07 to -0.15),
  `SG011327` (DTX42/43-like MATE; G, I, T), `SG016312` (2-oxoglutarate-dependent dioxygenase; L, T).

## lib/

| Script | Used by | Purpose |
|---|---|---|
| `junction_anchor_counts.sh` | 01 / 03 B | per-sample junction reads with >= 12 bp aligned on both sides |
| `junction_multimap_counts.sh` | 03 C | per-sample junction reads from multi-mapped primary alignments |
| `select_transcripts.py` | 03 C | paralog filter; sample-support selection with junction rescue |
| `bridging_ratio.py`, `readthrough_by_condition.py`, `readthrough_filter.py`, `split_disjoint_genes.py` | 01 stages 4-6 / 03 D | junction-level read-through filter |
| `readthrough_tx_filter.py` | 03 D2 | transcript-level read-through filter |
| `homology_bridges.py` | 03 D3 | Swiss-Prot homology of bridged reference genes |
| `assign_genes.py` | 03 E | gene ids by exon overlap, reference-aware; `ref_gene_id` |
| `rmats_junction_rules.py`, `prefilter_rmats_events.py`, `subset_rmats_counts.py` | 11, 12 | junction rules for rMATS events |
| `grase_anchor_filter.py` | 09 D | overhang filter of GrASE events |

## Environments

Conda envs under `ENVS` (set in `pipeline.config.sh`; tool paths are in `config.sh`). Used for the water study:
`scallop2` (Scallop2 1.1.2, build `h503566f_7`), `stringtie22` (gffcompare 0.12.10), `blast` (BLAST+ 2.17.0),
`rmats420` (rMATS 4.2.0, samtools 1.21), `grase` (R 4.4, GrASE 2.0.0 from `GrASE/Rpkg`, Python 3.11 + HTSeq,
pandas, scipy; recipe `GrASE/environment.yml`). `setup.sh` checks that each tool runs on the current CPU. Run rMATS through the env's
python (`$ENV/bin/python $ENV/rMATS/rmats.py`).

## Why junction support is counted from the BAMs, not from STAR's SJ.out.tab

1. **The downstream tools read the BAM**; STAR's `outSJfilter*` settings only decide what is written to
   `SJ.out.tab` (with `outFilterType Normal` the reads stay in the BAM).
2. **SJ.out.tab under-reports junctions that are real in the BAM** (first build: 1,218 introns absent from
   every sample's `SJ.out.tab`; all 12 spot-checked had spliced reads, up to 451).
3. **The quality criteria are applied explicitly instead**, uniformly across samples:

| | STAR `SJ.out.tab` (these runs' settings) | This pipeline |
|---|---|---|
| Splice motif | non-canonical allowed with >= 3 unique reads | canonical only (GT-AG/GC-AG/AT-AC) |
| Read support | >= 1 unique read (canonical), per sample | >= 3 unique reads in >= 3 samples |
| Multi-mapped reads | counted separately | excluded (MAPQ 255); junctions > 50% multi-mapped dropped |
| Intron length | longer introns need more reads (up to ~590 kb) | <= 20 kb |
| Overhang (aligned bases each side) | >= 12 bp for unannotated canonical | each supporting sample has >= 1 read with >= 12 bp each side (micro-exon introns exempt) |

The overhang rule catches untrimmed read ends forced into a splice (e.g. reads `90M65N3M7S` giving a 65 bp
"intron" with 12,305 reads and none anchored).

## Notes

- Do not edit a script while it is running: bash reads scripts incrementally.
- When stopping jobs, match processes by pid or exact argument; `pkill -f <pattern>` also matches the shell
  that runs it, and so does a watcher loop whose own command line contains the pattern.
- Kill an `xargs` parent before its workers, or it starts the next job for each one killed.
