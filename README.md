# Splicing pipeline (water-study *Metrosideros*)

Builds a transcript annotation for splicing analysis from 63 STAR BAMs (4 taxa G/I/L/T; 101 bp paired-end,
fr-firststrand / dUTP), then runs rMATS and GrASE within each taxon (drought vs control, waterlogging vs
control). Designed to be reusable for non-model organisms: the annotation is built from the RNA-seq reads;
a reference gene-model set and Swiss-Prot are used only to catch read-through/fusions and to name genes, and
both are optional. Python helpers are in `lib/`.

## Running it

1. Copy `pipeline/` into a project directory (scripts work relative to `pipeline/..`), with a clone of
   GrASE (`GrASE/`, or set `GRASE_DIR`).
2. Write the sample sheet `support/samples.tsv`: `sample <tab> group <tab> condition`, no header. Within each
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

Tested end to end in a fresh directory (`test_e2e/`: chr10, 2 groups x 3 conditions x 3 samples) on
2026-11-15. This project's first full build used `TAG=scallop2` (outputs `rmats_scallop2*/`, `grase_scallop2/`);
`run_all.sh` writes `rmats/`, `rmats_prefiltered/`, `grase/`.

| Step | Script | What it does | Main output |
|---|---|---|---|
| | `setup.sh` | Checks and one-time setup (BAM indexes, genome index, Swiss-Prot BLAST db) | |
| 02 | `02_assemble_scallop2.sh` | Scallop2 per sample on uniquely mapped reads | `assembly_unique/scallop2/<taxon>/<S>.scallop2.gtf` |
| 03 | `03_union_filter.sh` | Union of the 63 assemblies; junction rules; paralog filter; selection; read-through and homology checks; genes | `assembly_unique/annotation.final.gtf` |
| 08 | `08_run_rmats.sh` | rMATS: read BAMs, detect and count events (`--novelSS --statoff`) | `rmats_scallop2/` |
| 12 | `12_rmats_prefiltered.sh` | rMATS statistics on events passing the junction rules | `rmats_scallop2_prefiltered/` |
| 11 | `11_filter_rmats.py` | rMATS significance filter | `rmats_scallop2_prefiltered/*/*/filtered_jc10_fdr0.05_absdpsi0.1/` |
| 09 | `09_grase_annotation.sh` | GrASE per-gene files, splicing graphs, bipartitions, overhang filter | `grase_scallop2/` |
| 10 | `10_grase_count.sh` | GrASE Stage 1: DEXSeq exonic-part counts | `grase_scallop2/DEXSeq/count_files/` |
| 13 | `13_grase_test.sh` | GrASE Stages 2, 2b, 3: counts onto bipartitions, junction substitution, tests | `grase_scallop2/test_anc12/<taxon>/grase_exsj/` |
| | `run_downstream.sh` | Runs 08 -> 12 -> 11 and 09 -> 10 -> 13 in parallel on one annotation | `rmats[_TAG]/`, `grase[_TAG]/` |
| | `run_all.sh` | setup -> 02 -> 03 -> run_downstream | `run_all.log` |

`01_filter_support.sh` is not run on its own any more: step 03 calls it as the engine for the junction rules
(stage B) and the junction-level read-through filter (stage D). Its per-sample caches are shared by all steps:
`support/junction_counts/` (unique junction reads), `support/junction_counts_anchor/` (overhang),
`support/junction_counts_multi/` (multi-mapped reads). Sample groups: `support/samples.tsv`.

Current build: 2026-11-13 to 2026-11-15. Numbers below are from it.

**History.** The first build (2026-10 to 2026-11-12) assembled with StringTie 2.2.3 per sample (cluster run,
multi-mapped reads included, MAKER guide) plus an Aletsch meta-assembly. It fused neighbouring genes and
near-identical paralogs into single loci (e.g. a receptor-like-protein cluster, a MATE transporter pair, a
read-through across four unrelated genes), and most GrASE events that "replicated" across taxa sat in such
loci. Its scripts are in `obsolete/`; its outputs were deleted on 2026-11-15 (the original per-sample StringTie
assemblies in `taxon*/stringtie/` and the merges `merged.f0.01.gtf`, `merged.f0.05.gtf` are kept as source
data). Lessons that shaped the current design:
- Multi-mapped reads joined paralog copies; assembly now uses unique reads only (MAPQ 255).
- Strict assembler settings (StringTie 3 `-g 0 -c 5 -f 0.15 -j 3`) reduced fusions per sample but every
  merger re-fused them, and mergers (`stringtie --merge`, TACO) kept only 45-57% of well-supported junctions
  against 73% in the per-sample assemblies. The annotation is therefore a filtered union, not a merge.
- Fusions are handled explicitly after the union (read-through and protein-homology checks), not through
  isoform-fraction cut-offs.

---

## 02 `02_assemble_scallop2.sh` - Scallop2 on unique reads

- **Method:** Scallop2 1.1.2 per sample (bioconda build `h503566f_7`; newer builds need AVX2 and crash on this
  CPU): `--library_type first --min_mapping_quality 255 --min_transcript_coverage 3.0`. MAPQ 255 = unique for
  STAR (use 60 for HISAT2); applied while reading, so no filtered BAM copies are needed.
- **Comparison only (`COMPARE_MERGERS=1`, off by default):** `stringtie --merge -g 0 -c 3 -f 0.1` (StringTie
  3.0.0) and `taco_run --isoform-frac 0.1` of the Scallop2 GTFs -> `assembly_unique/merge/`.
- **Output:** 40,384-57,576 transcripts per sample (median 49,546); `assembly_unique/gtf_list.scallop2.txt`.
  ~3 h for 63 samples with 12 workers.

## 03 `03_union_filter.sh` - annotation by evidence-based union

| Stage | What | Result (this build) |
|---|---|---|
| A | `gffcompare -i`: union of the 63 assemblies, identical intron chains collapsed; `num_samples` = samples whose assembly has exactly that chain | 454,391 transcripts |
| B | Junction rules (step 01 engine): every intron canonical (GT-AG/GC-AG/AT-AC), <= 20 kb, >= 3 unique reads in >= 3 samples, each of those with >= 1 read anchored >= 12 bp on both sides (micro-exon introns exempt); single-exon transcripts dropped | 160,585 of 242,071 introns pass (19% non-canonical); 332,472 transcripts |
| C | `lib/select_transcripts.py`: paralog filter (drop transcripts using a junction whose reads are > 50% multi-mapped, `MAX_MULTI_FRAC`), then keep transcripts assembled in >= 5 samples (`MIN_SAMPLES`) plus rarer ones that add a junction no kept transcript has | 2,925 paralog junctions (9,328 transcripts); 101,160 + 31,168 rescued = 132,328 transcripts, no junction lost |
| D | Junction-level read-through (step 01 stages 7-9): junctions joining two reference genes with < 0.2 x the genes' own junction reads | 216 transcripts dropped |
| D2 | `lib/readthrough_tx_filter.py`: transcripts whose exons touch >= 2 reference genes (however linked) dropped if assembled in < 0.5 x the samples of the stronger gene's best own transcript (`RT_TX_RATIO`, `RT_TX_COMPARE=max`) | 1,344 dropped |
| D3 | `lib/homology_bridges.py`: for each pair of reference genes a remaining bridge joins, both reference proteins vs Swiss-Prot (local blastp, E <= 1e-5, top 20). Halves of one protein -> keep; unrelated proteins (fusion) or the same region (paralog chimera) -> drop; no hit -> keep | 2,213 pairs: 567 split genes, 291 different, 259 duplicate, 1,096 unknown; 1,706 transcripts dropped |
| E | `lib/assign_genes.py`: genes = transcripts with overlapping exons on the same strand, but groups touching different reference genes are not joined by exon overlap alone; `ref_gene_id` = overlapping reference genes or `novel` | 6,819 overlap joins blocked |

- **Output:** `assembly_unique/annotation.final.gtf`: **129,062 transcripts, 20,478 genes** (6.3 per gene,
  max 44), 150,813 well-supported junctions (63% of all junctions with >= 3 unique reads in >= 3 samples),
  1,330 genes (6.5%) still joining >= 2 reference genes (mostly gene models the reference split in two).
  Tables in `assembly_unique/union/`: `union.intron_qc.tsv`, `select.log`, `union.readthrough_tx.tsv`,
  `union.bridge_homology.tsv`. ~20 min (BLAST results are cached in `support/swissprot/`).
- **Reference ("MAKER") genes** = all 27,953 genes of `Metrosideros_polymrpha_incana.scaffolded.gtf`
  (21,669 `maker-*`, 3,606 `augustus_masked-*`, 2,678 `snap_masked-*`). Without a reference (`REF_GTF=""`)
  stages D-D3 are skipped and E groups by exon overlap only; without Swiss-Prot (`HOMOLOGY_DB=""`) D3 is
  skipped. Swiss-Prot: UniProt release in `support/swissprot/RELEASE.txt` (2026_04).
- **Stress exemption** (keep read-through elevated under stress) is off (`RT_STRESS=0`).
- **Check (40 sampled bridged pairs, web BLAST vs local):** 39 of 40 classified the same.

## 08 `08_run_rmats.sh` - rMATS detection and counting

- rMATS 4.2.0 per taxon x contrast, groups from `taxon*/metadata/<contrast>/b{1,2}.txt` (stress / control),
  `--readLength 101 --variable-read-length -t paired --libType fr-firststrand`, `--novelSS --statoff`.
- Output `rmats_scallop2/<taxon>/<contrast>/`; ~1 h per contrast with 4 threads.

## 12 `12_rmats_prefiltered.sh` - rMATS statistics on junction-rule events

- `lib/prefilter_rmats_events.py` keeps events whose every junction passes the junction rules (incl.
  overhang; `lib/rmats_junction_rules.py`), `lib/subset_rmats_counts.py` restricts step 08's tables,
  `rmats.py --task stat`; FDR is over kept events. Output `rmats_scallop2_prefiltered/`.

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

- GrASE's `dexseq_count.py` per BAM: `--format bam --order pos --paired yes --stranded reverse -a 10`
  (multi-mapped reads skipped); `count_files/<taxon>_<condition>/`. ~12 h for 63 samples with 12 workers.

## 13 `13_grase_test.sh` - GrASE tests

- Follows `GrASE/scripts/GSE188847/` (stages 30-60), with three groups per taxon and both contrasts per run:
  2. `exoncnt.R`; 2b. `bipartition_sjcnt.R` + `merge_exon_sj_counts.R` + `combine_exoncnt.py` (a side with no
  exclusive exonic part uses its exclusive junctions' reads, from the step 01 BAM counts written as
  `grase_scallop2/SJ/<taxon>_<condition>/<S>.SJ.out.tab`); 3. `exontest.R`: `betabinom_EBmap`,
  `--dpi_gate=both --padj_threshold=0.05 --padj_method=BH --independent_filtering`, per-base pi.
- GrASE's R scripts write their error/warning logs to the directory they are run from (the project root);
  this build's are in `grase_scallop2/root_logs/` (they also hold the first build's entries). Stage 2 failed
  for one gene (`SG010750` internal, "`x` and `y` must share the same src") in all four taxa; every other
  gene was counted. `phi.glmmtmb.errors.*` are convergence warnings from the dispersion fits.
- **Adjustment: plain BH (`PADJ_METHOD=BH`, default since 2026-11-14)**, i.e. FDR over all events of a
  contrast. GrASE's former default `nested_BH` (gene screen, then within-gene BH x m/R) controls a different
  rate and gave 0.6-1x as many calls here (few genes pass its screen, so m/R was 200-5,500).
- **This build** (significant / tested events; ~13 h total):

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
| `bridging_ratio.py`, `readthrough_by_condition.py`, `readthrough_filter.py`, `split_disjoint_genes.py` | 01 stages 7-9 / 03 D | junction-level read-through filter |
| `readthrough_tx_filter.py` | 03 D2 | transcript-level read-through filter |
| `homology_bridges.py` | 03 D3 | Swiss-Prot homology of bridged reference genes |
| `assign_genes.py` | 03 E | gene ids by exon overlap, reference-aware; `ref_gene_id` |
| `rmats_junction_rules.py`, `prefilter_rmats_events.py`, `subset_rmats_counts.py` | 11, 12 | junction rules for rMATS events |
| `grase_anchor_filter.py` | 09 D | overhang filter of GrASE events |

## Environments (`/srv/jupyter-homes/jupyter-eynizade/.conda/envs/`)

`scallop2` (Scallop2 1.1.2, build `h503566f_7`), `stringtie3` (StringTie 3.0.0, build `h29c0135_0`; 3.0.1-3.0.3
builds need AVX2 and crash on this CPU), `stringtie22` (StringTie 2.2.3, gffcompare 0.12.10), `taco` (TACO
0.7.3, comparison only), `blast` (BLAST+ 2.17.0), `rmats420` (rMATS 4.2.0, samtools), `grase` (R 4.4,
GrASE 2.0.0 from `GrASE/Rpkg` with plain BH as default, Python 3.11 + HTSeq; recipe `GrASE/environment.yml`),
`aletsch`, `agat_env`. Packages were installed with `CONDA_PKGS_DIRS=.conda/pkgs_mirahan` (the shared cache is
not writable). Run rMATS through the env's python (`$ENV/bin/python $ENV/rMATS/rmats.py`).

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
