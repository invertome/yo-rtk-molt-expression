# Receptor tyrosine kinase expression across the crustacean molt cycle

Code and data for the analysis of receptor tyrosine kinase (RTK) expression in the accompanying
paper. The repository holds the sample sheets, the RTK gene panels, the configuration of each
dataset, the analysis code, and the commands that regenerate every table and figure of the paper
(Figures 1 to 16, Figures S1 to S5, Tables S1 and S2) from the Salmon quantifications.

**The question.** RTKs of the insulin, epidermal growth factor, fibroblast growth factor and
PDGF/VEGF receptor families may carry molt-related signals to the tissues that drive the
crustacean molt cycle. For each RTK transcript in each tissue or experiment, the analysis asks
whether its expression changes across the molt cycle and, if it does, between which stages.

**The datasets.**

- *Gecarcinus lateralis* (the blackback land crab): Y-organ (YO, the molting gland), 115
  libraries in two experiments analyzed separately. In the multiple leg autotomy (MLA)
  experiment, 54 animals were sampled at six stages: intermolt (IM), early premolt (EP),
  mid-premolt (MP), late premolt (LP), postmolt (PM), and blocked animals (BL) that did not enter
  premolt. Stages were assigned from the regeneration index (R) and the hemolymph ecdysteroid
  titer. In the eyestalk ablation (ESA) experiment, 59 males were sampled at days 0, 1, 3, 5, 7, 14
  and 19 after ablation; the intact animals are Day 0 controls. The panel holds 17 RTK
  transcripts, one per gene.
- *Carcinus maenas* (the green crab): 25 YO and 25 central nervous system (CNS) libraries from
  paired tissues, at four molt stages: intermolt (IM), mid-premolt (MP), late premolt (LP) and
  postmolt (PM). The panel holds 29 RTK transcripts of 23 genes; the transcripts of a gene are
  summed and the gene is analyzed.

## Reproduce

Work from the repository root, with the quantifications in place (see "Inputs"). The R side runs
once per dataset:

    cp pipeline/site.env.example site.env      # optional: names the conda environment that holds R and pins the OpenBLAS kernel
    bash pipeline/run_pipeline.sh configs/glat.yaml
    bash pipeline/run_pipeline.sh configs/cmaenas.yaml

Each run ends with `VALIDATION PASSED` and `PIPELINE DONE` and writes its tables to
`results/glat/` or `results/cmaenas/`. A run takes about ten minutes on four cores. On a Slurm
cluster, `pipeline/run_pipeline.sbatch` wraps the same script (see `pipeline/README.md`,
"Running"). The tables and figures are then built from the results:

    export BASE_LOCAL=.
    python3 pipeline/table_panel.py configs/glat.yaml    results/glat    tables
    python3 pipeline/table_panel.py configs/cmaenas.yaml results/cmaenas tables
    python3 pipeline/fig_main.py    configs/glat.yaml    results/glat    figures/main
    python3 pipeline/fig_main.py    configs/cmaenas.yaml results/cmaenas figures/main
    python3 pipeline/fig_supp.py    configs/glat.yaml    results/glat    figures/supplementary
    python3 pipeline/fig_supp.py    configs/cmaenas.yaml results/cmaenas figures/supplementary

The expected molt-regulated counts are 4 of the 17 RTK transcripts in the MLA experiment, 6 of 17
in the ESA experiment, 4 of the 23 RTK genes in the *C. maenas* YO, and 0 of 23 in the
*C. maenas* CNS.

## How the analysis works

`pipeline/README.md` describes the general method step by step. The choices for these datasets
are set in `configs/glat.yaml` and `configs/cmaenas.yaml`, whose comments explain each value.

1. **Four analysis units**, each fitted on its own: MLA and ESA for *G. lateralis*, YO and CNS
   for *C. maenas*. The configuration gives the number of libraries expected at every stage, and
   the run checks it against the sample sheet before any model is fitted.
2. **Transcriptome-wide DESeq2 model** with one mean per stage, fitted to every quantified
   transcript (for *C. maenas*, to every gene after summing its transcripts). Panel transcripts
   stay in the fit even when barely expressed.
3. **Test for a molt-stage effect.** A likelihood-ratio test compares that model with a reduced
   model in which the stages share one mean. In the MLA experiment the blocked animals keep their
   own mean in the reduced model, so they stay in the fit while the test asks only whether
   expression varies across the five stages of the cycle. DESeq2's Cook's-distance outlier
   exclusion is switched off in every test, so a single high-count animal cannot remove a
   transcript from the analysis.
4. **Which transcripts are molt-regulated.** The likelihood-ratio p-values of the panel
   transcripts are adjusted across the panel only (Benjamini-Hochberg over 17 or 23 tests, the
   targeted false discovery rate); a transcript is molt-regulated when that value is below 0.05.
   The transcriptome-wide adjusted p-value is reported alongside but does not decide the call.
5. **Where in the cycle expression changes.** For each molt-regulated transcript, every pair of
   stages is compared with a Wald test, the p-values are adjusted within the transcript with the
   Holm method, and the result is drawn as compact letters over the bars: stages that share a
   letter do not differ. A transcript that is not molt-regulated has the same letter at every
   stage in the letters table and no letters on the figures.
6. **Supplementary comparisons.** For every step through the cycle (IM to EP, EP to MP, and so on,
   including PM back to IM), for each day after ablation against the previous one, and for the
   blocked animals against every MLA stage, the shrunken log2 fold changes of every transcript
   are tabulated and drawn as volcano and MA plots with the panel transcripts marked (Figures S1
   to S5). A principal component analysis of each unit is written alongside.
7. **Validation.** After every run the output tables are checked for internal consistency (one
   row per transcript and stage, sample counts equal to the sheet, calls equal to the adjusted
   p-values, every missing test result explained); a failed check stops the run.

## Layout

    configs/                 one configuration file per dataset (glat.yaml, cmaenas.yaml)
    data/glat/               samples.tsv (115 libraries) and glat_rtk_genes.tsv (the 17-transcript panel)
    data/cmaenas/            cmaenas_samplesheet.tsv (50 libraries) and cmaenas_panel_reps.tsv (29 transcripts, 23 genes)
    pipeline/                the analysis code (dge_panel_pipeline 1.0.4); see pipeline/README.md
    species_photos/          species silhouettes drawn on the main figures, with attribution
    results/glat/            DESeq2 output tables for MLA and ESA
    results/cmaenas/         DESeq2 output tables for YO and CNS
    figures/main/            Figures 1 to 16
    figures/supplementary/   Figures S1 to S5
    tables/                  Tables S1 and S2

## Inputs

**Quantifications.** The pipeline starts from one Salmon `quant.sf` file per library. These files
are not in the repository because of their size. Raw reads are in the NCBI Sequence Read Archive
under BioProject [SRA BioProject accession], and the quantifications are on Zenodo
([Zenodo DOI]). Place them as

    data/glat/salmon_output/<filename>/quant.sf     one directory per row of data/glat/samples.tsv (column filename)
    data/cmaenas/salmon_output/<run>/quant.sf       one directory per row of data/cmaenas/cmaenas_samplesheet.tsv (column run)

The *G. lateralis* libraries were quantified with Salmon 1.10.2 and the *C. maenas* libraries
with Salmon 1.10.3, both with automatic library-type detection.

**Sample sheets.** `data/glat/samples.tsv` has one row per *G. lateralis* library:

| Column | Holds |
| :--- | :--- |
| `sampleid` | the library identifier |
| `filename` | the name of the library's quantification directory |
| `experiment` | `MLA` or `ESA` |
| `stage` | MLA: `IM`, `EP`, `MP`, `LP`, `PM`, `BL`. ESA: `D0`, `D1`, `D3`, `D5`, `D7`, `D14`, `D19` (days after ablation; the intact animals are listed as `D0`), or `female` for the two female libraries, which stay in the sheet and are excluded by the configuration |

`data/cmaenas/cmaenas_samplesheet.tsv` has one row per *C. maenas* library:

| Column | Holds |
| :--- | :--- |
| `run` | the SRA run accession, which also names the quantification directory |
| `crab` | the animal identifier |
| `aiken_stage` | the molt stage in the finer staging of the original study (C3, C4, D1, D3, D4, A, B) |
| `tissue` | `YO` or `CNS` |
| `mykles_stage` | the four stages used in the analysis: `IM` (C3 and C4), `MP` (D1), `LP` (D3 and D4), `PM` (A and B) |

**Gene panels.** `data/glat/glat_rtk_genes.tsv` lists the 17 *G. lateralis* RTK transcripts with
their gene name, transcript identifier, family (`InsR`, `EGFR`, `FGFR`, `PVR`), and whether the
transcript is partial. `data/cmaenas/cmaenas_panel_reps.tsv` lists the 29 *C. maenas* RTK
transcripts of 23 genes by gene, transcript identifier and family; each of the six genes with
several transcripts is quantified as the sum of its transcripts.

## Outputs

**Result tables** in `results/glat/` (units `MLA` and `ESA`; the unit column is `experiment`) and
`results/cmaenas/` (units `YO` and `CNS`; the unit column is `tissue`). `<unit>` stands for one of
those names.

| File | Columns |
| :--- | :--- |
| `<unit>_rtk_letters.tsv` | one row per RTK transcript and stage: `gene`; the unit column; `stage`; `letter` (compact letter at that stage); `lrt_p` (raw p-value of the likelihood-ratio test); `q_targeted` (adjusted across the panel; the call is made on it); `q_transcriptomewide` (DESeq2's adjusted p-value across all transcripts, the transcriptome-wide false discovery rate); `molt_regulated` (`TRUE` when `q_targeted` is below 0.05). The test values repeat on every stage row of a transcript. |
| `<unit>_rtk_tpm_by_stage.tsv` | `gene`; `id` (the transcript, or the summed transcripts joined by `;`); the unit column; `stage`; `n` (libraries at that stage); `mean_tpm`; `sem_tpm` (standard error of the mean) |
| `glat_rtk_persample_tpm.tsv`, `cmaenas_rtk_persample_tpm.tsv` | the transcripts per million (TPM) of every RTK transcript in every library: `gene`; `id`; the unit column; `sampleid`; `stage`; `tpm` |
| `cmaenas_rtk_panel.tsv` | *C. maenas* only: `gene_id`; `rep_ids` (the transcripts summed for that gene, joined by `;`); `n_reps` |
| `<unit>_deseq2_LRT_all.tsv` | the likelihood-ratio test of every transcript: `feature`; `baseMean` (mean normalized count); `log2FoldChange` and `lfcSE` (the last stage against the first, a by-product of the fit that the test does not depend on); `stat` (likelihood-ratio statistic); `pvalue`; `padj` (transcriptome-wide Benjamini-Hochberg) |
| `<unit>_dge_<A>_vs_<B>.tsv` | one file per pairwise comparison, `A` relative to `B` (a positive fold change means higher expression in `A`), every transcript: `feature`; `baseMean`; `log2FoldChange` (shrunken); `lfcSE`; `stat` (Wald statistic); `pvalue`; `padj`; `gene` and `family` (RTK transcripts only); `is_rtk` (`TRUE` for an RTK transcript); `rtk_sig` (`TRUE` for an RTK transcript with `padj` below 0.05 in this comparison; these are the starred transcripts in Figures S1 to S5); the unit column; `contrast` |
| `<unit>_pca.tsv` | `sample`; `PC1`; `PC2`; `stage`: the first two principal components of the 500 most variable transcripts after variance-stabilizing transformation |
| `<unit>_pca_var.tsv` | `PC`; `pct`: the percentage of variance carried by each component |

Each run also writes `<unit>_dds.rds`, the fitted DESeq2 object, which the repository does not
track. The pairwise comparisons are: MLA `EP_vs_IM`, `MP_vs_EP`, `LP_vs_MP`, `PM_vs_LP`,
`IM_vs_PM` and `BL_vs_<stage>` for each of the five stages; ESA `D1_vs_D0` through `D19_vs_D14`;
YO and CNS `MP_vs_IM`, `LP_vs_MP`, `PM_vs_LP`, `IM_vs_PM`.

**Figures.** Each figure is written as PNG (300 dpi) and PDF (vector, fonts embedded, silhouette
embedded at 600 dpi).

| File | Paper |
| :--- | :--- |
| `figures/main/Fig01_Gl_InsR_MLA` to `Fig08_Gl_PVR_ESA` | Figures 1 to 8: *G. lateralis*, one figure per family (InsR, EGFR, FGFR, PVR) and experiment (MLA, then ESA) |
| `figures/main/Fig09_Cm_InsR_YO` to `Fig16_Cm_PVR_CNS` | Figures 9 to 16: *C. maenas*, one figure per family and tissue (YO, then CNS) |
| `figures/supplementary/FigS1_Gl_MLA_transitions` | Figure S1: MLA, consecutive molt stages |
| `figures/supplementary/FigS2_Gl_MLA_blocked` | Figure S2: MLA, blocked animals against each molt stage |
| `figures/supplementary/FigS3_Gl_ESA_days` | Figure S3: ESA, consecutive days after ablation |
| `figures/supplementary/FigS4_Cm_YO_transitions` | Figure S4: *C. maenas* YO, consecutive molt stages |
| `figures/supplementary/FigS5_Cm_CNS_transitions` | Figure S5: *C. maenas* CNS, consecutive molt stages |

In Figures 1 to 16, each panel shows one RTK transcript: pale bars give the mean TPM at each stage
with its standard error, one dot per library, and compact letters over the molt-regulated
transcripts; the ESA figures draw a line over days instead of bars. In Figures S1 to S5, each row
is one pairwise comparison, with a volcano plot on the left (log2 fold change against -log10
adjusted p-value) and an MA plot on the right (log2 fold change against mean expression); every
transcript is a gray point, RTK transcripts are red, and those significant in that comparison are
starred and named. `figures/supplementary/glat_stars.json` and `cmaenas_stars.json` list, for each
supplementary figure, the (comparison, gene) pairs that received a star, which can be checked
against the `rtk_sig` column of the result tables.

**Tables.** `tables/TableS1_Gl` (Table S1, *G. lateralis*) and `tables/TableS2_Cm` (Table S2,
*C. maenas*), each as `.tsv` and `.md`, list every RTK transcript in every unit with `gene`,
`family`, the unit column, `LRT_p`, `q_targeted`, `q_transcriptomewide` and `molt_regulated`. The
Markdown version writes the numbers in scientific notation.

## Software

R 4.5.3 with DESeq2 1.50.2, tximport 1.38.2, readr 2.2.0, data.table 1.18.6.1, multcompView
0.1-12, yaml 2.3.12 (`pipeline/env_r.yaml`); the exact resolved environment is listed in
`env_explicit.txt`. The table and figure scripts need Python 3 with pandas, matplotlib, numpy,
PyYAML and Pillow.

The tables in `results/` are computed with the OpenBLAS kernel pinned
(`OPENBLAS_CORETYPE=HASWELL`, set in `site.env`). Keep that setting to reproduce them to the last
digit: OpenBLAS builds that select a CPU-specific kernel at run time round differently on
different processor generations, which moves the last decimal places of the transcriptome-wide
coefficient estimates. The figures use the Liberation fonts; without them matplotlib falls back to
another typeface and the data drawn stay the same.

## Checking a rerun

To check that a rerun reproduces the tables in `results/`, write the rerun to a separate
directory and compare the two, table by table:

    PIPELINE_OUT=tmp/rerun/glat    bash pipeline/run_pipeline.sh configs/glat.yaml
    PIPELINE_OUT=tmp/rerun/cmaenas bash pipeline/run_pipeline.sh configs/cmaenas.yaml
    python3 pipeline/compare_outputs.py results/glat    tmp/rerun/glat    --units MLA,ESA --name glat    --tag rtk
    python3 pipeline/compare_outputs.py results/cmaenas tmp/rerun/cmaenas --units YO,CNS  --name cmaenas --tag rtk

`compare_outputs.py` compares text columns (letters, calls, sample names) exactly and numeric
columns within a relative tolerance of 1e-6; p-values must also agree within 0.01 on the log10
scale. It prints one line per table (PASS, or FAIL with the differing cells) and exits with
status 0 only when every expected table passes.

## Attribution

Species silhouettes are from PhyloPic under the CC0 1.0 public-domain dedication:
*Gecarcinus* by Margot Michaud and *Carcinus maenas* by Mario Quevedo
(`species_photos/silhouettes/*.txt`).

## License

Copyright (c) 2026 Jorge L. Pérez-Moreno.

The code (`pipeline/`, `configs/`) is released under the MIT License (`LICENSE`). The sample
sheets, result tables, tables and figures (`data/`, `results/`, `tables/`, `figures/`,
`env_explicit.txt`) are released under the Creative Commons Attribution 4.0 International
license (`LICENSE-data`). The species silhouettes keep their CC0 1.0 dedication (see
"Attribution").
