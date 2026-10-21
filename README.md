# Receptor tyrosine kinase expression across the crustacean molt cycle

This repository holds the data, code and results behind the accompanying paper on receptor tyrosine
kinase (RTK) expression in crabs. Everything needed to rebuild every table and figure of the paper
(Figures 1 to 16, Figures S1 to S5, Tables S1 and S2) is here, except the large per-sample input
files, which are on Zenodo (see "Step 3").

You do not need to run anything to look at the results: the finished tables are in `results/` and
`tables/`, and the figures are in `figures/`. The rest of this page explains the study in plain
terms, defines the technical words, and walks through rerunning the analysis one step at a time.

## The study in plain terms

Crabs grow by molting: they shed their hard shell and grow a new, larger one. Molting is
controlled by a small gland, the **Y-organ** (YO), which makes the molting hormones
(ecdysteroids). The YO's activity is controlled by signals from the rest of the body.

**Receptor tyrosine kinases** are proteins that sit in the cell membrane and pass signals from
outside the cell to the inside. Four RTK families are studied here: insulin receptors (InsR),
epidermal growth factor receptors (EGFR), fibroblast growth factor receptors (FGFR), and
PDGF/VEGF-related receptors (PVR). If the YO changes how much of these receptors it makes during
the molt cycle, that is a hint that they help control molting.

So for each RTK, in each tissue and experiment, the analysis asks two questions:

1. Does its expression change across the molt cycle?
2. If it does, between which stages?

### The animals and samples

- ***Gecarcinus lateralis*** (the blackback land crab): 115 Y-organ libraries from two
  experiments, analyzed separately.
  - **Multiple leg autotomy (MLA).** Removing several legs makes a crab start a new molt cycle.
    54 animals were sampled at six stages: intermolt (IM), early premolt (EP), mid-premolt (MP),
    late premolt (LP), postmolt (PM), and blocked animals (BL) that did not enter premolt. Stages
    were assigned from the regeneration index (how far the new limb buds have grown) and the
    amount of molting hormone in the blood (hemolymph ecdysteroid titer).
  - **Eyestalk ablation (ESA).** Removing the eyestalks removes a hormone that holds the YO back,
    so molting starts. 59 males were sampled at days 0, 1, 3, 5, 7, 14 and 19 after ablation; the
    intact animals are the Day 0 controls.
  - The *G. lateralis* gene panel has 17 RTK transcripts, one per gene.
- ***Carcinus maenas*** (the green crab): 25 YO and 25 central nervous system (CNS) libraries
  from paired tissues, at four molt stages: IM, MP, LP and PM. The panel has 29 RTK
  transcripts from 23 genes; when a gene has several transcripts, they are added together and the
  gene is analyzed as one.

## Key terms

| Term | Meaning |
| :--- | :--- |
| Gene, transcript | A gene is a stretch of DNA; a transcript is an RNA copy of it. One gene can give several slightly different transcripts. Measuring transcripts tells us how actively a gene is used. |
| RNA-seq, library | RNA sequencing reads millions of short pieces of RNA from a tissue sample. The set of reads from one sample is called a library. |
| Salmon | A program that counts how many reads came from each transcript. Its output, one `quant.sf` file per library, is where this analysis starts. |
| TPM | Transcripts per million: how much of a transcript there is, scaled so that libraries of different sizes can be compared. |
| Normalized TPM | TPM corrected for differences in overall makeup between libraries (each library's TPM is divided by a "size factor"). The bar charts in the paper show normalized TPM. |
| DESeq2 | A widely used R package for testing whether expression differs between groups of samples. It works on read counts and accounts for differences in library size. |
| Likelihood-ratio test (LRT) | The test for question 1. It compares a model in which each stage has its own average with one in which all stages share a single average. A small p-value means expression changes across the cycle. |
| Wald test | The test for question 2: it compares one stage with another. |
| False discovery rate, q-value | When many genes are tested at once, some will look significant by chance. Adjusted p-values (q-values) correct for this; a transcript is called molt-regulated when its q-value is below 0.05. |
| Compact letters | The letters above the bars. Stages that share a letter are not significantly different; stages with no letter in common are. |
| Log2 fold change | How many times higher or lower expression is in one group than another, on a log2 scale (1 means twice as high, -1 means half). |
| Volcano and MA plots | Two standard ways to show many genes at once in a comparison: fold change against significance (volcano), and fold change against average expression (MA). |
| conda | A tool that installs programs and their exact versions into a separate "environment", so they do not conflict with other software on your computer. |
| Configuration file | A plain-text file (`configs/*.yaml`) that tells the pipeline where the data are and which choices to use for each dataset. |

## What is in this repository

    configs/                 one settings file per species (glat.yaml, cmaenas.yaml)
    data/glat/               G. lateralis sample sheet (samples.tsv) and RTK list (glat_rtk_genes.tsv)
    data/cmaenas/            C. maenas sample sheet (cmaenas_samplesheet.tsv) and RTK list (cmaenas_panel_reps.tsv)
    pipeline/                the analysis code (dge_panel_pipeline 1.0.5); pipeline/README.md explains it in depth
    results/glat/            result tables for the MLA and ESA experiments
    results/cmaenas/         result tables for the YO and CNS
    figures/main/            Figures 1 to 16
    figures/supplementary/   Figures S1 to S5
    tables/                  Tables S1 and S2
    species_photos/          crab silhouettes used on the figures, with their credits

## Rerunning the analysis, step by step

The commands below are typed into a terminal (Linux or macOS), from the top folder of this
repository unless a step says otherwise. Each dataset takes about ten minutes on a computer with
four cores.

### Step 1. Get the repository

    git clone https://github.com/invertome/yo-rtk-molt-expression.git
    cd yo-rtk-molt-expression

### Step 2. Install the software

The statistics run in R and the figures are drawn in Python. If you do not have conda yet,
install Miniconda or Miniforge first (both are free). Then create the R environment, which
installs R and every package at the exact versions used for the paper:

    conda env create -f pipeline/env_r.yaml

This makes an environment named `dge_panel_r`. For the figures and tables you also need Python 3
with pandas, matplotlib, numpy, PyYAML and Pillow, for example:

    pip install pandas matplotlib numpy pyyaml pillow

### Step 3. Download the input files

The analysis starts from one Salmon `quant.sf` file per library. These files are too large for
GitHub and are on Zenodo ([Zenodo DOI]); the raw reads are in the NCBI Sequence Read Archive under
BioProject [SRA BioProject accession]. Put each file in its own folder, named as in the sample
sheet:

    data/glat/salmon_output/<filename>/quant.sf     <filename> = the "filename" column of data/glat/samples.tsv
    data/cmaenas/salmon_output/<run>/quant.sf       <run> = the "run" column of data/cmaenas/cmaenas_samplesheet.tsv

### Step 4. Tell the pipeline where R is

Copy the example settings file and open the copy in a text editor:

    cp pipeline/site.env.example site.env

Fill in two lines: `CONDA_SH`, the path of conda's start-up script (it is
`<conda folder>/etc/profile.d/conda.sh`, and `conda info --base` prints the conda folder), and
`CONDA_ENV=dge_panel_r`. The file also sets `OPENBLAS_CORETYPE=HASWELL`, which makes the numbers
come out identical on different computers (see "Software"). That setting applies to Intel and AMD
processors only; on any other processor (for example ARM), comment it out and remove the `#` in
front of `DGE_BLAS_PIN=none`.

You can skip this step if `Rscript` from the `dge_panel_r` environment is already on your `PATH`
(for example after `conda activate dge_panel_r`); the pipeline then says that the BLAS setting is
unpinned.

### Step 5. Run the statistics

Run each species once:

    bash pipeline/run_pipeline.sh configs/glat.yaml
    bash pipeline/run_pipeline.sh configs/cmaenas.yaml

A run that worked ends with `VALIDATION PASSED` and `PIPELINE DONE`, and its tables are in
`results/glat/` or `results/cmaenas/`. If it stops with an error, the last lines of the output
show what went wrong.

On a computing cluster that uses Slurm, `pipeline/run_pipeline.sbatch` runs the same script as a
job; `pipeline/README.md` ("Running") shows how.

### Step 6. Build the tables and figures

    export BASE_LOCAL=.
    python3 pipeline/table_panel.py configs/glat.yaml    results/glat    tables
    python3 pipeline/table_panel.py configs/cmaenas.yaml results/cmaenas tables
    python3 pipeline/fig_main.py    configs/glat.yaml    results/glat    figures/main
    python3 pipeline/fig_main.py    configs/cmaenas.yaml results/cmaenas figures/main
    python3 pipeline/fig_supp.py    configs/glat.yaml    results/glat    figures/supplementary
    python3 pipeline/fig_supp.py    configs/cmaenas.yaml results/cmaenas figures/supplementary

### Step 7. Check the result

The paper reports these numbers of molt-regulated RTKs; a correct rerun gives the same:
4 of the 17 transcripts in the MLA experiment, 6 of 17 in the ESA experiment, 4 of the 23 genes
in the *C. maenas* YO, and 0 of 23 in the *C. maenas* CNS. "Checking a rerun" below compares
every table, value by value.

## How the analysis works

This section summarizes the method; `pipeline/README.md` gives the full details, and the comments
in `configs/glat.yaml` and `configs/cmaenas.yaml` explain each setting.

1. **Four separate analyses.** MLA and ESA for *G. lateralis*, YO and CNS for *C. maenas*. Before
   any statistics, the pipeline checks that the number of samples at every stage matches the
   configuration.
2. **One model for all genes.** DESeq2 is fitted to every transcript in the sample (for
   *C. maenas*, every gene), not only the RTKs, with one average per stage. Using all genes gives
   DESeq2 better estimates of how variable expression is. RTKs stay in the model even when
   barely expressed.
3. **Does expression change over the cycle?** The likelihood-ratio test compares that model with
   one in which all stages share a single average. In the MLA experiment the blocked animals keep
   their own average in both models, so the test asks only whether expression varies across the
   five stages of the actual cycle. DESeq2 can automatically discard a gene in which one animal
   has an extreme value (detected with Cook's distance); that option is switched off here, so a
   single unusual animal cannot remove an RTK from the analysis.
4. **Which RTKs are molt-regulated.** The p-values of the RTKs are adjusted for multiple testing
   across the RTK panel only (Benjamini-Hochberg over the 17 or 23 RTKs). An RTK is molt-regulated
   when that adjusted value (`q_targeted`) is below 0.05. The value adjusted across all genes
   (`q_transcriptomewide`) is reported as well but does not decide the call.
5. **Which stages differ.** For each molt-regulated RTK, every pair of stages is compared with a
   Wald test, the p-values are adjusted within that RTK (Holm method), and the result is shown as
   compact letters above the bars. RTKs that are not molt-regulated get no letters on the figures.
6. **Supplementary comparisons.** Each step of the cycle (IM to EP, EP to MP, and so on, including
   PM back to IM), each day after ablation against the day before, and the blocked animals against
   every MLA stage are compared for all genes, and shown as volcano and MA plots with the RTKs
   marked (Figures S1 to S5). A principal component analysis (a summary of overall similarity
   between samples) is saved for each analysis too.
7. **Automatic checks.** After every run the output tables are checked for consistency (one row
   per transcript and stage, sample counts that match the sheet, calls that match the adjusted
   p-values, every missing value explained). If a check fails, the run stops with an error rather
   than producing doubtful tables.

## Reference: input files

The details below are for readers who want to reuse the data or adapt the pipeline.

**Quantifications.** One Salmon `quant.sf` file per library, placed as in "Step 3". The *G. lateralis* libraries were quantified with Salmon 1.10.2 and the *C. maenas* libraries
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

## Reference: output files

**Result tables** in `results/glat/` (units `MLA` and `ESA`; the unit column is `experiment`) and
`results/cmaenas/` (units `YO` and `CNS`; the unit column is `tissue`). `<unit>` stands for one of
those names.

| File | Columns |
| :--- | :--- |
| `<unit>_rtk_letters.tsv` | one row per RTK transcript and stage: `gene`; the unit column; `stage`; `letter` (compact letter at that stage); `lrt_p` (raw p-value of the likelihood-ratio test); `q_targeted` (adjusted across the panel; the call is made on it); `q_transcriptomewide` (DESeq2's adjusted p-value across all transcripts, the transcriptome-wide false discovery rate); `molt_regulated` (`TRUE` when `q_targeted` is below 0.05). The test values repeat on every stage row of a transcript. |
| `<unit>_rtk_tpm_by_stage.tsv` | `gene`; `id` (the transcript, or the summed transcripts joined by `;`); the unit column; `stage`; `n` (libraries at that stage); `mean_tpm`; `sem_tpm` (standard error of the mean); `mean_ntpm`; `sem_ntpm` (mean and standard error of normalized TPM) |
| `glat_rtk_persample_tpm.tsv`, `cmaenas_rtk_persample_tpm.tsv` | the transcripts per million (TPM) of every RTK transcript in every library: `gene`; `id`; the unit column; `sampleid`; `stage`; `tpm`; `tpm_sizefactor`; `ntpm` (normalized TPM, `tpm / tpm_sizefactor`) |
| `<unit>_rtk_tpm_sizefactors.tsv` | `sampleid`; `stage`; `tpm_sizefactor`: one row per library of the unit, the median-of-ratios size factor computed by DESeq2 (`estimateSizeFactorsForMatrix`) on the TPM of all transcripts |
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

In Figures 1 to 16, each panel shows one RTK transcript: gray bars give the mean normalized TPM at
each stage with its standard error, one dot per library, and compact letters over the molt-regulated
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
0.1-12 and yaml 2.3.12 (`pipeline/env_r.yaml`); `env_explicit.txt` lists the exact environment
used. The table and figure scripts need Python 3 with pandas, matplotlib, numpy, PyYAML and Pillow.

**Why the numbers can differ in the last decimal places.** The tables in `results/` were computed
with the OpenBLAS math library pinned to one processor setting (`OPENBLAS_CORETYPE=HASWELL`, set
in `site.env`). Without it, OpenBLAS picks a setting to suit each processor, and different
processors round slightly differently, which changes the last decimal places of some
transcriptome-wide estimates. The figures use the Liberation
fonts; if those are missing, matplotlib uses another typeface and the plotted data stay the same.

## Checking a rerun

To confirm that a rerun reproduces the tables in `results/`, write the rerun to a separate folder
and compare the two, table by table:

    PIPELINE_OUT=tmp/rerun/glat    bash pipeline/run_pipeline.sh configs/glat.yaml
    PIPELINE_OUT=tmp/rerun/cmaenas bash pipeline/run_pipeline.sh configs/cmaenas.yaml
    python3 pipeline/compare_outputs.py results/glat    tmp/rerun/glat    --units MLA,ESA --name glat    --tag rtk
    python3 pipeline/compare_outputs.py results/cmaenas tmp/rerun/cmaenas --units YO,CNS  --name cmaenas --tag rtk

`compare_outputs.py` prints one line per table: PASS, or FAIL with the cells that differ. Text
columns (letters, calls, sample names) must match exactly; numbers must agree within a relative
difference of 1e-6, and p-values also within 0.01 on the log10 scale. It exits with status 0 only
when every table passes. It does not compare the normalized-TPM columns or the size-factor files.

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
