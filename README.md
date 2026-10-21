# Receptor tyrosine kinase expression across the crustacean molt cycle

This repository holds the data, code and results behind the accompanying paper on receptor tyrosine
kinase (RTK) expression in the Y-organ and central nervous system of *Gecarcinus lateralis* and
*Carcinus maenas*. It contains everything needed to regenerate every table and figure of the paper
(Figures 1 to 16, Figures S1 to S5, Tables S1 and S2) except the per-library Salmon quantifications,
which are deposited on Zenodo (see "Step 3").

The result tables are in `results/` and `tables/`, and the figures are in `figures/`, so nothing
needs to be run to inspect the results. This page summarizes the study design, defines the
bioinformatic and statistical terms used, and gives the commands to rerun the analysis.

## Background

Crustaceans grow through a series of molts, in which the old exoskeleton is shed and replaced by a
larger one. In decapods, molting is controlled by the Y-organ (YO), the molting gland, which
synthesizes and secretes the steroid molting hormones (ecdysteroids) that initiate and coordinate
the molt cycle.

Insulin and growth factor receptors belong to the receptor tyrosine kinase (RTK) superfamily. Each
receptor has an extracellular ligand-binding domain, a single-pass transmembrane domain and an
intracellular tyrosine kinase domain. Ligand binding activates the kinase, and the main downstream
targets are the Ras/Raf/MAP kinase and PI3K/Akt/mTOR signaling pathways. Crustaceans express four
major types of these receptors: the insulin receptor (InsR), the epidermal growth factor receptor
(EGFR), the fibroblast growth factor receptor (FGFR) and the PDGF/VEGF-related receptor (PVR).

In insects, insulin-like peptide and growth factor RTK signaling stimulates ecdysteroid synthesis in
the prothoracic gland, the insect molting gland. The direct effects of these ligands on
ecdysteroidogenesis in the decapod YO have not been measured. The RTKs are expressed in the YO,
which suggests that growth factors contribute to YO regulation over the molt cycle.

This study uses RNA sequencing (RNA-seq) of the YO, and in *C. maenas* also of the central nervous
system (CNS), to test whether RTK transcript abundance changes over the molt cycle. For each RTK, in
each tissue and experiment, the analysis asks two questions.

1. Does its expression change across the molt cycle?
2. If it does, between which stages?

### The animals and samples

- ***Gecarcinus lateralis*** (the blackback land crab): 115 YO libraries from two molt-induction
  experiments, analyzed separately.
  - **Multiple leg autotomy (MLA).** Autotomy of the walking legs induces entry into premolt. 54
    animals were sampled at five molt stages, intermolt (IM), early premolt (EP), mid-premolt (MP),
    late premolt (LP) and postmolt (PM), and as blocked animals (BL), which were refractory to molt
    induction and did not enter premolt. Stages were assigned from the regeneration index (R, the
    length of the limb regenerate as a percentage of carapace width), the hemolymph ecdysteroid
    titer and the condition of the exoskeleton.
  - **Eyestalk ablation (ESA).** Bilateral eyestalk ablation removes the principal source of
    molt-inhibiting hormone (MIH), which activates the YO, and the animals enter premolt. 59 males
    were sampled at days 0, 1, 3, 5, 7, 14 and 19 after ablation. The intact animals are Day 0
    controls.
  - The *G. lateralis* panel has 17 RTK transcripts, one per gene.
- ***Carcinus maenas*** (the green shore crab): 25 YO and 25 CNS libraries from paired tissues of
  spontaneously molting animals, staged by setal development and grouped into four molt stages (IM,
  MP, LP and PM). The panel has 29 RTK transcripts from 23 genes. When a gene has several
  transcripts, their counts are summed and the gene is analyzed as one feature.

## Key terms

| Term | Meaning |
| :--- | :--- |
| Library | The sequencing sample made from one animal's tissue (for the YO, the paired YOs of one animal). Each library is one biological replicate in the analysis. |
| Salmon, `quant.sf` | Salmon estimates how many reads come from each transcript of a reference transcriptome, using selective alignment and correcting for sequence-composition and GC-content bias. It writes one `quant.sf` file per library, with estimated read counts and TPM for every transcript. The analysis starts from these files. |
| tximport | The R package that reads the `quant.sf` files into matrices of estimated counts and TPM (transcripts by libraries) and adjusts each library for differences in average transcript length. For *C. maenas* it also sums the transcripts of each gene. |
| TPM | Transcripts per million, an abundance measure corrected for transcript length and library size. TPM values are proportions of each library, so a large change in a few abundant transcripts shifts the TPM of every other transcript in that library. |
| Size factor (median-of-ratios) | A scaling factor per library. Each transcript's value is divided by its geometric mean across libraries, and the factor is the median of these ratios. On counts it corrects for sequencing depth. On TPM it corrects for the composition effect described above. |
| Normalized TPM | TPM divided by the library's median-of-ratios size factor, computed from the TPM of all transcripts. The bar charts in the paper show normalized TPM. |
| DESeq2, negative binomial GLM | DESeq2 models the read counts of each transcript with a negative binomial generalized linear model (GLM), a regression model for count data that allows for variation among animals beyond sampling noise. Molt stage is the single explanatory factor. |
| Dispersion | The variation in expression among animals at the same stage, beyond sampling noise. DESeq2 estimates it for every transcript and shrinks each estimate toward the values of transcripts with similar mean expression, which stabilizes it when each stage has few animals. |
| Likelihood-ratio test (LRT) | The test for question 1. It compares a full model, with one mean per stage, with a reduced model in which the stages share one mean. The degrees of freedom are the number of stages in the test minus one. A small p-value means that the stage means are not all equal. |
| Wald test | The test for question 2. It tests the difference between two stage means, the log2 fold change, against its standard error. |
| False discovery rate (FDR), Benjamini-Hochberg, q | The Benjamini-Hochberg procedure adjusts a set of p-values so that the expected proportion of false positives among the tests called significant (the FDR) stays at or below the threshold. The adjusted p-values are reported as q. An RTK is molt-regulated when `q_targeted` is below 0.05. |
| Targeted and transcriptome-wide FDR | `q_targeted` adjusts across the RTK panel only and decides the call. `q_transcriptomewide` is DESeq2's adjustment across all transcripts in the model, reported for comparison. |
| Holm procedure, family-wise error rate | The Holm procedure adjusts the pairwise p-values of one RTK so that the probability of even one false positive among its comparisons (the family-wise error rate) stays at or below 5%. |
| Compact letter display | The letters above the bars. Stages that share a letter are not significantly different, and stages with no letter in common are. |
| Log2 fold change, shrinkage | The log2 ratio of expression between two groups (1 means twice as high, -1 means half). Shrinkage with a normal prior pulls imprecise fold changes of low-count transcripts toward zero and leaves well-measured ones almost unchanged. |
| Cook's distance | A measure of how strongly a single library's count influences the fitted means of a transcript. DESeq2 uses it to detect outlier counts. |
| Variance-stabilizing transformation (VST), PCA | VST transforms the counts so that their variance depends little on mean expression. The principal component analysis (PCA) of the VST counts summarizes how similar the libraries are overall. |
| Volcano and MA plots | Two standard ways to show every transcript in one comparison: log2 fold change against significance (volcano), and log2 fold change against mean expression (MA). |
| Configuration file | A plain-text YAML file (`configs/*.yaml`) that gives the pipeline the data locations and the analysis settings of each dataset. |

## What is in this repository

    configs/                 one settings file per species (glat.yaml, cmaenas.yaml)
    data/glat/               G. lateralis sample sheet (samples.tsv) and RTK list (glat_rtk_genes.tsv)
    data/cmaenas/            C. maenas sample sheet (cmaenas_samplesheet.tsv) and RTK list (cmaenas_panel_reps.tsv)
    pipeline/                the analysis code (dge_panel_pipeline 1.0.6); pipeline/README.md explains it in depth
    results/glat/            result tables for the MLA and ESA experiments
    results/cmaenas/         result tables for the YO and CNS
    figures/main/            Figures 1 to 16
    figures/supplementary/   Figures S1 to S5
    tables/                  Tables S1 and S2
    species_photos/          crab silhouettes used on the figures, with their credits

## Rerunning the analysis, step by step

Run the commands below in a terminal (Linux or macOS) from the top folder of this repository unless
a step says otherwise. Each dataset takes about ten minutes on a computer with four cores.

### Step 1. Get the repository

    git clone https://github.com/invertome/yo-rtk-molt-expression.git
    cd yo-rtk-molt-expression

### Step 2. Install the software

The statistics run in R and the figures are drawn in Python. Install conda (Miniconda or Miniforge)
if it is not already available. Then create the R environment, which installs R and every package at
the exact versions used for the paper:

    conda env create -f pipeline/env_r.yaml

This creates an environment named `dge_panel_r`. The table and figure scripts also need Python 3
with pandas, matplotlib, numpy, PyYAML and Pillow, for example:

    pip install pandas matplotlib numpy pyyaml pillow

### Step 3. Download the input files

The analysis starts from one Salmon `quant.sf` file per library. These files are too large for
GitHub and are on Zenodo ([Zenodo DOI]). The raw reads are in the NCBI Sequence Read Archive (SRA)
under BioProject [SRA BioProject accession]. Place each file in its own directory, named as in the
sample sheet:

    data/glat/salmon_output/<filename>/quant.sf     <filename> = the "filename" column of data/glat/samples.tsv
    data/cmaenas/salmon_output/<run>/quant.sf       <run> = the "run" column of data/cmaenas/cmaenas_samplesheet.tsv

### Step 4. Point the pipeline to the R environment

Copy the template of the machine-specific settings and edit the copy:

    cp pipeline/site.env.example site.env

Set `CONDA_SH` to the path of conda's initialization script, which is `<conda
folder>/etc/profile.d/conda.sh` (`conda info --base` prints the conda folder), and set
`CONDA_ENV=dge_panel_r`. The file also sets `OPENBLAS_CORETYPE=HASWELL`, which pins one OpenBLAS
kernel so that the numbers are identical on different computers (see "Software"). That setting
applies to Intel and AMD processors only. On any other processor (for example ARM), comment it out
and remove the `#` in front of `DGE_BLAS_PIN=none`.

This step can be skipped if `Rscript` from the `dge_panel_r` environment is already on your `PATH`
(for example after `conda activate dge_panel_r`). The pipeline then runs with the BLAS kernel
unpinned and says so.

### Step 5. Run the statistics

Run each species once:

    bash pipeline/run_pipeline.sh configs/glat.yaml
    bash pipeline/run_pipeline.sh configs/cmaenas.yaml

A successful run ends with `VALIDATION PASSED` and `PIPELINE DONE` and writes its tables to
`results/glat/` or `results/cmaenas/`. A failed check stops the run with a message that names the
problem.

On a computing cluster that uses Slurm, `pipeline/run_pipeline.sbatch` runs the same script as a
job, and `pipeline/README.md` ("Running") shows how.

### Step 6. Build the tables and figures

    export BASE_LOCAL=.
    python3 pipeline/table_panel.py configs/glat.yaml    results/glat    tables
    python3 pipeline/table_panel.py configs/cmaenas.yaml results/cmaenas tables
    python3 pipeline/fig_main.py    configs/glat.yaml    results/glat    figures/main
    python3 pipeline/fig_main.py    configs/cmaenas.yaml results/cmaenas figures/main
    python3 pipeline/fig_supp.py    configs/glat.yaml    results/glat    figures/supplementary
    python3 pipeline/fig_supp.py    configs/cmaenas.yaml results/cmaenas figures/supplementary

### Step 7. Check the result

A correct rerun gives the same numbers of molt-regulated RTKs as the paper: 4 of the 17 RTK
transcripts in the MLA experiment, 6 of 17 in the ESA experiment, 4 of the 23 RTK genes in the *C.
maenas* YO, and 0 of 23 in the *C. maenas* CNS. "Checking a rerun" below compares every table, value
by value.

## How the analysis works

`pipeline/README.md` gives the full technical description of the pipeline, and the comments in
`configs/glat.yaml` and `configs/cmaenas.yaml` explain each setting for these datasets.

### Four analysis units

The data form four analysis units, each fitted and tested independently: the MLA and ESA experiments
in *G. lateralis*, and the YO and CNS in *C. maenas*. Before any model is fitted, the pipeline
checks the number of libraries at every stage against the number given in the configuration and
stops if they differ.

### From Salmon estimates to a count matrix

tximport reads the Salmon `quant.sf` files into a matrix of estimated counts with one row per
transcript and one column per library, and adjusts each library for differences in average
transcript length. For *C. maenas*, the transcripts of each gene are summed at this step so that
each gene is one row.

Transcripts with few reads are removed before the model is fitted. A transcript is retained when at
least as many libraries as the smallest stage group have 10 or more counts for it. This filter
improves the dispersion estimates and reduces the number of transcripts in the transcriptome-wide
correction for multiple testing. The RTK transcripts are exempt from the filter, so that no RTK is
excluded for low expression.

### Counts for testing, normalized TPM for display

The statistical tests use read counts, because the negative binomial model in DESeq2 is built for
counts. A count of 5 reads carries much less precision than a count of 5,000, and the model uses
that information. TPM values lose it, so they are not used for testing.

Counts depend on sequencing depth, which makes them hard to compare on a graph, so the figures show
normalized TPM. TPM corrects for library size and transcript length, but it is a proportion of each
library. When a few transcripts are extremely abundant in some libraries, they take a larger share
of the "per million" and lower the TPM of every other transcript. To remove this composition effect,
each library receives a median-of-ratios size factor computed with DESeq2 from the TPM of all
transcripts, and each TPM value is divided by its library's factor.

### The DESeq2 model

DESeq2 is fitted transcriptome-wide, to every transcript that passed the filter, and not only to the
RTKs. The counts of each transcript are modeled with a negative binomial GLM in which molt stage is
the single factor. Median-of-ratios size factors on the counts correct for differences in sequencing
depth between libraries. The full model gives every stage its own mean.

Each stage has few animals, too few to estimate the dispersion of one transcript reliably on its
own. DESeq2 therefore shrinks each transcript's dispersion estimate toward the values of transcripts
with similar expression levels. Fitting the whole transcriptome provides the thousands of
transcripts that this sharing of information requires.

### Question 1: does expression change across the molt cycle?

The likelihood-ratio test (LRT) compares the full model with a reduced model in which all stages
share one mean. If the full model fits the counts significantly better, expression differs among
stages somewhere in the cycle. The degrees of freedom equal the number of stages in the test minus
one, and the pipeline stops if the two models do not differ by that number of parameters. The test
gives one p-value per transcript and does not identify which stages differ, which is question 2.

In the MLA experiment the blocked animals are not part of the molt cycle, so they keep their own
mean in both the full and the reduced model. The LRT then asks only whether expression varies among
the five stages of the cycle (IM, EP, MP, LP, PM), while the blocked animals still contribute to the
dispersion estimates.

DESeq2 uses Cook's distance to detect a count in one animal that lies so far from the other animals
of its stage that it alone could drive the result. By default, in stages with seven or more animals
the outlier count is replaced with a typical value for that transcript estimated from all animals,
and in stages with fewer animals the transcript is left out of the test. The replacement is kept
here and the exclusion is switched off (`cooksCutoff=FALSE`), so no RTK goes untested because of one
animal.

### Correcting for multiple testing

When many transcripts are tested, some p-values fall below 0.05 by chance alone. The
Benjamini-Hochberg procedure adjusts the p-values to control the false discovery rate (FDR), the
expected proportion of false positives among the transcripts called significant. The adjusted values
are reported as q.

Two adjustments are reported. `q_targeted` applies Benjamini-Hochberg across the RTK panel only (17
tests for *G. lateralis*, 23 for *C. maenas*). The panel was defined before the analysis,
independently of these data, and the RTKs are the only transcripts the study asks about, so this
targeted FDR decides the call. An RTK is **molt-regulated** when `q_targeted` is below 0.05.
`q_transcriptomewide` is DESeq2's Benjamini-Hochberg adjustment across all transcripts in the model.
It covers many more tests and is therefore usually larger. It is reported alongside, in the result
tables and in Tables S1 and S2, and does not affect the call.

### Question 2: which stages differ?

For each molt-regulated RTK, every pair of stages, including the blocked animals in the MLA
experiment, is compared with a Wald test on the full model. Because one RTK gets many pairwise
comparisons, the p-values are adjusted again within that RTK with the Holm procedure. Holm controls
the family-wise error rate, keeping the probability of even one false positive among that RTK's
comparisons at or below 5%.

The results are summarized as a compact letter display above the bars. Two stages that share a
letter are not significantly different (Holm-adjusted p of 0.05 or more), and two stages with no
letter in common are. For example, stages labeled "a", "ab" and "b" mean that the first and last
differ and the middle one differs from neither. An RTK that is not molt-regulated has the same
letter at every stage in the letters table and no letters on the figures.

### Supplementary comparisons

For Figures S1 to S5, the pipeline compares neighboring time points across the transcriptome with
Wald tests: each step of the cycle (IM to EP, EP to MP, and so on, including PM back to IM), each
day after ablation against the previous day, and the blocked animals against every MLA stage. Each
comparison is shown as a volcano plot and an MA plot, with the RTKs in red and starred when their
transcriptome-wide adjusted p-value in that comparison is below 0.05.

The log2 fold changes in these tables and plots are shrunken with DESeq2's `lfcShrink` and a normal
prior (`shrinkage: normal` in the configurations). A transcript with few reads can show a large fold
change from a small random difference in counts. Shrinkage pulls such imprecise estimates toward
zero and leaves well-measured ones almost unchanged, so low-count transcripts do not dominate the
plots. The p-values are those of the Wald test. The LRT and the compact letters do not use the
shrunken fold changes.

Each run also writes a principal component analysis (PCA) of each unit. It uses the 500 most
variable transcripts after a variance-stabilizing transformation (VST) of the counts, computed
without reference to the stages, and its first two components summarize how similar the libraries
are overall. The PCA is descriptive and plays no part in the tests.

### Reading the main figures

Each panel of Figures 1 to 16 is one RTK. The gray bars are the mean normalized TPM at each stage,
the error bars are one standard error of the mean (SEM, the standard deviation divided by the square
root of the number of animals), and each dot is one animal. The ESA figures draw a line through the
days after ablation instead of bars. Letters appear only above molt-regulated RTKs.

### Automatic checks

After every run, a validation script checks the output tables for internal consistency: one row per
transcript and stage, sample counts that match the sample sheet, molt-regulated calls that match the
adjusted p-values, and an explanation for every missing test result. If any check fails, the run
stops with an error instead of producing tables that might be wrong.

## Reference: input files

The details below are for readers who want to reuse the data or adapt the pipeline.

**Quantifications.** One Salmon `quant.sf` file per library, placed as in "Step 3". The *G.
lateralis* libraries were quantified with Salmon 1.10.2 and the *C. maenas* libraries with Salmon
1.10.3, both with automatic library-type detection.

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
transcripts of 23 genes by gene, transcript identifier and family. Each of the six genes with
several transcripts is quantified as the sum of its transcripts.

## Reference: output files

**Result tables** in `results/glat/` (units `MLA` and `ESA`, in the unit column `experiment`) and
`results/cmaenas/` (units `YO` and `CNS`, in the unit column `tissue`). `<unit>` stands for one of
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
track. The pairwise comparisons are MLA `EP_vs_IM`, `MP_vs_EP`, `LP_vs_MP`, `PM_vs_LP`,
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

In Figures 1 to 16, each panel shows one RTK transcript. Gray bars give the mean normalized TPM at
each stage with its standard error, dots show one library each, and compact letters sit over the
molt-regulated transcripts. The ESA figures draw a line over days instead of bars. In Figures S1 to
S5, each row is one pairwise comparison, with a volcano plot on the left (log2 fold change against
-log10 adjusted p-value) and an MA plot on the right (log2 fold change against mean expression).
Every transcript is a gray point, RTK transcripts are red, and those significant in that comparison
are starred and named. `figures/supplementary/glat_stars.json` and `cmaenas_stars.json` list, for
each supplementary figure, the (comparison, gene) pairs that received a star, which can be checked
against the `rtk_sig` column of the result tables.

**Tables.** `tables/TableS1_Gl` (Table S1, *G. lateralis*) and `tables/TableS2_Cm` (Table S2,
*C. maenas*), each as `.tsv` and `.md`, list every RTK transcript in every unit with `gene`,
`family`, the unit column, `LRT_p`, `q_targeted`, `q_transcriptomewide` and `molt_regulated`. The
Markdown version writes the numbers in scientific notation.

## Software

R 4.5.3 with DESeq2 1.50.2, tximport 1.38.2, readr 2.2.0, data.table 1.18.6.1, multcompView 0.1-12
and yaml 2.3.12 (`pipeline/env_r.yaml`), and `env_explicit.txt` lists the exact environment used.
The table and figure scripts need Python 3 with pandas, matplotlib, numpy, PyYAML and Pillow.

**Numerical reproducibility.** The tables in `results/` were computed with the OpenBLAS kernel
pinned (`OPENBLAS_CORETYPE=HASWELL`, set in `site.env`). Keep that setting to reproduce them to the
last digit. OpenBLAS builds that select a CPU-specific kernel at run time round differently on
different processor generations, which changes the last decimal places of some transcriptome-wide
coefficient estimates. The figures use the Liberation fonts. If those are missing, matplotlib uses
another typeface and the plotted data stay the same.

## Checking a rerun

To confirm that a rerun reproduces the tables in `results/`, write the rerun to a separate folder
and compare the two, table by table:

    PIPELINE_OUT=tmp/rerun/glat    bash pipeline/run_pipeline.sh configs/glat.yaml
    PIPELINE_OUT=tmp/rerun/cmaenas bash pipeline/run_pipeline.sh configs/cmaenas.yaml
    python3 pipeline/compare_outputs.py results/glat    tmp/rerun/glat    --units MLA,ESA --name glat    --tag rtk
    python3 pipeline/compare_outputs.py results/cmaenas tmp/rerun/cmaenas --units YO,CNS  --name cmaenas --tag rtk

`compare_outputs.py` prints one line per table: PASS, or FAIL with the cells that differ. Text
columns (letters, calls, sample names) must match exactly. Numbers must agree within a relative
difference of 1e-6, and p-values also within 0.01 on the log10 scale. It exits with status 0 only
when every table passes. It does not compare the normalized-TPM columns or the size-factor files.

## Attribution

Species silhouettes are from PhyloPic under the CC0 1.0 public-domain dedication:
*Gecarcinus* by Margot Michaud and *Carcinus maenas* by Mario Quevedo
(`species_photos/silhouettes/*.txt`).

## License

Copyright (c) 2026 Jorge L. Pérez-Moreno, Ph.D.

The code (`pipeline/`, `configs/`) is released under the MIT License (`LICENSE`). The sample
sheets, result tables, tables and figures (`data/`, `results/`, `tables/`, `figures/`,
`env_explicit.txt`) are released under the Creative Commons Attribution 4.0 International
license (`LICENSE-data`). The species silhouettes keep their CC0 1.0 dedication (see
"Attribution").
