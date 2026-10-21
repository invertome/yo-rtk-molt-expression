# DGE-panel pipeline

A small, configuration-driven pipeline for one recurring question in expression studies: does
each gene in a panel of genes of interest change its expression across the stages of an
experiment? The stages may be the phases of a cycle, the days after a treatment, or any ordered
set of conditions. The pipeline fits DESeq2 models transcriptome-wide, tests the stage effect, reports
the panel genes in detail, and draws the figures and tables that go into a paper. One YAML file
per dataset holds every dataset-specific value; the same code runs for every dataset.

## How the analysis works

The steps below are what `run_pipeline.sh` executes, in order, for every analysis unit (a tissue,
an experiment) named in the configuration.

1. **Inputs are checked before anything is computed.** The sample sheet must list every library
   with its stage and unit; every stage value must be one of the configured levels or one of the
   values deliberately excluded; the number of libraries per stage must match the expected counts
   when they are given; every quantification file must exist; the panel must list each transcript
   once and each gene in one family. A problem stops the run with a message naming it.
2. **Quantifications are imported** with tximport from one Salmon `quant.sf` file per library.
   Depending on the configuration, each panel gene is either one transcript or the sum of several
   transcripts.
3. **A transcriptome-wide DESeq2 model is fitted** with one mean per stage. Features are kept when at
   least as many libraries as the smallest stage group have ten or more counts; panel features are
   always kept, so every panel gene is reported even when barely expressed.
4. **The stage effect is tested** with a likelihood-ratio test that compares the fitted model with
   a reduced model in which all the stages share one mean. Levels listed in `outside_stage_test`
   (for example animals prevented from progressing through a cycle) keep their own mean in the
   reduced model: they stay in the fit and inform the dispersion estimates, but the test asks
   only whether expression varies across the remaining stages. DESeq2's Cook's-distance outlier
   exclusion is switched off, so that a single high-count library cannot silently remove a gene
   from the test.
5. **Panel genes are screened at a targeted false discovery rate.** The raw p-values of the panel
   genes are adjusted with Benjamini-Hochberg across the panel only, and a gene is called
   regulated when that adjusted value is below `alpha`. The transcriptome-wide adjusted p-value is
   reported alongside for reference but does not decide the call.
6. **Regulated genes get compact letters.** Every pair of stages is compared with a Wald test, the
   p-values are adjusted within the gene (Holm by default), and the outcome is summarized as
   letters: two stages that share a letter are not significantly different. A gene that is not
   regulated carries the same letter at every stage.
7. **Supplementary tables** give, for every consecutive pair of stages (and the closing pair of a
   cycle, and any extra pairs), the transcriptome-wide shrunken log2 fold changes with their tests, plus a
   principal component analysis of the variance-stabilized counts.
8. **The output is validated**: every gene appears once at every stage, the sample counts agree
   with the sheet, the regulated calls agree with the adjusted p-values, and every gene without a
   test result has an identifiable reason. The size-factor files cover every library and agree
   with the per-sample factors. Normalized TPM equals plain TPM divided by its factor, and stage
   means and SEM agree with the per-sample normalized values (SEM is NA for a single library).
   The run is declared done only when these checks pass.
9. **Figures and tables** are then drawn from the output tables by the Python scripts: the main
   figures (mean and standard error of normalized TPM per stage, with the letters), the supplementary volcano
   and MA plots, and the significance table.

For each unit, normalized TPM uses median-of-ratios size factors computed by
`DESeq2::estimateSizeFactorsForMatrix` on the transcript-level TPM matrix of all transcripts.
With gene aggregation enabled, a separate transcript-level import supplies that matrix. Each
library's exported panel TPM is divided by its factor, matched by sample ID, before computing
stage means and SEM (`sd / sqrt(n)`). The main figures plot these normalized values; plain TPM
is still exported alongside them.

## Files

| File | Role |
| :--- | :--- |
| `run_pipeline.sh` | runs the three R scripts below for one configuration, on the current machine |
| `run_pipeline.sbatch` | Slurm wrapper around `run_pipeline.sh` |
| `dge_panel.R` | steps 1 to 6: checks, import, fits, test, targeted FDR, letters, TPM tables |
| `supp_contrasts.R` | step 7: pairwise tables and PCA, from the fits saved by `dge_panel.R` |
| `validate_outputs.R` | step 8: consistency checks on the output tables |
| `lib_dge.R` | shared functions: configuration checks, TPM size factors, the reduced model, the targeted FDR, the letters |
| `fig_main.py` | main figures, one per gene family and unit, in two styles |
| `fig_supp.py` | supplementary volcano and MA plots |
| `table_panel.py` | the significance table (TSV and Markdown) |
| `compare_outputs.py` | compares two output directories cell by cell, for example a rerun against an existing set of tables |
| `env_r.yaml` | conda environment for the R side, with every package pinned to a verified version |
| `site.env.example` | template for the machine-specific settings (R environment, BLAS kernel) |
| `configs/example.yaml` | template configuration with every key explained |
| `tests/` | unit tests of `lib_dge.R` (`Rscript tests/test_lib_dge.R`, from the package directory) |

## Requirements

- R with `DESeq2`, `tximport`, `data.table`, `multcompView`, and `yaml`. The conda environment in
  `env_r.yaml` provides them at the versions the pipeline was verified with:
  `conda env create -f env_r.yaml` creates it under the name `dge_panel_r`.
- Python 3 with `pandas`, `matplotlib`, `numpy`, `PyYAML`, and `Pillow` for the figures and tables.
- The Liberation fonts (Liberation Serif and Liberation Sans) reproduce the figures exactly;
  without them matplotlib falls back to Times New Roman, Arial, or DejaVu, and the typeface
  changes while the data drawn do not.

## Input files

**Sample sheet.** A tab-separated file with one row per sequenced library. The pipeline does not
fix the column names; the configuration says which column plays which role:

| Configuration key | The column holds |
| :--- | :--- |
| `sample_id_col` | the library identifier (used as the column name of the count matrix) |
| `quant_subdir_col` | the name of the library's subdirectory under `quant_dir`, which holds its `quant.sf` |
| `level_col` | the stage of the library |
| `unit_col` | the analysis unit the library belongs to |

Column names are lower-cased when the sheet is read. Each unit is fitted from the rows whose unit
column equals the unit's name. A row whose stage is listed in `exclude_levels` is dropped and
reported; a row whose stage is neither a configured level nor an excluded value stops the run.

**Panel.** A tab-separated file with the columns `gene`, `transcript_id`, and `family` (further
columns are ignored). With `panel_aggregate: false` every gene has exactly one row and is analyzed
as that transcript. With `panel_aggregate: true` a gene may have several rows, one per transcript,
and the transcripts are summed before the analysis; the transcripts summed for each gene are
written to `<name>_<panel_tag>_panel.tsv`.

## Configuration keys

`configs/example.yaml` carries every key with a comment. In brief:

### Top level

| Key | Meaning |
| :--- | :--- |
| `name` | short name of the dataset; used in the filenames of the tables that span all units, and as `--name` for `compare_outputs.py` |
| `base` | directory under which every relative path resolves (the `BASE` environment variable overrides it) |
| `quant_dir` | directory of Salmon outputs, relative to `base` |
| `sample_sheet` | the sample sheet, relative to `base` |
| `sample_id_col`, `quant_subdir_col`, `level_col`, `unit_col` | the roles of the sample-sheet columns (table above) |
| `panel` | the panel file, relative to `base` |
| `panel_aggregate` | `true` sums the transcripts of a gene; `false` uses one transcript per gene |
| `panel_tag` | one word used in the names of the panel tables (`<unit>_<panel_tag>_letters.tsv` and friends); default `panel` |
| `output_dir` | where the R scripts write, relative to `base` (`PIPELINE_OUT` overrides it) |
| `alpha` | significance threshold for the targeted false discovery rate and for the letters |
| `pairwise_adjust` | adjustment of the pairwise p-values within each regulated gene: `holm` or `BH` |
| `shrinkage` | fold-change shrinkage for the supplementary tables: `normal` (DESeq2's normal prior, the only method the environment installs) |
| `units` | one entry per analysis unit (next table) |
| `figures`, `tables` | layout of the figures and of the significance table (tables below) |

### Per unit (`units.<unit>.*`)

| Key | Meaning |
| :--- | :--- |
| `levels` | the stages, in the order they are tested and drawn |
| `outside_stage_test` | stages kept in the model with their own mean but left out of the stage test |
| `closing_contrast` | `true` adds the comparison of the first stage against the last (a cycle) to the supplement |
| `extra_contrasts` | further `[later, earlier]` pairs for the supplement |
| `exclude_levels` | stage values that are in the sheet on purpose but not analyzed |
| `expected_n` | libraries expected per stage; a mismatch with the sheet stops the run |
| `x_axis` | `categorical` (one bar per stage) or `numeric` (a line over a numeric axis) |
| `x_label` | axis title of the main figures |
| `numeric_from` | with a numeric axis: a regular expression whose first group gives the number in each stage name |

### `figures.*`

| Key | Meaning |
| :--- | :--- |
| `family_order` | gene families in the order the figures are numbered |
| `panel_label` | how the panel genes are called in the supplementary legends (default `panel gene`) |
| `styles` | figure styles to draw: `manuscript` (serif, gray bars, few ticks) and `enhanced` (sans-serif, one dot per animal) |
| `style_dirs` | output subdirectory of each style; `""` writes that style directly into the figure directory |
| `filename` | main-figure filename with the placeholders `n`, `family`, `unit`; `n` also accepts a Python format spec, for example for zero padding |
| `ncols` | panels per row for a family (default 2) |
| `figure_numbers` | figure number of each `family/unit` combination |
| `silhouette`, `silhouette_label` | image drawn in a corner of every main figure (relative to `base`) and its italic caption |
| `fixed_ylim` | y-axis maxima for the manuscript style, one per panel of a `family/unit`; panels without an entry scale to their data |
| `display_names` | gene names to show differently from the panel file |
| `supplement` | one supplementary figure per entry: `{unit, contrasts, tag, title}` and optionally `filename`, where `contrasts` is `transitions` (consecutive stages, plus the closing pair when configured) or `extra` (the unit's `extra_contrasts`), and `filename` names the output files without extension (default `Fig<tag>_<unit>_dge`) |

### `tables.*`

| Key | Meaning |
| :--- | :--- |
| `order` | row order of the significance table: `panel` (order of the panel file) or `family_name` (family, then display name) |
| `name` | filename stem of the significance table (`.tsv` and `.md`) |

## Running

**On one machine, without a scheduler.** Copy `site.env.example` to `site.env` in the directory you
run from and fill it in (which conda environment holds R, and the BLAS setting), or run without a
site file: the `Rscript` found on `PATH` is then used and the BLAS kernel is left unpinned, which
the runner says out loud. Then:

    bash run_pipeline.sh /path/to/config.yaml
    PIPELINE_OUT=fresh_dir bash run_pipeline.sh /path/to/config.yaml   # write the tables elsewhere

The runner prints `PIPELINE DONE` only after the validator printed `VALIDATION PASSED`. A run of a
dataset with 50 to 115 libraries takes about ten minutes on four cores.

**On a Slurm cluster.** Submit from the login host, with the data and the configuration already
under the configuration's `base` (the package never copies data). Account, partition and log
paths go on the command line; the variables the runner reads travel with `--export`:

    P=/path/to/dge_panel_pipeline
    sbatch --account=<account> --partition=<partition> \
      --output=<base>/logs/%x-%j.out --error=<base>/logs/%x-%j.err \
      --export=ALL,DGE_PIPELINE_DIR=$P,CONDA_SH=<path to conda.sh>,CONDA_ENV=dge_panel_r,OPENBLAS_CORETYPE=HASWELL \
      $P/run_pipeline.sbatch <base>/configs/<config>.yaml

Add `PIPELINE_OUT=<dir>` to the `--export` list to write elsewhere, and `BASE=<dir>` to run the
configuration against data kept elsewhere.

**Figures and tables**, after a run, from the project directory (`BASE_LOCAL` names it, so that
the panel file and the silhouette are found):

    export BASE_LOCAL=.
    python3 table_panel.py <config.yaml> <output_dir> <table_dir>
    python3 fig_main.py    <config.yaml> <output_dir> <figure_dir>
    python3 fig_supp.py    <config.yaml> <output_dir> <figure_dir>/supplementary

## Output files

The R scripts write these tables into `output_dir`, for every unit `<unit>` (with `<tag>` the
`panel_tag` and `<name>` the dataset name):

| File | Columns |
| :--- | :--- |
| `<unit>_<tag>_letters.tsv` | `gene`; the unit column; `stage`; `letter` (the compact letter of the gene at that stage); `lrt_p` (raw p-value of the likelihood-ratio test); `q_targeted` (that p-value adjusted across the panel, the number that decides the call); `q_transcriptomewide` (DESeq2's transcriptome-wide adjusted p-value); `molt_regulated` (`TRUE` when `q_targeted` is below `alpha`). The test statistics are repeated on every stage row of the gene. |
| `<unit>_<tag>_tpm_by_stage.tsv` | `gene`; `id` (the transcript, or the summed transcripts joined by `;`); the unit column; `stage`; `n` (libraries at that stage); `mean_tpm`; `sem_tpm`; `mean_ntpm`; `sem_ntpm` (mean and SEM of normalized TPM) |
| `<name>_<tag>_persample_tpm.tsv` | `gene`; `id`; the unit column; `sampleid`; `stage`; `tpm` of every library, all units together; `tpm_sizefactor`; `ntpm` (`tpm / tpm_sizefactor`) |
| `<unit>_<tag>_tpm_sizefactors.tsv` | `sampleid`; `stage`; `tpm_sizefactor`: one row per library of the unit |
| `<name>_<tag>_panel.tsv` | only with `panel_aggregate: true`: `gene_id`; `rep_ids` (the transcripts summed, joined by `;`); `n_reps` |
| `<unit>_deseq2_LRT_all.tsv` | every feature of the transcriptome-wide test: `feature`; `baseMean` (mean normalized count); `log2FoldChange` and `lfcSE` (the last stage against the first, a by-product of the fit that the test does not depend on); `stat` (likelihood-ratio statistic); `pvalue`; `padj` (transcriptome-wide Benjamini-Hochberg) |
| `<unit>_dge_<A>_vs_<B>.tsv` | one file per pairwise comparison, `A` relative to `B`, every feature: `feature`; `baseMean`; `log2FoldChange` (shrunken); `lfcSE`; `stat` (Wald statistic); `pvalue`; `padj`; `gene` and `family` (panel features only); `is_rtk` (`TRUE` for a panel feature; the column name is fixed whatever the panel is); `rtk_sig` (`TRUE` for a panel feature with `padj` below `alpha` in this comparison; these are the starred genes of the supplement); the unit column; `contrast` |
| `<unit>_pca.tsv` | `sample`; `PC1`; `PC2`; `stage`: the first two principal components of the 500 most variable features of the variance-stabilized counts |
| `<unit>_pca_var.tsv` | `PC`; `pct`: the percentage of variance carried by each component |
| `<unit>_dds.rds` | the fitted DESeq2 object, the variance-stabilized counts and the sample table, for `supp_contrasts.R` |

The Python scripts write:

| File | Content |
| :--- | :--- |
| `<figure_dir>/<style_dir>/<filename>.png` and `.pdf` | main figures, one per family and unit: mean and standard error of normalized TPM per stage (bars or a line), the letters over regulated genes, one dot per animal in the enhanced style. PNG at 300 dpi; PDF with vector graphics and embedded fonts |
| `<figure_dir>/supplementary/<filename>.png` and `.pdf` (default `Fig<tag>_<unit>_dge`) | supplementary figures: one row per comparison, volcano plot left, MA plot right; panel genes in red, significant ones starred and named. `<name>_stars.json` next to them lists the starred (comparison, gene) pairs |
| `<table_dir>/<name>.tsv` and `.md` | the significance table: `gene`, `family`, the unit column, `LRT_p`, `q_targeted`, `q_transcriptomewide`, `molt_regulated`; the Markdown version formats the numbers in scientific notation |

## Checking a run

`validate_outputs.R` runs at the end of every `run_pipeline.sh` and exits non-zero on any failed
check; `run_pipeline.sh` then stops without printing `PIPELINE DONE`. The header of that script
lists the checks.

`compare_outputs.py` checks a rerun against an existing set of tables:

    python3 compare_outputs.py <reference_dir> <new_dir> --units U1,U2 --name <name> --tag <panel_tag>

It compares the v1.0.4 columns of the expected tables row by row; it does not compare the
normalized-TPM columns or the size-factor file. It compares text columns exactly, numeric columns
within a relative tolerance (default 1e-6) that absorbs floating-point noise between machines, and
p-values also within 0.01 on the log10 scale. It prints one line per table (PASS, FAIL with the
offending cells, SKIPPED, EXTRA, N/A) and exits with status 0 only when every expected table
passed. The docstring at the top of the script has the details.

## Reproducibility notes

- Results are deterministic: the same inputs and software give the same tables. The only
  randomness anywhere is the jitter of the dots in the enhanced figures, which uses a fixed seed.
- OpenBLAS builds that select a CPU-specific kernel at run time round differently on different
  processor generations, which moves the last decimal places of the transcriptome-wide estimates. Setting
  `OPENBLAS_CORETYPE` (for example `HASWELL`) in the site file pins one kernel and makes the tables
  reproducible to the last digit on any x86 machine; `DGE_BLAS_PIN=none` opts out on purpose. The R
  scripts refuse to run when neither is set; `run_pipeline.sh` sets `DGE_BLAS_PIN=none` itself, and
  says so, when there is no site file.
- `env_r.yaml` pins every package to the versions the pipeline was verified with.

## Adding a dataset

1. Write the sample sheet (one row per library) with the columns your configuration will name in
   `sample_id_col`, `quant_subdir_col`, `level_col`, and `unit_col`.
2. Write the panel file with `gene`, `transcript_id`, and `family`.
3. Copy `configs/example.yaml`, point `name`, `base`, `quant_dir`, `sample_sheet`, `panel`, and
   `output_dir` at the new dataset, and set `panel_tag` if you want a word other than `panel` in the
   table names.
4. For each unit, set `levels` and `outside_stage_test`, plus `exclude_levels`, `expected_n`, and
   the figure keys as needed.
5. Run the pipeline, confirm `VALIDATION PASSED` and `PIPELINE DONE`, then build the tables and
   figures with `table_panel.py`, `fig_main.py`, and `fig_supp.py` on the output directory.
