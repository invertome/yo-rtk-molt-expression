#!/usr/bin/env Rscript
# Check the output tables of one run for internal consistency.
#
# Usage: Rscript validate_outputs.R <config.yaml> [output directory, if not the one in the config]
#
# This script does not recompute anything. It reads the letters and TPM tables that dge_panel.R
# wrote and confirms that they say what they should: every gene appears once at every stage, the
# sample counts agree with the sample sheet, a gene is marked regulated exactly when its targeted
# false discovery rate is below alpha, only regulated genes carry more than one letter, and every
# gene without a test result has an identifiable reason (no counts). Any failed check prints a FAIL
# line and makes the script exit with status 1, which makes run_pipeline.sh stop before declaring
# the run done. WARN lines report oddities worth a look that do not invalidate the run.
suppressMessages({library(DESeq2); library(data.table); library(yaml)})
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) %in% 1:2) stop("usage: Rscript validate_outputs.R <config.yaml> [output_dir]")
PIPE <- Sys.getenv("DGE_PIPELINE_DIR")   # package root; falls back to this script's own directory
if (!nzchar(PIPE)) { f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); PIPE <- if (length(f)) dirname(normalizePath(f[1])) else "." }
source(file.path(PIPE, "lib_dge.R"))
cfg <- read_config(args[1]); BASE <- Sys.getenv("BASE", cfg$base)
ODIR <- if (length(args) == 2) args[2] else file.path(BASE, cfg$output_dir)
ss <- fread(file.path(BASE, cfg$sample_sheet)); setnames(ss, tolower(names(ss)))
panel <- fread(file.path(BASE, cfg$panel)); genes <- unique(panel$gene)
fails <- character(0); warns <- character(0)
FAIL <- function(...) fails <<- c(fails, sprintf(...)); WARN <- function(...) warns <<- c(warns, sprintf(...))

for (uname in names(cfg$units)) {
  unit <- cfg$units[[uname]]; lv <- unit$levels
  lef <- file.path(ODIR, paste0(uname, "_", cfg$panel_tag, "_letters.tsv")); stf <- file.path(ODIR, paste0(uname, "_", cfg$panel_tag, "_tpm_by_stage.tsv"))
  for (f in c(lef, stf)) if (!file.exists(f)) { FAIL("%s: missing", basename(f)); next }
  if (!all(file.exists(c(lef, stf)))) next
  le <- fread(lef); st <- fread(stf)
  # 1. Both tables hold exactly one row per gene and stage, no more and no fewer.
  for (nm in c("letters", "tpm")) {
    d <- if (nm == "letters") le else st
    have <- d[, paste(gene, stage)]; want <- as.vector(outer(genes, lv, paste))
    for (m in setdiff(want, have)) FAIL("%s %s: missing row gene x level '%s'", uname, nm, m)
    for (x in setdiff(have, want)) FAIL("%s %s: unexpected row '%s'", uname, nm, x)
    if (anyDuplicated(have)) FAIL("%s %s: duplicated rows: %s", uname, nm, paste(have[duplicated(have)], collapse=","))
  }
  # The sample count reported at each stage equals the number of rows the sample sheet has for it.
  meta <- ss[get(cfg$unit_col) == uname & get(cfg$level_col) %in% lv]
  n_expected <- as.integer(table(factor(meta[[cfg$level_col]], levels=lv)))
  for (i in seq_along(lv)) {
    n_have <- unique(st[stage == lv[i], n])
    if (length(n_have) != 1 || n_have != n_expected[i]) FAIL("%s tpm: n at level %s is %s, sample sheet says %d", uname, lv[i], paste(n_have, collapse="/"), n_expected[i])
  }
  # 2. The regulated call is exactly "targeted FDR below alpha", and a gene that is not regulated
  #    carries the same letter everywhere. A regulated gene whose stages all share one letter is
  #    possible (the overall test is significant but no single pair survives the within-gene
  #    adjustment); it is reported as a warning so that a reader can look at it.
  g <- le[, .(q=q_targeted[1], reg=molt_regulated[1], p=lrt_p[1], letters=paste(letter, collapse="")), by=gene]
  for (i in seq_len(nrow(g))) {
    exp_reg <- !is.na(g$q[i]) && g$q[i] < cfg$alpha
    if (g$reg[i] != exp_reg) FAIL("%s letters: %s molt_regulated=%s but q_targeted=%s (alpha %g)", uname, g$gene[i], g$reg[i], g$q[i], cfg$alpha)
    lets <- le[gene == g$gene[i], letter]
    if (!g$reg[i] && any(lets != "a")) FAIL("%s letters: %s is not molt-regulated but carries letters %s", uname, g$gene[i], paste(lets, collapse=","))
    if (g$reg[i] && length(unique(lets)) == 1) WARN("%s letters: %s is molt-regulated (q=%.3g) but no pair of levels differs (all '%s')", uname, g$gene[i], g$q[i], lets[1])
  }
  # 3. A gene with no test result (NA p-value) must have a reason. With Cook's-distance exclusion
  #    switched off, the only legitimate reason is that the gene had no counts at all; anything else
  #    is a defect and fails the run.
  na_genes <- g[is.na(p), gene]
  cat(sprintf("[%s] targeted-BH denominator: %d of %d genes testable\n", uname, sum(!is.na(g$p)), nrow(g)))
  if (length(na_genes)) {
    obj <- readRDS(file.path(ODIR, paste0(uname, "_dds.rds"))); dds <- obj$dds
    feat_of <- if (isTRUE(cfg$panel_aggregate)) setNames(genes, genes) else setNames(panel$transcript_id, panel$gene)
    for (gn in na_genes) {
      ft <- feat_of[[gn]]
      if (!ft %in% rownames(dds)) { FAIL("%s: %s has NA LRT p and is absent from the fitted object", uname, gn); next }
      bm <- mcols(dds)[ft, "baseMean"]; mc <- mcols(dds)[ft, "maxCooks"]   # reported only
      reason <- if (is.na(bm) || bm == 0) "no counts" else "unexplained"   # Cook's exclusion is off (cooksCutoff=FALSE), so an NA with counts is a defect
      if (reason == "unexplained") FAIL("%s: %s has NA LRT p with no identifiable reason (baseMean %.3g, maxCooks %s)", uname, gn, bm, mc)
      else cat(sprintf("[%s] NA LRT p: %s -> %s\n", uname, gn, reason))
    }
  }
}
for (w in warns) cat("WARN ", w, "\n"); for (f in fails) cat("FAIL ", f, "\n")
cat(sprintf("validation: %d FAIL, %d WARN\n", length(fails), length(warns)))
if (length(fails)) { cat("VALIDATION FAILED\n"); quit(status=1) } else cat("VALIDATION PASSED\n")
