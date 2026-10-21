#!/usr/bin/env Rscript
# Check the output tables of one run for internal consistency.
#
# Usage: Rscript validate_outputs.R <config.yaml> [output directory, if not the one in the config]
#
# This script reads the letters and TPM tables that dge_panel.R wrote and confirms that they say
# what they should: every gene appears once at every stage, the sample counts agree with the
# sample sheet, a gene is marked regulated exactly when its targeted false discovery rate is below
# alpha, only regulated genes carry more than one letter, and every gene without a test result
# has an identifiable reason (no counts). Normalized TPM must agree with plain TPM divided by
# the library's factor in the size-factor file. Stage means and SEM must agree with sample values,
# and SEM is NA only for a stage with one library.
# Any failed check prints a FAIL line and makes the script exit with status 1, which makes
# run_pipeline.sh stop before declaring the run done. WARN lines report oddities worth a look
# that do not invalidate the run.
suppressMessages({library(DESeq2); library(data.table); library(yaml)})
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) %in% 1:2) stop("usage: Rscript validate_outputs.R <config.yaml> [output_dir]")
PIPE <- Sys.getenv("DGE_PIPELINE_DIR")   # package root (falls back to this script's own directory)
if (!nzchar(PIPE)) { f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); PIPE <- if (length(f)) dirname(normalizePath(f[1])) else "." }
source(file.path(PIPE, "lib_dge.R"))
cfg <- read_config(args[1]); BASE <- Sys.getenv("BASE", cfg$base)
ODIR <- if (length(args) == 2) args[2] else file.path(BASE, cfg$output_dir)
ss <- fread(file.path(BASE, cfg$sample_sheet)); setnames(ss, tolower(names(ss)))
panel <- fread(file.path(BASE, cfg$panel)); genes <- unique(panel$gene)
fails <- character(0); warns <- character(0)
FAIL <- function(...) fails <<- c(fails, sprintf(...)); WARN <- function(...) warns <<- c(warns, sprintf(...))

psf <- file.path(ODIR, paste0(cfg$name, "_", cfg$panel_tag, "_persample_tpm.tsv"))
ps <- NULL
if (!file.exists(psf)) {
  FAIL("%s: missing", basename(psf))
} else {
  ps <- fread(psf, colClasses=list(character="sampleid"))
  missing <- setdiff(c("gene", cfg$unit_col, "sampleid", "stage", "tpm", "tpm_sizefactor", "ntpm"), names(ps))
  if (length(missing)) {
    FAIL("%s: missing column(s): %s", basename(psf), paste(missing, collapse=","))
    ps <- NULL
  } else {
    ps[, stage := as.character(stage)]
    extra <- setdiff(unique(ps[[cfg$unit_col]]), names(cfg$units))
    if (length(extra)) FAIL("%s: unexpected unit(s): %s", basename(psf), paste(extra, collapse=","))
  }
}

for (uname in names(cfg$units)) {
  unit <- cfg$units[[uname]]; lv <- unit$levels
  lef <- file.path(ODIR, paste0(uname, "_", cfg$panel_tag, "_letters.tsv")); stf <- file.path(ODIR, paste0(uname, "_", cfg$panel_tag, "_tpm_by_stage.tsv"))
  sff <- file.path(ODIR, paste0(uname, "_", cfg$panel_tag, "_tpm_sizefactors.tsv"))
  for (f in c(lef, stf, sff)) if (!file.exists(f)) { FAIL("%s: missing", basename(f)); next }
  if (!all(file.exists(c(lef, stf, sff)))) next
  le <- fread(lef); st <- fread(stf)
  le[, stage := as.character(stage)]; st[, stage := as.character(stage)]
  sft <- fread(sff, colClasses=list(character="sampleid"))
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
  sample_ids <- as.character(meta[[cfg$sample_id_col]])
  n_expected <- as.integer(table(factor(meta[[cfg$level_col]], levels=lv)))
  for (i in seq_along(lv)) {
    n_have <- unique(st[stage == lv[i], n])
    if (length(n_have) != 1 || n_have != n_expected[i]) FAIL("%s tpm: n at level %s is %s, sample sheet says %d", uname, lv[i], paste(n_have, collapse="/"), n_expected[i])
  }
  # The factor file must hold one finite, positive factor for every library in the unit.
  sf_ok <- TRUE
  missing <- setdiff(c("sampleid", "stage", "tpm_sizefactor"), names(sft))
  if (length(missing)) {
    FAIL("%s: missing column(s): %s", basename(sff), paste(missing, collapse=",")); sf_ok <- FALSE
  } else {
    sft[, stage := as.character(stage)]
    have <- sft$sampleid
    if (!setequal(have, sample_ids) || anyDuplicated(have)) {
      FAIL("%s: sample ids must match the sample sheet once each; missing: %s; unexpected: %s; duplicated: %s",
           basename(sff), paste(setdiff(sample_ids, have), collapse=","), paste(setdiff(have, sample_ids), collapse=","),
           paste(unique(have[duplicated(have)]), collapse=","))
      sf_ok <- FALSE
    }
    if (!is.numeric(sft$tpm_sizefactor) || any(!is.finite(sft$tpm_sizefactor) | sft$tpm_sizefactor <= 0)) {
      FAIL("%s: tpm_sizefactor must be finite and positive, with no NA values", basename(sff)); sf_ok <- FALSE
    }
  }
  # Normalized summaries must be present, finite and nonnegative, except for the undefined SEM
  # of a stage with one library. Means and SEM must describe the samples of that gene and stage.
  # The tolerance is relative to the expected value. An expected zero must match exactly.
  st_ok <- TRUE
  for (col in c("mean_ntpm", "sem_ntpm")) {
    if (!col %in% names(st)) {
      FAIL("%s tpm: missing column %s", uname, col); st_ok <- FALSE
    } else {
      values <- st[[col]]
      # fread can infer a logical column when every value is NA, as with all-singleton units.
      if (!is.numeric(values) && !(is.logical(values) && all(is.na(values)))) {
        FAIL("%s tpm: %s must be numeric", uname, col); st_ok <- FALSE; next
      }
      allowed_na <- col == "sem_ntpm" & !is.na(st$n) & st$n == 1 & is.na(values) & !is.nan(values)
      if (any((!is.finite(values) | values < 0) & !allowed_na)) {
        FAIL("%s tpm: %s must be finite and nonnegative%s", uname, col,
             if (col == "sem_ntpm") "; NA is allowed only when n == 1" else ", with no NA values")
        st_ok <- FALSE
      }
    }
  }
  if (!is.null(ps)) {
    pu <- ps[get(cfg$unit_col) == uname]
    for (gn in genes) {
      have <- pu[gene == gn, sampleid]
      if (!setequal(have, sample_ids) || anyDuplicated(have))
        FAIL("%s persample: %s must have one row per library; missing: %s; unexpected: %s; duplicated: %s",
             uname, gn, paste(setdiff(sample_ids, have), collapse=","), paste(setdiff(have, sample_ids), collapse=","),
             paste(unique(have[duplicated(have)]), collapse=","))
    }
    extra <- setdiff(unique(pu$gene), genes)
    if (length(extra)) FAIL("%s persample: unexpected gene(s): %s", uname, paste(extra, collapse=","))
    stages <- setNames(as.character(meta[[cfg$level_col]]), sample_ids)[pu$sampleid]
    bad <- which(is.na(pu$stage) | is.na(stages) | as.character(pu$stage) != stages)
    if (length(bad)) FAIL("%s persample: stage does not match sample sheet for library: %s", uname,
                          paste(unique(pu$sampleid[bad]), collapse=","))
    ps_ok <- TRUE
    for (col in c("tpm", "tpm_sizefactor", "ntpm")) {
      values <- pu[[col]]
      if (!is.numeric(values) || any(!is.finite(values) | values < 0) ||
          (col == "tpm_sizefactor" && any(values == 0))) {
        FAIL("%s persample: %s must be finite and %s, with no NA values", uname, col,
             if (col == "tpm_sizefactor") "positive" else "nonnegative")
        ps_ok <- FALSE
      }
    }
    if (ps_ok) {
      if (sf_ok) {
        expected_sf <- setNames(sft$tpm_sizefactor, sft$sampleid)[pu$sampleid]
        bad <- which(!is.finite(expected_sf) | abs(pu$tpm_sizefactor - expected_sf) > 1e-9 * abs(expected_sf))
        if (length(bad)) FAIL("%s persample: tpm_sizefactor differs from %s for library: %s", uname,
                              basename(sff), paste(unique(pu$sampleid[bad]), collapse=","))
      }
      expected <- pu$tpm / pu$tpm_sizefactor
      bad <- which(!is.finite(expected) | abs(pu$ntpm - expected) > 1e-9 * abs(expected))
      if (length(bad)) FAIL("%s persample: ntpm differs from tpm / tpm_sizefactor for gene/library: %s", uname,
                            paste(paste(pu$gene[bad], pu$sampleid[bad], sep="/"), collapse=","))
      if (st_ok) {
        summaries <- pu[, .(expected_ntpm=mean(ntpm), expected_sem_ntpm=sd(ntpm)/sqrt(.N)), by=.(gene, stage)]
        check <- merge(st[, .(gene, stage, mean_ntpm, sem_ntpm)], summaries, by=c("gene", "stage"), all=TRUE)
        bad <- which(!is.finite(check$mean_ntpm) | !is.finite(check$expected_ntpm) |
                     abs(check$mean_ntpm - check$expected_ntpm) > 1e-9 * abs(check$expected_ntpm))
        if (length(bad)) FAIL("%s tpm: mean_ntpm differs from persample mean for gene/stage: %s", uname,
                              paste(paste(check$gene[bad], check$stage[bad], sep="/"), collapse=","))
        bad <- which(is.na(check$sem_ntpm) != is.na(check$expected_sem_ntpm) |
                     (!is.na(check$expected_sem_ntpm) &
                      (!is.finite(check$sem_ntpm) | !is.finite(check$expected_sem_ntpm) |
                       abs(check$sem_ntpm - check$expected_sem_ntpm) > 1e-9 * abs(check$expected_sem_ntpm))))
        if (length(bad)) FAIL("%s tpm: sem_ntpm differs from persample SEM for gene/stage: %s", uname,
                              paste(paste(check$gene[bad], check$stage[bad], sep="/"), collapse=","))
      }
    }
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
