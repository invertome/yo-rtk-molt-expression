# Shared helper functions for the DGE panel pipeline.
#
# The pipeline asks one question for each gene in a small panel of genes of interest: does its
# expression change across the stages of an experiment (for example the stages of a molt cycle, or
# the days after a treatment)? It answers it with DESeq2, fitted transcriptome-wide, and reports the panel
# genes in detail. The functions below are used by dge_panel.R, supp_contrasts.R and
# validate_outputs.R. None of them reads or writes files; they only check the configuration, build
# the statistical models, compute TPM size factors, and turn p-values into the summary columns of
# the output tables.
#
# Vocabulary used in the comments:
#   unit      one analysis run: a tissue or an experiment with its own set of samples
#   level     one value of the stage column of the sample sheet (a molt stage, a day, a condition)
#   panel     the genes of interest, listed in the panel TSV
#   feature   one row of the quantification (a transcript, or a gene after summing transcripts)
suppressMessages({library(data.table)})

# Estimate one TPM size factor for each library from all transcripts in the unit.
# DESeq2 compares each transcript with its geometric mean across libraries, then takes the median
# of those ratios for each library. Keep the library names so callers can match factors to samples
# by name. Missing or invalid values stop the run instead of producing misleading normalized TPM.
tpm_size_factors <- function(abundance) {
  libraries <- colnames(abundance)
  if (is.null(libraries)) stop("TPM matrix has no column names for libraries (columns ",
                                paste(seq_len(ncol(abundance)), collapse=","), ")")
  bad <- which(is.na(libraries) | !nzchar(libraries))
  if (length(bad)) stop("TPM matrix has a missing column name for library column ", paste(bad, collapse=","))
  dup <- unique(libraries[duplicated(libraries)])
  if (length(dup)) stop("TPM matrix has duplicated library names: ", paste(dup, collapse=","))
  bad <- libraries[colSums(is.na(abundance)) > 0]
  if (length(bad)) stop("TPM abundance has NA values for library: ", paste(bad, collapse=","))
  bad <- libraries[colSums(abundance < 0) > 0]
  if (length(bad)) stop("TPM abundance has negative values for library: ", paste(bad, collapse=","))
  bad <- libraries[colSums(!is.finite(abundance)) > 0]
  if (length(bad)) stop("TPM abundance has non-finite values for library: ", paste(bad, collapse=","))
  bad <- libraries[colSums(abundance > 0) == 0]
  if (length(bad)) stop("TPM abundance is all zero for library: ", paste(bad, collapse=","))
  sf <- tryCatch(DESeq2::estimateSizeFactorsForMatrix(abundance), error=function(e)
    stop("TPM size factors for libraries ", paste(libraries, collapse=","), ": ", conditionMessage(e)))
  bad <- libraries[!is.finite(sf) | sf <= 0]
  if (length(bad)) stop("TPM size factor must be finite and positive for library: ", paste(bad, collapse=","))
  sf
}

# Read one YAML configuration file and check it before anything is computed.
# Every key the pipeline reads is listed in configs/example.yaml. A missing key, an unknown key
# (usually a typo), or a value the pipeline cannot use stops the run here with a message that names
# the problem, so that a mistake in the configuration never turns into a silently wrong result.
read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  req <- c("name","base","quant_dir","sample_sheet","sample_id_col","quant_subdir_col","level_col","unit_col",
           "panel","panel_aggregate","output_dir","alpha","pairwise_adjust","shrinkage","units")
  miss <- setdiff(req, names(cfg))
  if (length(miss)) stop("config is missing: ", paste(miss, collapse=", "))
  if (!cfg$pairwise_adjust %in% c("holm","BH")) stop("pairwise_adjust must be holm or BH, got ", cfg$pairwise_adjust)
  if (!identical(cfg$shrinkage, "normal")) stop("shrinkage must be normal (the only shrinkage method the environment installs), got ", cfg$shrinkage)
  if (!length(cfg$units)) stop("config has no units")
  # panel_tag is the word used in the names of the panel tables (<unit>_<tag>_letters.tsv and
  # friends). It defaults to "panel"; a project may choose another word (for example kin).
  cfg$panel_tag <- if (is.null(cfg$panel_tag)) "panel" else as.character(cfg$panel_tag)
  if (length(cfg$panel_tag) != 1L || !grepl("^[A-Za-z0-9_]+$", cfg$panel_tag))
    stop("panel_tag must be one word of letters, digits or underscores, got '", paste(cfg$panel_tag, collapse=","), "'")
  # Every top-level key must be one the pipeline reads. A misspelled key (say "shrinkge") would
  # otherwise be ignored and the default used without anyone noticing, so it is an error instead.
  known <- c(req, "panel_tag", "figures", "tables")
  extra <- setdiff(names(cfg), known)
  if (length(extra)) stop("config has unknown top-level key(s): ", paste(extra, collapse=", "))
  unit_keys <- c("levels","outside_stage_test","closing_contrast","extra_contrasts","exclude_levels","expected_n",
                 "x_axis","x_label","numeric_from")
  # The same checks for each unit: the levels must be named, spelled with plain characters (the
  # pairwise results are keyed as "A-B", so a dash inside a level name would be ambiguous), and the
  # optional settings must refer to levels that exist.
  for (u in names(cfg$units)) {
    unit <- cfg$units[[u]]
    extra <- setdiff(names(unit), unit_keys)
    if (length(extra)) stop("unit ", u, ": unknown key(s): ", paste(extra, collapse=", "))
    if (length(unit$levels) < 2) stop("unit ", u, ": levels must list at least 2 levels")
    unit$levels <- as.character(unit$levels)
    bad <- unit$levels[!grepl("^[A-Za-z0-9_.]+$", unit$levels)]
    if (length(bad)) stop("unit ", u, ": level names may use only letters, digits, underscore or dot (pairwise letters are keyed as A-B): ", paste(bad, collapse=","))
    unit$outside_stage_test <- as.character(unlist(unit$outside_stage_test))
    if (is.null(unit$closing_contrast)) unit$closing_contrast <- FALSE
    unit$extra_contrasts <- lapply(unit$extra_contrasts, as.character)
    if (is.null(unit$x_axis)) unit$x_axis <- "categorical"
    bad <- setdiff(unit$outside_stage_test, unit$levels)
    if (length(bad)) stop("unit ", u, ": outside_stage_test level(s) not in levels: ", paste(bad, collapse=","))
    if (length(unit$levels) - length(unit$outside_stage_test) < 2) stop("unit ", u, ": fewer than 2 levels remain in the stage test")
    for (ct in unit$extra_contrasts) if (length(ct) != 2 || !all(ct %in% unit$levels))
      stop("unit ", u, ": extra contrast must be two known levels, got ", paste(ct, collapse=","))
    # Two optional safety settings. exclude_levels lists values of the stage column that are known to
    # be in the sample sheet but deliberately left out of the analysis (for example a sex or a
    # condition outside the design); any other unexpected value stops the run. expected_n gives the
    # number of samples expected at each level; dge_panel.R compares it with what the sheet actually
    # provides, so a wrong or outdated sheet stops the run before any model is fitted.
    unit$exclude_levels <- as.character(unlist(unit$exclude_levels))
    both <- intersect(unit$exclude_levels, unit$levels)
    if (length(both)) stop("unit ", u, ": exclude_levels overlap levels: ", paste(both, collapse=","))
    if (!is.null(unit$expected_n)) {
      en <- unit$expected_n
      miss_n <- setdiff(unit$levels, names(en)); extra_n <- setdiff(names(en), unit$levels)
      if (length(miss_n)) stop("unit ", u, ": expected_n lacks level(s): ", paste(miss_n, collapse=","))
      if (length(extra_n)) stop("unit ", u, ": expected_n names level(s) not in levels: ", paste(extra_n, collapse=","))
      ok <- vapply(en[unit$levels], function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v) && v >= 1,
                   logical(1))
      if (!all(ok)) stop("unit ", u, ": expected_n must be positive integers, got ",
                         paste(sprintf("%s=%s", unit$levels[!ok], sapply(en[unit$levels][!ok], format)), collapse=","))
      unit$expected_n <- setNames(as.integer(unlist(en[unit$levels])), unit$levels)
    }
    cfg$units[[u]] <- unit
  }
  cfg
}

# Build the two models that the likelihood-ratio test compares.
#
# The full model gives every level its own mean expression. The reduced model gives one shared
# mean to all the levels that belong to the stage test, and keeps a separate mean for each level
# listed in `outside` (for example animals that were blocked from progressing through the cycle:
# they stay in the fit, so they contribute to the dispersion estimates, but the test does not ask
# whether they differ from the cycle). The likelihood-ratio test then asks whether letting the
# stage levels have different means explains the counts better than one shared mean does. When
# `outside` is empty the reduced model is the usual intercept-only model.
#
# Both models are returned as design matrices, which DESeq2 accepts directly; the degrees of
# freedom of the test is the difference in the number of columns. The last check confirms that the
# reduced model is genuinely a special case of the full one (otherwise the test would be invalid).
reduced_matrix <- function(cd, level_col, outside) {
  lev <- cd[[level_col]]
  if (!is.factor(lev)) stop(level_col, " must be a factor")
  bad <- setdiff(outside, levels(lev))
  if (length(bad)) stop("outside level(s) not in the factor: ", paste(bad, collapse=","))
  full <- model.matrix(as.formula(paste0("~", level_col)), data=cd)
  if (!length(outside)) return(list(full=full, reduced=model.matrix(~1, data=cd), df=ncol(full) - 1L))
  cd$outside_ <- factor(ifelse(as.character(lev) %in% outside, as.character(lev), "cycle"), levels=c("cycle", outside))
  red <- model.matrix(~outside_, data=cd)
  proj <- full %*% solve(crossprod(full), crossprod(full, red))
  if (max(abs(proj - red)) > 1e-10) stop("reduced model is not nested in the full model")
  list(full=full, reduced=red, df=ncol(full) - ncol(red))
}

# Summarize the likelihood-ratio test for the panel features.
#
# `res` holds the transcriptome-wide DESeq2 result; `ids` are the panel features. For each panel feature
# the function reports the raw p-value of the test, the false discovery rate computed across the
# panel only ("targeted", Benjamini-Hochberg over a handful of genes, which is the number that
# decides whether a gene is called regulated), and the transcriptome-wide false discovery rate from
# DESeq2 for reference. A feature that DESeq2 dropped (too few counts, or filtered) gets NA in all
# three and is not called regulated.
targeted_stats <- function(res, ids, alpha) {
  p <- setNames(rep(NA_real_, length(ids)), ids); qg <- p
  inr <- ids[ids %in% rownames(res)]
  p[inr] <- res[inr, "pvalue"]; qg[inr] <- res[inr, "padj"]
  qt <- p.adjust(p, method="BH")
  data.table(id=ids, lrt_p=unname(p), q_targeted=unname(qt), q_transcriptomewide=unname(qg),
             molt_regulated=!is.na(qt) & qt < alpha)
}

# Turn one gene's pairwise comparisons into compact letters.
#
# `praw` holds the raw p-values of every pairwise Wald test between two stages, named "A-B". They
# are adjusted for multiple testing within the gene (Holm by default) and then summarized as the
# familiar letter groupings: two stages that share a letter are not significantly different, two
# stages with no letter in common are. The letters are returned in the order of `lv`, which is the
# order the stages are drawn in the figures.
letters_for <- function(praw, lv, method, alpha, label) {
  padj <- p.adjust(praw, method=method)
  if (any(is.na(padj))) stop(label, ": NA pairwise p-values for ", paste(names(padj)[is.na(padj)], collapse=","),
                             "; letters cannot be built")
  multcompView::multcompLetters(padj, threshold=alpha)$Letters[lv]
}

# List the pairwise comparisons that go into the supplementary tables and figures.
#
# "Transitions" are the consecutive steps through the stage test, each later stage against the
# earlier one, plus the closing step from the last stage back to the first when the configuration
# asks for it (a cycle). "Extra" comparisons are taken from the configuration as written, for
# example each stage against the animals kept outside the stage test.
contrasts_for_unit <- function(unit) {
  cycle <- setdiff(as.character(unit$levels), as.character(unlist(unit$outside_stage_test)))
  trans <- lapply(seq_len(length(cycle) - 1L), function(i) c(cycle[i + 1L], cycle[i]))
  if (isTRUE(unit$closing_contrast)) trans <- c(trans, list(c(cycle[1L], cycle[length(cycle)])))
  list(transitions=trans, extra=lapply(unit$extra_contrasts, as.character))
}

# Make the numerical setup explicit before any model is fitted.
#
# R does its linear algebra through a BLAS library. Some builds of OpenBLAS pick a CPU-specific
# code path at run time, and different paths round in slightly different ways, so the same data
# can give estimates that differ in the last decimal places from one machine to another. Setting
# OPENBLAS_CORETYPE pins one code path and makes a run reproducible everywhere it can run;
# DGE_BLAS_PIN=none says explicitly that this does not matter for the current run. With neither
# set the run stops here rather than produce numbers that cannot be reproduced exactly.
assert_blas_pin <- function() {
  ct <- Sys.getenv("OPENBLAS_CORETYPE"); opt <- Sys.getenv("DGE_BLAS_PIN")
  cat(sprintf("BLAS=%s | OPENBLAS_CORETYPE=%s | DGE_BLAS_PIN=%s\n", basename(extSoftVersion()[["BLAS"]]),
              if (nzchar(ct)) ct else "(unset)", if (nzchar(opt)) opt else "(unset)"))
  if (!nzchar(ct) && opt != "none")
    stop("set OPENBLAS_CORETYPE to pin one BLAS kernel (reproduces a run across machines) or DGE_BLAS_PIN=none to run unpinned")
  invisible(if (nzchar(ct)) ct else "none")
}
