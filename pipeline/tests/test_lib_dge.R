#!/usr/bin/env Rscript
# Tests for lib_dge.R and the output validator. The configuration checks reject malformed input
# with a message that names the problem, and the model builder and the letter builder must return
# what the statistics require. The size-factor group checks agreement with DESeq2 and rejection of
# invalid libraries. Output checks cover the factor file, normalized summaries and singleton stages.
# Each test builds a small example and checks its result or error text. Run it from the package
# directory with `Rscript tests/test_lib_dge.R`; it prints one line per passed group and stops at
# the first failure.
suppressMessages({library(data.table); library(yaml); library(multcompView); library(DESeq2)})
source(file.path(Sys.getenv("DGE_PIPELINE_DIR", "."), "lib_dge.R"))
expect_error <- function(expr, pattern) {
  msg <- tryCatch({ expr; "" }, error=function(e) conditionMessage(e))
  if (!grepl(pattern, msg)) stop(sprintf("expected error matching %s, got: %s", pattern, msg))
}
tmp <- tempfile(fileext=".yaml")
base_cfg <- list(name="t", base="/x", quant_dir="q", sample_sheet="s.tsv", sample_id_col="sampleid",
  quant_subdir_col="filename", level_col="stage", unit_col="experiment", panel="p.tsv", panel_aggregate=FALSE,
  output_dir="o", alpha=0.05, pairwise_adjust="holm", shrinkage="normal",
  units=list(U1=list(levels=c("A","B","C","D","E","X"), outside_stage_test="X", closing_contrast=TRUE,
                      extra_contrasts=list(c("X","A"), c("X","B"), c("X","C"), c("X","D"), c("X","E")), x_axis="categorical"),
             U2=list(levels=c("T0","T1","T3","T5","T7","T14","T19"), x_axis="numeric")))
write_yaml(base_cfg, tmp); cfg <- read_config(tmp)
stopifnot(identical(cfg$units$U2$outside_stage_test, character(0)), isFALSE(cfg$units$U2$closing_contrast))
bad <- base_cfg; bad$units <- NULL; write_yaml(bad, tmp); expect_error(read_config(tmp), "units")
bad <- base_cfg; bad$units$U1$outside_stage_test <- "XX"; write_yaml(bad, tmp); expect_error(read_config(tmp), "XX")
bad <- base_cfg; bad$shrinkage <- "apeglm"; write_yaml(bad, tmp); expect_error(read_config(tmp), "shrinkage")
bad <- base_cfg; bad$shrinkage <- "ashr"; write_yaml(bad, tmp); expect_error(read_config(tmp), "shrinkage")

cd <- data.frame(stage=factor(rep(c("A","B","C","D","E","X"), times=c(12,10,6,7,12,7)),
                              levels=c("A","B","C","D","E","X")))
mm <- reduced_matrix(cd, "stage", "X")
stopifnot(mm$df == 4, ncol(mm$full) == 6, ncol(mm$reduced) == 2, all(mm$reduced[cd$stage == "X", 2] == 1),
          all(mm$reduced[cd$stage != "X", 2] == 0))
mm0 <- reduced_matrix(cd, "stage", character(0)); stopifnot(mm0$df == 5, ncol(mm0$reduced) == 1)
expect_error(reduced_matrix(cd, "stage", "ZZ"), "ZZ")

res <- data.frame(pvalue=c(0.001, 0.02, 0.5, 0.9), padj=c(0.01, 0.2, 0.8, 0.95), row.names=c("g1","g2","g3","g4"))
ts <- targeted_stats(res, c("g1","g2","g3","g5"), 0.05)
stopifnot(nrow(ts) == 4, is.na(ts[id == "g5", lrt_p]),
          isTRUE(all.equal(ts[id != "g5", q_targeted], p.adjust(c(0.001, 0.02, 0.5), "BH"))),
          identical(ts$molt_regulated, c(TRUE, TRUE, FALSE, FALSE)),
          isTRUE(all.equal(ts$q_transcriptomewide, c(0.01, 0.2, 0.8, NA_real_))))

lv <- c("A","B","C"); praw <- c("A-B"=0.001, "A-C"=0.4, "B-C"=0.002)
let <- letters_for(praw, lv, "holm", 0.05, "test gene")
stopifnot(identical(names(let), lv), let[["A"]] == let[["C"]], let[["A"]] != let[["B"]])
expect_error(letters_for(c("A-B"=NA, "A-C"=0.4, "B-C"=0.002), lv, "holm", 0.05, "geneX"), "geneX.*A-B")

ct <- contrasts_for_unit(cfg$units$U1)
stopifnot(identical(ct$transitions, list(c("B","A"), c("C","B"), c("D","C"), c("E","D"), c("A","E"))),
          identical(ct$extra, list(c("X","A"), c("X","B"), c("X","C"), c("X","D"), c("X","E"))))
ct <- contrasts_for_unit(cfg$units$U2)
stopifnot(identical(ct$transitions, list(c("T1","T0"), c("T3","T1"), c("T5","T3"), c("T7","T5"), c("T14","T7"), c("T19","T14"))),
          length(ct$extra) == 0)
ct <- contrasts_for_unit(list(levels=c("A","C","D","E"), outside_stage_test=character(0), closing_contrast=TRUE, extra_contrasts=list()))
stopifnot(identical(ct$transitions, list(c("C","A"), c("D","C"), c("E","D"), c("A","E"))))
## optional per-unit keys exclude_levels / expected_n: accepted, normalized, or rejected with the offender named
full <- base_cfg; full$units$U2$exclude_levels <- "other"
full$units$U2$expected_n <- list(T19=5L, T0=13L, T1=8L, T3=9L, T5=7L, T7=8L, T14=9L)   # out of level order on purpose
write_yaml(full, tmp); cfgf <- read_config(tmp)
stopifnot(identical(cfgf$units$U2$exclude_levels, "other"),
          identical(cfgf$units$U2$expected_n, c(T0=13L, T1=8L, T3=9L, T5=7L, T7=8L, T14=9L, T19=5L)),
          identical(cfgf$units$U1$exclude_levels, character(0)), is.null(cfgf$units$U1$expected_n))
stopifnot(identical(cfg$units$U2$exclude_levels, character(0)), is.null(cfg$units$U2$expected_n))
bad <- full; bad$units$U2$exclude_levels <- c("other", "T3"); write_yaml(bad, tmp); expect_error(read_config(tmp), "U2: exclude_levels overlap levels: T3")
bad <- full; bad$units$U2$expected_n$T19 <- NULL; write_yaml(bad, tmp); expect_error(read_config(tmp), "U2: expected_n lacks level.*: T19")
bad <- full; bad$units$U2$expected_n$XX <- 4L; write_yaml(bad, tmp); expect_error(read_config(tmp), "U2: expected_n names level.* not in levels: XX")
bad <- full; bad$units$U2$expected_n$T1 <- 8.5; write_yaml(bad, tmp); expect_error(read_config(tmp), "U2: expected_n must be positive integers, got T1=8\\.5")
bad <- full; bad$units$U2$expected_n$T5 <- 0L; write_yaml(bad, tmp); expect_error(read_config(tmp), "U2: expected_n must be positive integers, got T5=0")

## unknown keys and bad level names are rejected, with the offender named
bad <- base_cfg; bad$outputdir <- "o"; write_yaml(bad, tmp); expect_error(read_config(tmp), "unknown top-level key.*outputdir")
bad <- base_cfg; bad$units$U1$expectd_n <- list(A=1L); write_yaml(bad, tmp); expect_error(read_config(tmp), "U1: unknown key.*expectd_n")
bad <- base_cfg; bad$units$U1$levels <- c("A-B","C","D","E","F","X"); write_yaml(bad, tmp); expect_error(read_config(tmp), "U1: level names.*A-B")
## panel_tag and the BLAS guard
stopifnot(identical(cfg$panel_tag, "panel"))
tg <- base_cfg; tg$panel_tag <- "kin"; write_yaml(tg, tmp); stopifnot(identical(read_config(tmp)$panel_tag, "kin"))
tg <- base_cfg; tg$panel_tag <- "a b"; write_yaml(tg, tmp); expect_error(read_config(tmp), "panel_tag")
Sys.setenv(OPENBLAS_CORETYPE="HASWELL"); Sys.unsetenv("DGE_BLAS_PIN"); stopifnot(identical(assert_blas_pin(), "HASWELL"))
Sys.unsetenv("OPENBLAS_CORETYPE"); Sys.setenv(DGE_BLAS_PIN="none"); stopifnot(identical(assert_blas_pin(), "none"))
Sys.unsetenv("DGE_BLAS_PIN"); expect_error(assert_blas_pin(), "DGE_BLAS_PIN=none")

## tsf_matches_deseq2: doubling most transcripts clearly changes the library's factor.
abundance <- cbind(library1=rep(10, 50), library2=101:150, library3=rep(20, 50), library4=rep(10, 50))
abundance[1:30, "library4"] <- abundance[1:30, "library1"] * 2
sf <- tpm_size_factors(abundance)
stopifnot(identical(sf, DESeq2::estimateSizeFactorsForMatrix(abundance)),
          identical(names(sf), colnames(abundance)), abs(sf[["library4"]] / sf[["library1"]] - 1) > 0.10)
reordered <- abundance[, c(4, 2, 1, 3)]
stopifnot(identical(tpm_size_factors(reordered), DESeq2::estimateSizeFactorsForMatrix(reordered)),
          identical(names(tpm_size_factors(reordered)), colnames(reordered)))
# A zero in one transcript is allowed when other transcripts support the estimate.
with_zero <- abundance; with_zero[1, "library2"] <- 0
stopifnot(identical(tpm_size_factors(with_zero), DESeq2::estimateSizeFactorsForMatrix(with_zero)))
cat("tsf_matches_deseq2 PASSED\n")

## tsf_rejects: errors identify the invalid value and the affected library.
bad <- abundance; colnames(bad) <- NULL
expect_error(tpm_size_factors(bad), "column names.*librar")
bad <- abundance; colnames(bad)[2] <- ""
expect_error(tpm_size_factors(bad), "column name.*library column 2")
bad <- abundance; colnames(bad)[2] <- colnames(bad)[1]
expect_error(tpm_size_factors(bad), "duplicated.*library1")
bad <- abundance; bad[1, "library2"] <- -1
expect_error(tpm_size_factors(bad), "negative.*library2")
bad <- abundance; bad[1, "library3"] <- NA_real_
expect_error(tpm_size_factors(bad), "NA.*library3")
bad <- abundance; bad[1, "library2"] <- Inf
expect_error(tpm_size_factors(bad), "non-finite.*library2")
bad <- abundance; bad[, "library4"] <- 0
expect_error(tpm_size_factors(bad), "zero.*library4")
# Finite input can still produce an infinite or zero factor through numerical overflow.
bad <- abundance; bad[,] <- 1e-300; bad[, "library4"] <- 1e300
expect_error(tpm_size_factors(bad), "size factor.*finite.*positive.*library4")
bad <- abundance; bad[,] <- 1e300; bad[, "library4"] <- 1e-300
expect_error(tpm_size_factors(bad), "size factor.*finite.*positive.*library4")
# Each library has expression, but no transcript is positive in every library.
bad <- abundance; bad[1:25, 1:2] <- 0; bad[26:50, 3:4] <- 0
expect_error(tpm_size_factors(bad), "libraries.*library1.*library4.*zero")
cat("tsf_rejects PASSED\n")

## Normalized output checks run the validator against small tables, without fitting any models.
# Numeric stage labels in U1 and text labels in U2 make fread infer different stage types in the
# per-unit and combined tables. U2 also tests a SEM column containing only singleton NA values.
output_dir <- tempfile("ntpm_outputs_"); dir.create(output_dir)
output_cfg <- base_cfg; output_cfg$base <- output_dir; output_cfg$output_dir <- "."
output_cfg$sample_sheet <- "samples.tsv"; output_cfg$panel <- "panel.tsv"
output_cfg$units <- list(U1=list(levels=c("0", "1")), U2=list(levels=c("A", "B")))
write_yaml(output_cfg, file.path(output_dir, "config.yaml"))
samples <- data.table(sampleid=paste0("library", 1:5), filename=paste0("quant", 1:5),
                      experiment=c("U1", "U1", "U1", "U2", "U2"), stage=c("0", "0", "1", "A", "B"))
fwrite(samples, file.path(output_dir, "samples.tsv"), sep="\t")
fwrite(data.table(gene="g1", transcript_id="tx1", family="f1"), file.path(output_dir, "panel.tsv"), sep="\t")
ps_base <- samples[, .(gene="g1", id="tx1", experiment, sampleid, stage)]
ps_base[, tpm_sizefactor := c(2, 0.25, 4, 1, 2)]
ps_base[, ntpm := c(2, 4, 3, 5, 3) * 1e-4]
ps_base[, tpm := ntpm * tpm_sizefactor]
st_base <- ps_base[, .(id=id[1], n=.N, mean_tpm=mean(tpm), sem_tpm=sd(tpm)/sqrt(.N),
                       mean_ntpm=mean(ntpm), sem_ntpm=sd(ntpm)/sqrt(.N)), by=.(gene, experiment, stage)]
sf_base <- ps_base[, .(sampleid, stage, tpm_sizefactor, experiment)]
ps_path <- file.path(output_dir, "t_panel_persample_tpm.tsv")
sf_path <- file.path(output_dir, "U1_panel_tpm_sizefactors.tsv")
st_path <- file.path(output_dir, "U1_panel_tpm_by_stage.tsv")
write_output_fixture <- function() {
  fwrite(ps_base, ps_path, sep="\t")
  for (u in names(output_cfg$units)) {
    fwrite(st_base[experiment == u], file.path(output_dir, paste0(u, "_panel_tpm_by_stage.tsv")), sep="\t")
    # Reverse the factor rows so a positional comparison cannot pass.
    sft <- sf_base[experiment == u, .(sampleid, stage, tpm_sizefactor)]
    fwrite(sft[rev(seq_len(nrow(sft)))], file.path(output_dir, paste0(u, "_panel_tpm_sizefactors.tsv")), sep="\t")
    le <- st_base[experiment == u, .(gene, experiment, stage, letter="a", lrt_p=0.5,
                                    q_targeted=0.5, q_transcriptomewide=0.5, molt_regulated=FALSE)]
    fwrite(le, file.path(output_dir, paste0(u, "_panel_letters.tsv")), sep="\t")
  }
}
expect_validation <- function(status, pattern) {
  script <- normalizePath(file.path(Sys.getenv("DGE_PIPELINE_DIR", "."), "validate_outputs.R"))
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  c("--vanilla", shQuote(script), shQuote(file.path(output_dir, "config.yaml")),
                                    shQuote(output_dir)), stdout=TRUE, stderr=TRUE,
                                  env=paste0("BASE=", shQuote(output_dir))))
  got <- attr(out, "status"); if (is.null(got)) got <- 0L
  if (got != status || !any(grepl(pattern, out)))
    stop(sprintf("expected validator status %d and %s, got %d:\n%s", status, pattern, got, paste(out, collapse="\n")))
}
write_output_fixture(); expect_validation(0L, "VALIDATION PASSED")
cat("ntpm_accepts_singletons_and_numeric_stages PASSED\n")

write_output_fixture(); unlink(sf_path)
expect_validation(1L, "FAIL .*tpm_sizefactors.*missing")
for (col in c("sampleid", "stage", "tpm_sizefactor")) {
  write_output_fixture(); bad <- fread(sf_path); bad[, (col) := NULL]
  fwrite(bad, sf_path, sep="\t")
  expect_validation(1L, "FAIL .*tpm_sizefactors.*missing column")
}
for (change in c("missing", "unexpected", "duplicate")) {
  write_output_fixture(); bad <- fread(sf_path)
  if (change == "missing") bad <- bad[-1]
  if (change == "unexpected") bad[1, sampleid := "unexpected_library"]
  if (change == "duplicate") bad <- rbind(bad, bad[1])
  fwrite(bad, sf_path, sep="\t")
  expect_validation(1L, "FAIL .*tpm_sizefactors.*sample")
}
for (value in c(NA_real_, Inf, -Inf, 0, -1)) {
  write_output_fixture(); bad <- fread(sf_path); bad[1, tpm_sizefactor := value]
  fwrite(bad, sf_path, sep="\t")
  expect_validation(1L, "FAIL .*tpm_sizefactors.*finite.*positive")
}
write_output_fixture(); bad <- fread(sf_path)
bad[sampleid == "library2", tpm_sizefactor := tpm_sizefactor * (1 + 2e-9)]
fwrite(bad, sf_path, sep="\t")
expect_validation(1L, "FAIL .*persample.*tpm_sizefactor.*library2")
write_output_fixture(); bad <- fread(sf_path)
bad[sampleid == "library2", tpm_sizefactor := tpm_sizefactor * (1 + 1e-10)]
fwrite(bad, sf_path, sep="\t")
expect_validation(0L, "VALIDATION PASSED")
cat("ntpm_checks_sizefactor_file PASSED\n")

write_output_fixture(); bad <- fread(st_path)
bad[stage == 0, sem_ntpm := sem_ntpm * (1 + 2e-9)]; fwrite(bad, st_path, sep="\t")
expect_validation(1L, "FAIL .*sem_ntpm.*persample")
write_output_fixture(); bad <- fread(st_path)
bad[stage == 0, sem_ntpm := sem_ntpm * (1 + 1e-10)]; fwrite(bad, st_path, sep="\t")
expect_validation(0L, "VALIDATION PASSED")
write_output_fixture(); bad <- fread(st_path)
bad[stage == 0, sem_ntpm := NA_real_]; fwrite(bad, st_path, sep="\t")
expect_validation(1L, "FAIL .*sem_ntpm")
for (value in c(0, Inf, -1)) {
  write_output_fixture(); bad <- fread(st_path); bad[stage == 1, sem_ntpm := value]
  fwrite(bad, st_path, sep="\t")
  expect_validation(1L, "FAIL .*sem_ntpm")
}
write_output_fixture()
zero_ps <- copy(ps_base); zero_ps[experiment == "U1" & stage == "0", `:=`(tpm=0, ntpm=0)]
zero_st <- copy(st_base[experiment == "U1"])
zero_st[stage == "0", `:=`(mean_tpm=0, sem_tpm=0, mean_ntpm=0, sem_ntpm=0)]
fwrite(zero_ps, ps_path, sep="\t"); fwrite(zero_st, st_path, sep="\t")
expect_validation(0L, "VALIDATION PASSED")
zero_st[stage == "0", sem_ntpm := 1e-20]; fwrite(zero_st, st_path, sep="\t")
expect_validation(1L, "FAIL .*sem_ntpm.*persample")
unlink(output_dir, recursive=TRUE)
cat("ntpm_checks_sem PASSED\n")
cat("ALL lib_dge TESTS PASSED\n")
