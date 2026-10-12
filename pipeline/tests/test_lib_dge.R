#!/usr/bin/env Rscript
# Unit tests for lib_dge.R. The configuration checks must reject every malformed configuration
# with a message that names the problem, and the model builder and the letter builder must return
# what the statistics require. Each test either writes a configuration that must be rejected and
# checks the error text, or builds a small example and checks the result. Run it from the package
# directory with `Rscript tests/test_lib_dge.R`; it prints one line per passed group and stops at
# the first failure.
suppressMessages({library(data.table); library(yaml); library(multcompView)})
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
cat("ALL lib_dge TESTS PASSED\n")
