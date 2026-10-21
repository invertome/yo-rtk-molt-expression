#!/usr/bin/env Rscript
# Supplementary tables: pairwise comparisons between stages, and a principal component analysis.
#
# Usage: Rscript supp_contrasts.R <config.yaml>   (after dge_panel.R, which saves the fitted models)
#
# For every unit, and for every comparison listed by contrasts_for_unit (each stage against the
# previous one, the closing step of a cycle when configured, and the extra comparisons of the
# configuration), this script writes one transcriptome-wide table with the log2 fold change of every
# feature, its standard error, the Wald test p-value and the transcriptome-wide adjusted p-value, plus the
# panel annotation. Fold changes are shrunken (DESeq2's lfcShrink) so that genes with few counts do
# not show inflated changes on the volcano and MA plots of the supplement. It also writes the first
# two principal components of the variance-stabilized counts, one point per sample, for a quick
# view of how the samples group.
suppressMessages({library(DESeq2); library(data.table); library(yaml)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop("usage: Rscript supp_contrasts.R <config.yaml>")
PIPE <- Sys.getenv("DGE_PIPELINE_DIR")   # package root (falls back to this script's own directory)
if (!nzchar(PIPE)) { f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); PIPE <- if (length(f)) dirname(normalizePath(f[1])) else "." }
source(file.path(PIPE, "lib_dge.R"))
assert_blas_pin()
cfg <- read_config(args[1]); BASE <- Sys.getenv("BASE", cfg$base); ODIR <- file.path(BASE, cfg$output_dir)
panel <- fread(file.path(BASE, cfg$panel))
# Panel annotation keyed by feature: the gene itself when transcripts were summed per gene, the
# transcript otherwise. (The `else` must stay on the first line, because at the top level R treats
# an `if` as complete at a line break.)
ann <- if (isTRUE(cfg$panel_aggregate)) unique(panel[, .(feature=gene, gene, family)]) else
       panel[, .(feature=transcript_id, gene, family)]
if (nrow(ann) != uniqueN(ann$feature)) stop("panel annotation is not one row per feature (a gene listed with two families?)")

for (uname in names(cfg$units)) {
  f <- file.path(ODIR, paste0(uname, "_dds.rds")); if (!file.exists(f)) stop("missing fit: ", f, " (run dge_panel.R first)")
  obj <- readRDS(f); dds <- obj$dds; vsd <- obj$vst; meta <- as.data.table(obj$meta)
  cts <- contrasts_for_unit(cfg$units[[uname]])
  # One table per comparison. ct[1] is the later (or treated) stage, ct[2] the reference, so a
  # positive log2 fold change means higher expression in ct[1].
  for (ct in c(cts$transitions, cts$extra)) {
    res <- lfcShrink(dds, contrast=c("stage", ct[1], ct[2]),
                     res=results(dds, contrast=c("stage", ct[1], ct[2]), cooksCutoff=FALSE), type=cfg$shrinkage)
    dt <- as.data.table(as.data.frame(res), keep.rownames="feature")
    dt <- merge(dt, ann, by="feature", all.x=TRUE)
    # is_rtk marks the panel features (the column name is fixed whatever the panel is).
    # rtk_sig marks the panel features whose transcriptome-wide adjusted p is below alpha in THIS comparison.
    # These are the genes starred and labeled on the supplementary plots.
    dt[, is_rtk := !is.na(gene)]
    dt[, rtk_sig := is_rtk & !is.na(padj) & padj < cfg$alpha]
    dt[, unit := uname]; dt[, contrast := paste0(ct[1], "_vs_", ct[2])]
    setnames(dt, "unit", cfg$unit_col)
    setcolorder(dt, c("feature","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","gene","family","is_rtk","rtk_sig", cfg$unit_col, "contrast"))
    fwrite(dt, file.path(ODIR, sprintf("%s_dge_%s_vs_%s.tsv", uname, ct[1], ct[2])), sep="\t")
    cat(sprintf("[%s] %-4s vs %-4s: %2d panel genes significant (padj < %g) -> %s\n", uname, ct[1], ct[2],
                dt[rtk_sig == TRUE, .N], cfg$alpha, paste(dt[rtk_sig == TRUE, gene], collapse=", ")))
  }
  # PCA on the 500 most variable features of the variance-stabilized counts (the usual DESeq2
  # practice). <unit>_pca.tsv has the coordinates of each sample on PC1 and PC2 with its stage;
  # <unit>_pca_var.tsv has the percentage of variance carried by each component.
  m <- assay(vsd); rv <- apply(m, 1, var)
  sel <- order(rv, decreasing=TRUE)[seq_len(min(500, nrow(m)))]
  pc <- prcomp(t(m[sel, ])); pctv <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  stg <- meta$stage[match(colnames(m), meta$sampleid)]
  if (anyNA(stg)) stop(uname, ": PCA: sample(s) absent from the fit's metadata: ", paste(colnames(m)[is.na(stg)], collapse=","))
  fwrite(data.table(sample=colnames(m), PC1=pc$x[, 1], PC2=pc$x[, 2],
                    stage=as.character(stg)),
         file.path(ODIR, paste0(uname, "_pca.tsv")), sep="\t")
  fwrite(data.table(PC=paste0("PC", seq_along(pctv)), pct=pctv), file.path(ODIR, paste0(uname, "_pca_var.tsv")), sep="\t")
}
cat("SUPP CONTRASTS DONE\n")
