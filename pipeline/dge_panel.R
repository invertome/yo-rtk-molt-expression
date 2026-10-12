#!/usr/bin/env Rscript
# Main analysis script: fit the DESeq2 models and write the panel tables.
#
# Usage: Rscript dge_panel.R <config.yaml>
#
# What happens, in order, for every unit named in the configuration:
#   1. Read the sample sheet and the panel, and check both against the configuration (a wrong sheet
#      stops the run before any fit).
#   2. Import the Salmon quantifications with tximport, keeping the transcript-level estimates or
#      summing transcripts per gene when the configuration says so.
#   3. Fit a transcriptome-wide DESeq2 model with one mean per stage.
#   4. Ask, gene by gene, whether expression changes across the stages with a likelihood-ratio test
#      (full model: one mean per stage; reduced model: one shared mean, with any level listed in
#      outside_stage_test keeping its own mean). Adjust the p-values across the panel only (the
#      targeted false discovery rate) and call a panel gene regulated when that rate is below alpha.
#   5. For each regulated gene, compare every pair of stages with a Wald test, adjust within the gene
#      (Holm by default), and summarize the result as compact letters.
#   6. Write the per-stage mean and SEM of TPM, the per-sample TPM, the transcriptome-wide test table and the
#      fitted model object.
# Output filenames use the unit name and the panel_tag of the configuration, for example
# T1_panel_letters.tsv when the unit is T1 and the tag is panel.
suppressMessages({library(tximport); library(DESeq2); library(data.table); library(multcompView); library(yaml)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) stop("usage: Rscript dge_panel.R <config.yaml>")
PIPE <- Sys.getenv("DGE_PIPELINE_DIR")   # directory holding lib_dge.R; defaults to this script's own directory
if (!nzchar(PIPE)) { f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); PIPE <- if (length(f)) dirname(normalizePath(f[1])) else "." }
source(file.path(PIPE, "lib_dge.R"))
assert_blas_pin()                          # refuse to run with an unpinned, unreproducible BLAS setup
cfg  <- read_config(args[1])
BASE <- Sys.getenv("BASE", cfg$base)      # every relative path in the configuration resolves under BASE
QDIR <- file.path(BASE, cfg$quant_dir); ODIR <- file.path(BASE, cfg$output_dir); TAG <- cfg$panel_tag
dir.create(ODIR, showWarnings=FALSE, recursive=TRUE)
cat(sprintf("DESeq2 %s | tximport %s | multcompView %s | config %s | output %s\n", packageVersion("DESeq2"),
            packageVersion("tximport"), packageVersion("multcompView"), args[1], ODIR))

## ---- inputs --------------------------------------------------------------------------------
# The sample sheet has one row per library. Column names are lower-cased so that the configuration
# can name them without worrying about capitalization.
ss <- fread(file.path(BASE, cfg$sample_sheet)); setnames(ss, tolower(names(ss)))
for (col in c(cfg$sample_id_col, cfg$quant_subdir_col, cfg$level_col, cfg$unit_col))
  if (!col %in% names(ss)) stop("sample sheet lacks column ", col, "; has: ", paste(names(ss), collapse=","))
# The panel lists the genes of interest: one row per transcript, with the gene it belongs to and the
# gene family. The same transcript cannot serve two genes, and a gene cannot sit in two families.
panel <- fread(file.path(BASE, cfg$panel))
for (col in c("gene","transcript_id","family")) if (!col %in% names(panel)) stop("panel lacks column ", col)
dup_tx <- unique(panel$transcript_id[duplicated(panel$transcript_id)])
if (length(dup_tx)) stop("panel: duplicated transcript_id: ", paste(sapply(dup_tx, function(t)
  sprintf("%s (genes %s)", t, paste(unique(panel$gene[panel$transcript_id == t]), collapse=","))), collapse="; "))
genes <- unique(panel$gene)
multi_fam <- panel[, uniqueN(family), by=gene][V1 > 1, gene]
if (length(multi_fam)) stop("panel: gene(s) listed with more than one family: ", paste(multi_fam, collapse=","))
# Two ways of matching panel genes to quantified features. With panel_aggregate true, the counts of
# all transcripts of a gene are summed and the feature is the gene (used when the assembly holds
# several transcripts per gene); the transcripts that were summed are written to <name>_<tag>_panel.tsv.
# With panel_aggregate false, each panel gene is one transcript and the feature is that transcript.
if (isTRUE(cfg$panel_aggregate)) {
  first_q <- file.path(QDIR, ss[[cfg$quant_subdir_col]][1], "quant.sf")
  if (!file.exists(first_q)) stop("missing quant.sf: ", first_q)
  all_tx <- fread(first_q, select="Name")$Name
  absent <- setdiff(panel$transcript_id, all_tx)
  if (length(absent)) stop("panel transcripts absent from quant.sf: ", paste(absent, collapse=","))
  t2g <- data.table(tx=all_tx, gene=all_tx); t2g[panel, gene := i.gene, on=c(tx="transcript_id")]
  id_of <- panel[, .(id=paste(sort(unique(transcript_id)), collapse=";"), n_reps=uniqueN(transcript_id)), by=gene]
  feat_of <- setNames(genes, genes)                                  # feature name in the count matrix = gene
  fwrite(id_of[, .(gene_id=gene, rep_ids=id, n_reps)], file.path(ODIR, paste0(cfg$name, "_", TAG, "_panel.tsv")), sep="\t")
} else {
  if (anyDuplicated(panel$gene)) stop("panel_aggregate is false but a gene has several transcript_id rows: ",
                                      paste(unique(panel$gene[duplicated(panel$gene)]), collapse=","))
  id_of <- panel[, .(gene, id=transcript_id)]
  feat_of <- setNames(panel$transcript_id, panel$gene)
}
cat(sprintf("panel: %d genes, %d transcripts, aggregate=%s\n", length(genes), nrow(panel), cfg$panel_aggregate))

## ---- sample sheet per unit (every guard runs for every unit before the first fit) --------------
# Select the rows of the sample sheet that belong to one unit and check them: every stage value must
# be either one of the configured levels or one of the deliberately excluded values; every level must
# have samples; the counts per level must match expected_n when it is given; sample ids and
# quantification paths must be unique; and every quant.sf must exist. The stage column becomes a
# factor whose level order is the order given in the configuration (the order of the figures).
unit_meta <- function(uname) {
  unit <- cfg$units[[uname]]; lv <- unit$levels
  meta <- ss[get(cfg$unit_col) == uname]
  if (!nrow(meta)) stop("unit ", uname, ": no samples with ", cfg$unit_col, " == ", uname)
  odd <- table(meta[[cfg$level_col]][!(meta[[cfg$level_col]] %in% c(lv, unit$exclude_levels))], useNA="ifany")
  if (length(odd)) stop("unit ", uname, ": ", cfg$level_col, " value(s) neither in levels nor in exclude_levels: ",
                        paste(sprintf("%s (n = %d)", names(odd), as.integer(odd)), collapse=", "))
  dropped <- meta[!(get(cfg$level_col) %in% lv)]
  cat(sprintf("[%s] %d samples; dropping %d whose %s is not in levels (%s)\n", uname, nrow(meta), nrow(dropped),
              cfg$level_col, if (nrow(dropped)) paste(unique(dropped[[cfg$level_col]]), collapse=",") else "none"))
  meta <- meta[get(cfg$level_col) %in% lv]
  if (!nrow(meta)) stop("unit ", uname, ": no samples remain after keeping levels ", paste(lv, collapse=","))
  meta[, stage := factor(get(cfg$level_col), levels=lv)]
  got <- as.integer(table(meta$stage))                                 # in level order, zeros included
  if (any(got == 0L)) stop("unit ", uname, ": level(s) with zero samples: ", paste(lv[got == 0L], collapse=","))
  if (!is.null(unit$expected_n) && !identical(got, unname(unit$expected_n)))
    stop("unit ", uname, ": per-level sample counts differ from expected_n: ",
         paste(sprintf("%s got %d expected %d", lv, got, unit$expected_n), collapse=", "))
  meta[, sampleid := as.character(get(cfg$sample_id_col))]
  meta[, qpath := file.path(QDIR, get(cfg$quant_subdir_col), "quant.sf")]
  dup <- unique(meta$sampleid[duplicated(meta$sampleid)])
  if (length(dup)) stop("unit ", uname, ": duplicated ", cfg$sample_id_col, ": ", paste(dup, collapse=","))
  dup <- unique(meta$qpath[duplicated(meta$qpath)])
  if (length(dup)) stop("unit ", uname, ": duplicated quant path: ", paste(dup, collapse=","))
  miss <- meta[!file.exists(qpath)]
  if (nrow(miss)) stop(uname, ": missing quant.sf for ", paste(miss$sampleid, collapse=","))
  meta
}

## ---- one unit --------------------------------------------------------------------------------
persample_all <- list()                    # per-sample TPM of every unit, written as one table at the end
run_unit <- function(uname, meta) {
  unit <- cfg$units[[uname]]; lv <- unit$levels; outside <- unit$outside_stage_test
  files <- setNames(meta$qpath, meta$sampleid)
  cd <- as.data.frame(meta); rownames(cd) <- meta$sampleid
  txi <- if (isTRUE(cfg$panel_aggregate))
           tximport(files, type="salmon", tx2gene=as.data.frame(t2g), txOut=FALSE, dropInfReps=TRUE)
         else tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE)
  feats <- unname(feat_of)
  absent <- feats[!feats %in% rownames(txi$abundance)]
  if (length(absent)) stop(uname, ": panel features absent from the quantification: ", paste(absent, collapse=","))

  # Transcriptome-wide fit. Features are kept when at least as many samples as the smallest stage group
  # have 10 or more counts; the panel features are always kept so that every panel gene is reported
  # even when it is barely expressed. DESeq() estimates size factors, dispersions and the per-stage
  # means for the full model.
  dds <- DESeqDataSetFromTximport(txi, colData=cd, design=~stage)
  minN <- min(table(meta$stage))
  keep <- (rowSums(counts(dds) >= 10) >= minN) | (rownames(dds) %in% feats)   # force-retain the panel
  dds  <- dds[keep, ]; dds <- DESeq(dds)
  # Likelihood-ratio test for a stage effect (see reduced_matrix in lib_dge.R for the two models).
  # Cook's-distance outlier exclusion is switched off on purpose: with few samples per stage a single
  # high-count animal would otherwise remove a gene from the test without any warning in the tables.
  mm <- reduced_matrix(as.data.frame(colData(dds)), "stage", outside)
  if (mm$df != length(lv) - length(outside) - 1L)
    stop("unit ", uname, ": LRT df is ", mm$df, ", expected ", length(lv) - length(outside) - 1L, " (", length(lv),
         " levels, ", length(outside), " outside the stage test)")
  ddsL <- if (length(outside)) DESeq(dds, test="LRT", full=mm$full, reduced=mm$reduced, quiet=TRUE)
          else DESeq(dds, test="LRT", reduced=~1, quiet=TRUE)
  res_lrt <- results(ddsL, cooksCutoff=FALSE)                       # Cook's-distance exclusion deliberately off
  cat(sprintf("[%s] LRT: %d levels tested%s (df = %d); %d samples, %d features after filtering\n", uname,
              length(lv) - length(outside), if (length(outside)) paste0(", ", paste(outside, collapse=","), " kept with a separate mean") else "",
              mm$df, ncol(dds), nrow(dds)))

  # Panel summary of the test: raw p, targeted FDR (across the panel), transcriptome-wide FDR, regulated call.
  st <- targeted_stats(res_lrt, feats, cfg$alpha)
  if (any(is.na(st$lrt_p))) cat(sprintf("[%s] NOTE: LRT p is NA (zero counts) for %s; targeted BH across %d testable genes\n",
                                        uname, paste(names(feat_of)[is.na(st$lrt_p)], collapse=","), sum(!is.na(st$lrt_p))))
  # Raw p-values of every pairwise Wald comparison between two stages, for every panel feature.
  # They are adjusted later, within each regulated gene, when the letters are built.
  pairs <- combn(lv, 2, simplify=FALSE)
  pmat <- matrix(NA_real_, nrow=length(feats), ncol=length(pairs),
                 dimnames=list(feats, sapply(pairs, function(p) paste(p, collapse="-"))))
  for (j in seq_along(pairs)) {
    pr <- pairs[[j]]; r <- results(dds, contrast=c("stage", pr[2], pr[1]), cooksCutoff=FALSE)
    inr <- feats[feats %in% rownames(r)]; pmat[inr, j] <- r[inr, "pvalue"]
  }

  # Per gene: mean and SEM of TPM at each stage, TPM of every sample, and the letters. A gene that
  # is not regulated gets the same letter at every stage (no pairwise comparison is claimed for it).
  tpm <- txi$abundance; stats <- list(); letters <- list(); persample <- list()
  for (i in seq_along(genes)) {
    gene <- genes[i]; feat <- feat_of[[gene]]; id <- id_of[gene == genes[i], id]
    for (s in lv) {
      cols <- meta$sampleid[meta$stage == s]
      vals <- as.numeric(tpm[feat, cols])
      stats[[length(stats)+1]] <- data.table(gene=gene, id=id, unit=uname, stage=s, n=length(cols),
                                             mean_tpm=mean(vals), sem_tpm=sd(vals)/sqrt(length(vals)))
      persample[[length(persample)+1]] <- data.table(gene=gene, id=id, unit=uname, sampleid=cols, stage=s, tpm=vals)
    }
    reg <- st[id == feat, molt_regulated]
    let <- if (isTRUE(reg)) letters_for(pmat[feat, ], lv, cfg$pairwise_adjust, cfg$alpha, paste(uname, gene))
           else setNames(rep("a", length(lv)), lv)
    letters[[gene]] <- data.table(gene=gene, unit=uname, stage=lv, letter=unname(let),
                                  lrt_p=st[id == feat, lrt_p], q_targeted=st[id == feat, q_targeted],
                                  q_transcriptomewide=st[id == feat, q_transcriptomewide], molt_regulated=isTRUE(reg))
  }
  # Write the tables. The "unit" column is renamed to the unit column of the sample sheet (for
  # example "experiment" or "tissue") so that the outputs use the project's own vocabulary.
  stt <- rbindlist(stats); le <- rbindlist(letters); ps <- rbindlist(persample)
  setnames(stt, "unit", cfg$unit_col); setnames(le, "unit", cfg$unit_col); setnames(ps, "unit", cfg$unit_col)
  persample_all[[uname]] <<- ps
  fwrite(stt, file.path(ODIR, paste0(uname, "_", TAG, "_tpm_by_stage.tsv")), sep="\t")
  fwrite(le,  file.path(ODIR, paste0(uname, "_", TAG, "_letters.tsv")), sep="\t")
  fwrite(as.data.table(res_lrt, keep.rownames="feature"), file.path(ODIR, paste0(uname, "_deseq2_LRT_all.tsv")), sep="\t")
  # The fitted model, a variance-stabilized version of the counts (used for the PCA), and the
  # sample table are saved for supp_contrasts.R, which computes the pairwise tables from the same fit.
  saveRDS(list(dds=dds, vst=vst(dds, blind=TRUE), meta=meta), file.path(ODIR, paste0(uname, "_dds.rds")))
  regs <- unique(le[molt_regulated == TRUE, gene])
  cat(sprintf("[%s] molt-regulated (targeted FDR < %g): %d/%d -> %s\n", uname, cfg$alpha, length(regs), length(genes),
              paste(regs, collapse=", ")))
  invisible(NULL)
}
metas <- lapply(setNames(nm=names(cfg$units)), unit_meta)             # every sample-sheet guard before the first fit
for (u in names(cfg$units)) run_unit(u, metas[[u]])
fwrite(rbindlist(persample_all), file.path(ODIR, paste0(cfg$name, "_", TAG, "_persample_tpm.tsv")), sep="\t")
cat("DGE PANEL DONE\n")
