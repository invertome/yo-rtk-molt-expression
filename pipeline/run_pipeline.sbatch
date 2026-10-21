#!/bin/bash
#SBATCH --job-name=dge_pipeline
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
# Slurm wrapper around run_pipeline.sh: it only hands the config path on to that script.
# The resource lines above are a starting point (4 cores and 48 GB suffice for datasets of 50 to
# 115 libraries, and a run takes about ten minutes). The account, the partition and the log
# paths differ between clusters and go on the sbatch command line. The variables run_pipeline.sh
# reads (DGE_PIPELINE_DIR, CONDA_SH, CONDA_ENV, OPENBLAS_CORETYPE, BASE, PIPELINE_OUT) travel with
# --export. See README.md, "Running".
exec bash "${DGE_PIPELINE_DIR:?pass DGE_PIPELINE_DIR with --export}/run_pipeline.sh" "$@"
