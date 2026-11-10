#!/bin/bash
#SBATCH --job-name=merge_expr_meta
#SBATCH --output=slurm_output/slurm-%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=4
#SBATCH --account=mwinslow
#SBATCH --output=slurm_output/slurm-%j.out

source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate protein
PARENT_DIR="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
OUT_DIR="${PARENT_DIR}/data/tcga_xena/processed"

# --test_rows 100
python ${PARENT_DIR}/scripts/python/04_zscore_and_filter_genes.py \
--input ${PARENT_DIR}/data/tcga_xena/intermediate/tcga_expr_annotated_core.csv \
--outdir ${OUT_DIR}\
