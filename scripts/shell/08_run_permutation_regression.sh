#!/bin/bash
#SBATCH --job-name=tcga_regression
#SBATCH --output=slurm_output/tcga_regression_perm_%j.out
#SBATCH --time=96:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --account=mwinslow

set -e

# --- Load R environment ---
source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate protein

module load R

# --- Paths ---
PARENT_DIR="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
R_SCRIPT_DIR="${PARENT_DIR}/scripts/r"
DATA_DIR="${PARENT_DIR}/data/tcga_xena/data_ready"

RESULTS_DIR="${PARENT_DIR}/results/regression/permuted"

# OUTPUT_GLOBAL="${RESULTS_DIR}/results_global.csv"
OUTPUT_GENEWISE_PERM="${RESULTS_DIR}/results_by_gene_permuted.csv"
OUTPUT_DATA_LONG_PERM="${PARENT_DIR}/data/tcga_xena/data_ready/null/tcga_expr_long_for_regression_imputed_permuted.parquet"

INPUT_DATA_LONG="${PARENT_DIR}/data/tcga_xena/data_ready/long/tcga_expr_long_for_regression_imputed.parquet"

echo "🚀 Starting TCGA permutation regression job..."
date

#   --genes TP53,BRCA1,PTEN
Rscript "${R_SCRIPT_DIR}/08_run_permutation_regression.R" \
  --input_data_long "$INPUT_DATA_LONG" \
  --output_data_long_perm "$OUTPUT_DATA_LONG_PERM" \
  --output_genewise_perm "$OUTPUT_GENEWISE_PERM" \


date
echo "✅ Job complete: results saved to ${RESULTS_DIR}"
