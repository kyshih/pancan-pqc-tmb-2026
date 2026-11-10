#!/bin/bash
#SBATCH --job-name=tcga_regression
#SBATCH --output=slurm_output/tcga_regression_%j.out
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

RESULTS_DIR="${PARENT_DIR}/results/regression/observed"

INPUT_FILE="${DATA_DIR}/tcga_regression_ready_imputed.csv"
OUTPUT_GLOBAL="${RESULTS_DIR}/results_global.csv"
OUTPUT_GENEWISE="${RESULTS_DIR}/results_by_gene.csv"
OUTPUT_DATA_LONG="${PARENT_DIR}/data/tcga_xena/data_ready/long/tcga_expr_long_for_regression_imputed.parquet"

echo "🚀 Starting TCGA regression job..."
date

#   --genes TP53,BRCA1,PTEN
Rscript "${R_SCRIPT_DIR}/07_run_regressions.R" \
  --input "${INPUT_FILE}" \
  --output_data_long "${OUTPUT_DATA_LONG}" \
  --output_global "${OUTPUT_GLOBAL}" \
  --output_genewise "${OUTPUT_GENEWISE}"


date
echo "✅ Job complete: results saved to ${RESULTS_DIR}"
