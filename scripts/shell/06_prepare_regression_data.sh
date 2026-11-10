#!/bin/bash
#SBATCH --job-name=prepare_regression_data
#SBATCH --output=slurm_output/prepare_regression_data_%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --account=mwinslow

# --- Load environment ---
module load python

# --- Paths ---
PARENT_DIR="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
SCRIPT_DIR="${PARENT_DIR}/scripts/python"
PROCESSED_DIR="${PARENT_DIR}/data/tcga_xena/processed"
DATA_DIR="/labs/mwinslow/Karen/TMBandTS"

# INPUT_EXPR="${PARENT_DIR}/data/tcga_xena/processed/tcga_expr_filtered_zscore.csv" # z-scored expression 
INPUT_EXPR="${PARENT_DIR}/data/tcga_xena/intermediate/tcga_expr_annotated_core.csv" # raw expression
INPUT_TMB="${PROCESSED_DIR}/tcga_tmb_summary.csv"
INPUT_IMPURITY="${DATA_DIR}/Data/TCGA/MetaData/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
OUT_DIR="${PARENT_DIR}/data/tcga_xena/"

echo "--- prepares TCGA Pan-Cancer expression data for regression analyses ---"

python3 "${SCRIPT_DIR}/06_prepare_regression_data.py" \
    --input_expr "${INPUT_EXPR}" \
    --input_tmb "${INPUT_TMB}" \
    --input_purity "${INPUT_IMPURITY}" \
    --outdir "${OUT_DIR}"