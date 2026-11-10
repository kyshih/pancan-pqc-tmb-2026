#!/bin/bash
#SBATCH --job-name=calculate_tmb
#SBATCH --output=slurm_output/calculate_tmb_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account=mwinslow

set -e

# --- Load environment ---
module load python  # adjust to match your SCG environment

# --- Paths ---
PARENT_DIR="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
SCRIPT_DIR="${PARENT_DIR}/scripts/python"
OUT_DIR="${PARENT_DIR}/data/tcga_xena/processed"
DATA_DIR="/labs/mwinslow/Karen/TMBandTS"

INPUT_MAF="${DATA_DIR}/Data/TCGA/Mutations/mc3.v0.2.8.PUBLIC.maf.gz"
OUTPUT_FILE="${OUT_DIR}/tcga_tmb_summary.csv"

echo "🚀 Starting TMB calculation job..."

python3 "${SCRIPT_DIR}/05_calculate_tmb.py" \
  --input "${INPUT_MAF}" \
  --output "${OUTPUT_FILE}"

date
echo "✅ Job complete: results saved to ${OUTPUT_FILE}"
