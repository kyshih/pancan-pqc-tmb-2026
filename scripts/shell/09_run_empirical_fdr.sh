#!/bin/bash
#SBATCH --job-name=calculate_empirical_fdr
#SBATCH --output=slurm_output/calculate_empirical_fdr_%j.out
#SBATCH --time=4:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=2
#SBATCH --account=mwinslow

set -e

# --- Load R environment ---
source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate protein

module load R

# --- Paths ---
PARENT_DIR="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
R_SCRIPT_DIR="${PARENT_DIR}/scripts/r"
RESULTS_DIR="${PARENT_DIR}/results/regression"
FDR_DIR="${RESULTS_DIR}/fdr"
PLOT_DIR="${FDR_DIR}/plots"
mkdir -p "$FDR_DIR" "$PLOT_DIR"

# --- Input files ---
GENEWISE_OBSERVED="${RESULTS_DIR}/observed/results_by_gene.csv"
GENEWISE_PERMUTED="${RESULTS_DIR}/permuted/results_by_gene_permuted.csv"

# --- Output files ---
GENEWISE_FDR="${FDR_DIR}/results_by_gene_with_empirical_fdr.csv"

echo "🚀 Starting empirical FDR calculation..."
date
echo "📥 Observed:  $GENEWISE_OBSERVED"
echo "📥 Permuted:  $GENEWISE_PERMUTED"
echo "📤 Output:    $GENEWISE_FDR"
echo "📊 Plots →    $PLOT_DIR"
echo "---------------------------------------------"

# --- Run FDR calculation ---
Rscript "${R_SCRIPT_DIR}/09_calculate_empirical_fdr.R" \
  --observed "$GENEWISE_OBSERVED" \
  --permuted "$GENEWISE_PERMUTED" \
  --output "$GENEWISE_FDR"

echo "✅ Empirical FDR completed successfully!"
date