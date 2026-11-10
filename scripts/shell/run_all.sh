#!/bin/bash
# run_all.sh — orchestrates TCGA regression pipeline

set -e  # stop on error
date
echo "🚀 Starting full TCGA regression pipeline..."

# Step 1. Data preparation
bash 01_download_and_load_raw.sh

# Step 2. Metadata processing
bash 02_process_sample_metadata.sh

# Step 3. Merge expression + metadata
bash 03_merge_expr_and_metadata.sh

# Step 4. Z-score normalization and filtering
bash 04_normalize_and_filter.sh

# Step 5. Imputation and regression ready data
bash 05_prepare_regression_input.sh

# Step 6. Run regressions
bash 06_run_regressions.sh

echo "✅ All steps completed!"
date



mkdir -p "${RESULTS_DIR}" slurm_output