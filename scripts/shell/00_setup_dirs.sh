#!/bin/bash
set -euo pipefail

# Base directory (adjust to your path)
BASE_DIR="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"

# Define subdirectories
DATA_DIR="${BASE_DIR}/data/tcga_xena"
RESULTS_DIR="${BASE_DIR}/results"

# Create subfolders
mkdir -p "${DATA_DIR}/raw"
mkdir -p "${DATA_DIR}/intermediate"
mkdir -p "${DATA_DIR}/processed"
mkdir -p "${DATA_DIR}/data_ready"
mkdir -p "${RESULTS_DIR}/regression"


echo "✅ Directory structure created under: ${BASE_DIR}"
