#!/bin/bash
#SBATCH --job-name=merge_expr_meta
#SBATCH --output=slurm_output/slurm-%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=kyshih@stanford.edu
#SBATCH --account=mwinslow

source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate protein
parent_dir="/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"

  #--expr ${parent_dir}/data/tcga_xena/raw/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena \
python ${parent_dir}/scripts/python/03_merge_expr_and_metadata.py \
  --expr ${parent_dir}/data/tcga_xena/intermediate/tcga_expr_truncated_barcodes.csv.gz \
  --meta ${parent_dir}/data/tcga_xena/intermediate/tcga_sample_metadata_for_merge.csv \
  --out ${parent_dir}/data/tcga_xena/processed/tcga_expr_annotated_core.csv \
  # --test_rows 100
