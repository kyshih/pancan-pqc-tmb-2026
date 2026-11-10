#!/usr/bin/env python3
"""
Step 4:
Filter low-variance genes (acorss samples) and z-score normalize (within each gene) TCGA Pan-Cancer expression data.

Input:
    tcga_expr_annotated_core.csv
Output:
    tcga_expr_filtered.tsv
    tcga_expr_filtered_zscore.tsv
"""

import argparse
import pandas as pd
from scipy.stats import zscore
import numpy as np
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Filter and z-score TCGA Pan-Cancer expression data.")
    parser.add_argument("--input", required=True, help="Path to annotated expression matrix.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--var_threshold", type=float, default=0.1, help="Minimum variance threshold per gene.")
    parser.add_argument("--test_rows", type=int, default=None, help="Read only first N rows (for testing)")
    args = parser.parse_args()

    # --------------------------------------------------------------------------
    # 1. Load expression matrix
    # --------------------------------------------------------------------------
    print(f"Loading expression data from {args.input} ...")
    df = pd.read_csv(args.input, nrows=args.test_rows)
    print(f"Loaded matrix: {df.shape[0]} samples × {df.shape[1]} genes")

    # Identify metadata columns to exclude
    metadata_cols = ["cancer_type", "sample_type_id", "sample_type", "barcode"]
    print(f"Excluding metadata columns: {metadata_cols}")
    
    # Subset to expression-only columns
    expr = df.drop(columns=metadata_cols, errors="ignore")
    
    # --------------------------------------------------------------------------
    # 2. Filter genes by variance
    # --------------------------------------------------------------------------
    gene_variance = expr.var(axis=0) # variance per gene across samples
    keep = gene_variance[gene_variance > args.var_threshold].index
    expr_filtered = expr.loc[:,keep]
    print(f"Filtered out low-variance genes (threshold={args.var_threshold}); kept {len(keep)} genes")
    print(f"Filtered out {expr.shape[1] - len(keep)} low-variance genes")

    # --------------------------------------------------------------------------
    # 3. Z-score normalize each gene (row-wise)
    # --------------------------------------------------------------------------

    print("Performing z-score normalization across samples for each gene ...")
    expr_z = expr_filtered.apply(lambda x: zscore(x, ddof=1), axis=0)
    expr_z = pd.DataFrame(expr_z, index=expr_filtered.index, columns=expr_filtered.columns)
    expr_z = expr_z.replace([np.inf, -np.inf], np.nan).fillna(0)
    
    # Reattach metadata (if any)
    df_filtered = pd.concat([df[metadata_cols], expr_filtered], axis=1)
    df_z = pd.concat([df[metadata_cols], expr_z], axis=1)

    # --------------------------------------------------------------------------
    # 4. Save outputs
    # --------------------------------------------------------------------------
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    filtered_path = outdir / "tcga_expr_filtered.csv"
    zscore_path = outdir / "tcga_expr_filtered_zscore.csv"

    df_filtered.to_csv(filtered_path, index=None)
    df_z.to_csv(zscore_path, index=None)
    print(f"✅ Saved filtered data to: {filtered_path}")
    print(f"✅ Saved z-scored data to: {zscore_path}")

if __name__ == "__main__":
    main()