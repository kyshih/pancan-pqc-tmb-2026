#!/usr/bin/env python3
"""
Step 6:
prepares TCGA Pan-Cancer expression data for regression analyses.
    It performs the following steps:
        1. Loads filtered and z-scored expression data.
        2. Keeps only primary tumor samples.
        3. Merges tumor mutation burden (TMB) data.
        4. Merges ABSOLUTE purity estimates.
        5. Summarizes missing purity fractions across cancer types.
        6. Creates both non-imputed and median-imputed (per cancer type) purity versions.
        7. Outputs regression-ready datasets.

Input:
    tcga_expr_filtered_zscore.csv
    TCGA_mastercalls.abs_tables_JSedit.fixed.txt
    tcga_tmb_summary

Output:
    tcga_regression_ready_no_impute.csv
    tcga_regression_ready_imputed.csv
    tcga_purity_missing_summary.csv
"""

import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser(
        description="Filter, merge, and prepare TCGA Pan-Cancer expression data for regression analysis."
    )
    parser.add_argument("--input_expr", required=True, help="Path to the filtered and z-scored TCGA Pan-Cancer expression CSV.")
    parser.add_argument("--input_tmb", required=True, help="Path to TMB CSV file.")
    parser.add_argument("--input_purity", required=True, help="Path to sample purity file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--test_rows", type=int, default=None, help="Read only first N rows (for testing).")
    args = parser.parse_args()

    # 1. Load data
    # Load the filtered and z-scored expression data
    df = pd.read_csv(args.input_expr, nrows=args.test_rows)
    print(f"Loaded matrix: {df.shape[0]} samples × {df.shape[1]} genes")
    # Load tmb per tumor sample
    tmb = pd.read_csv(args.input_tmb)
    tmb_slim = tmb[['barcode', 'log_mut']].copy()
    # Load purity
    purity = pd.read_csv(args.input_purity, sep='\t')
    purity_slim = purity[['array', 'purity']].rename(columns={'array': 'barcode'})
    
    # 2. Only keep primary tumors
    df_primary = df[df['sample_type_id']==1.0].copy()
    print(f"Sample types retained: {df_primary['sample_type'].unique()}")
    
    # 3. Remove these columns to get ready for regression
    metadata_cols = ["sample_type_id", "sample_type"]
    df_primary = df_primary.drop(columns=metadata_cols)
    print(f"Excluding metadata columns: {metadata_cols}")
 
    # 5. Merge tmb with expression data
    df_primary = df_primary.merge(tmb_slim, on='barcode', how='inner')
    print(f"After merging TMB: {df_primary.shape[0]} samples remain")
    
    # 6. merge purity with expresion data
    df_primary = df_primary.merge(purity_slim, on='barcode', how='left')
    
    # 7. impute missing purity with median
    df_primary_imputed = df_primary.copy()
    df_primary_imputed['purity'] = df_primary_imputed.groupby('cancer_type')['purity'].transform(lambda x: x.fillna(x.median()))

    # 8. summrize missing purity across cancer types
    purity_missing_summary = df_primary.groupby("cancer_type")["purity"].apply(lambda x: x.isna().mean()).reset_index(name="missing_fraction").sort_values("missing_fraction", ascending=False) # fraction missing
    # 9. save the results
    data_ready_dir = os.path.join(args.outdir, "data_ready")
    processed_dir = os.path.join(args.outdir, "processed")
    
    os.makedirs(data_ready_dir, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)

    # out_non_imputed = os.path.join(data_ready_dir, "tcga_regression_ready_no_impute.csv") # z-scored expression 
    # out_imputed = os.path.join(data_ready_dir, "tcga_regression_ready_imputed.csv") # z-scored expression 
    out_non_imputed = os.path.join(processed_dir, "tcga_expression_raw_sample_by_gene_no_impute.csv") # raw expression
    out_imputed = os.path.join(processed_dir, "tcga_expression_raw_sample_by_gene_imputed.csv") # raw expression
    out_purity_summary = os.path.join(processed_dir, "tcga_purity_missing_summary.csv")

    df_primary.to_csv(out_non_imputed, index=False)
    df_primary_imputed.to_csv(out_imputed, index=False)
    # purity_missing_summary.to_csv(out_purity_summary, index=False)

    print(f"✅ Saved non-imputed file: {out_non_imputed}")
    print(f"✅ Saved imputed file:     {out_imputed}")
    print(f"✅ Saved purity summary:   {out_purity_summary}")

    print(f"✅ Saved regression-ready data to: {data_ready_dir}")
    print(f"✅ Saved purity missing summary to: {processed_dir}")

if __name__ == "__main__":
    main()
    