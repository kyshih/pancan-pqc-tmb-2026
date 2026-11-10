#!/usr/bin/env python3
"""
Step 3:
Merge TCGA expression matrix with sample metadata.
- Reads large Xena expression file (optionally in chunks)
- Merges metadata
- Saves merged output

Preprocessing:
- Barcodes truncated to first 15 characters to ensure patient-level matching
  between TCGA expression matrix (sample-level IDs) and metadata (patient-level IDs).
"""
import pandas as pd
import argparse
from pathlib import Path

# ----------------------------
# Argument parser
# ----------------------------
def get_args():
    parser = argparse.ArgumentParser(description="Merge TCGA expression and metadata.")
    parser.add_argument("--expr", required=True, help="Path to EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena file")
    parser.add_argument("--meta", required=True, help="Path to tcga_sample_metadata_for_merge.csv file")
    parser.add_argument("--out", required=True, help="Output merged CSV path")
    parser.add_argument("--test_rows", type=int, default=None, help="Read only first N rows (for testing)")
    return parser.parse_args()

# ----------------------------
# Main
# ----------------------------
def main():
    args = get_args()
    expr_path = Path(args.expr)
    meta_path = Path(args.meta)
    out_path = Path(args.out)

    # expr = pd.read_csv(expr_path, sep='\t', nrows=args.test_rows).copy()
    expr = pd.read_csv(expr_path, nrows=args.test_rows).copy()

    df_type = pd.read_csv(meta_path)

    # # Ensure correct column names
    # if 'gene' not in expr.columns:
    #     expr.rename(columns={expr.columns[0]: 'gene'}, inplace=True)

    # # 1 Rename first column to "gene"
    # # (since Xena files often use 'sample' or no header for gene names)
    # expr = expr.rename(columns={"sample": "gene"})
    
    # # 2 Set '' as index and transpose → samples become rows, genes become columns
    # expr_t = expr.set_index("gene").T
    # expr_t.index.name = None  # optional: remove name for cleaner look
    
    # # 3 Move sample IDs (currently the index) into a column called "Barcode"
    # expr_t = expr_t.reset_index().rename(columns={"index": "barcode"})
    
    # # 4 Ensure Barcode is the first column
    # cols = ["barcode"] + [c for c in expr_t.columns if c != "barcode"]
    # expr_t = expr_t[cols]
    
    # 5 merge with cacner type metadata
    # expr_t=expr_t.merge(df_type, on='barcode', how='left')
    expr=expr.merge(df_type, on='barcode', how='left')
    
    # 6 Quick check
    print(expr.shape)
    print(expr.head(3))
    print(expr.cancer_type.unique())
    # print(expr_t.shape)
    # print(expr_t.head(3))
    # print(expr_t.cancer_type.unique())

    # expr_t.to_csv(out_path, index=None)
    expr.to_csv(out_path, index=None)

if __name__ == "__main__":
    main()

