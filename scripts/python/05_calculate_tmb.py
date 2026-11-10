"""
Step 5:
Compute Tumor Mutational Burden (TMB) from TCGA MC3 MAF file.

TMB is defined here as the total number of missense, nonsense,
and silent mutations in protein-coding genes per tumor sample.

Input:
    mc3.v0.2.8.PUBLIC.maf.gz

Output:
    tcga_tmb_summary.csv
"""

import pandas as pd
import argparse
import numpy as np
import os

def main():
    # === 1. Parse arguments ===
    parser = argparse.ArgumentParser(description="Compute Tumor Mutational Burden (TMB) from MC3 MAF file.")
    parser.add_argument("--input", required=True, help="Path to MC3 MAF file (mc3.v0.2.8.PUBLIC.maf.gz)")
    parser.add_argument("--output", required=True, help="Path to output CSV file (e.g., tcga_tmb_summary.csv)")
    parser.add_argument("--test_rows", type=int, default=None, help="Optional: read only first N rows for testing.")
    args = parser.parse_args()

    # === 2. Load relevant mutation columns ===
    cols_to_use = ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"]
    print(f"🔹 Loading MAF file: {args.input}")
    mut_df = pd.read_csv(args.input, sep='\t', usecols=cols_to_use, nrows=args.test_rows)
    print(f"   Loaded {len(mut_df):,} rows")

    # === 3. Subset to coding mutation types for TMB ===
    mut_classes = ["Missense_Mutation", "Silent", "Nonsense_Mutation"]
    mut_df_filtered = mut_df[mut_df["Variant_Classification"].isin(mut_classes)].copy()
    print(f"   Retained {len(mut_df_filtered):,} coding mutations after filtering")

    # === 4. Rename columns to cleaner snake_case ===
    mut_df_filtered.rename(
        columns={
            "Hugo_Symbol": "gene",
            "Variant_Classification": "variant_classification",
            "Tumor_Sample_Barcode": "barcode",
        },
        inplace=True,
    )
    mut_df_filtered["barcode"] = mut_df_filtered["barcode"].str[:15]

    # === 5. Summarize per tumor sample ===
    tmb = (
        mut_df_filtered.groupby("barcode")
        .size()
        .reset_index(name="num_mut")
    )

    tmb['log_mut'] = np.log10(tmb['num_mut'] + 1)
    tmb['num_mut_per_mb'] = tmb['num_mut'] / 30 # assume 30 Mb coding region

    # === 6. Save results ===
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    tmb.to_csv(args.output, index=False)

    # === 7. Summary log ===
    print(f"✅ TMB table saved to: {args.output}")
    print(f"   {len(tmb):,} tumor samples processed.")
    print(f"   Mean TMB: {tmb['num_mut_per_mb'].mean():.2f} mutations/Mb")
    print(f"   Median TMB: {tmb['num_mut_per_mb'].median():.2f} mutations/Mb")
                    
if __name__ == "__main__":
    main()
    