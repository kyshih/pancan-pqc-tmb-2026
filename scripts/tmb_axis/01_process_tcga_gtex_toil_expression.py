#!/usr/bin/env python3
"""
Prepare TCGA + GTEx expression dataset for downstream analysis.

This script processes the UCSC Xena TCGA-TARGET-GTEx dataset and produces
clean expression and phenotype tables for analysis.

Inputs
------
1. TcgaTargetGtex_rsem_gene_tpm.gz
    Gene expression matrix (TPM)
    Format:
        rows   = genes
        columns = samples

2. TcgaTargetGTEX_phenotype.txt.gz
    Sample metadata from UCSC Xena

3. tcga_cancer_type_lookup_original.csv
    Mapping between TCGA study names and cancer type abbreviations
    (e.g., BRCA, LUAD, etc.)

Outputs
-------
phenotype.parquet
    Cleaned sample metadata table

expression_all.parquet
    Expression matrix with phenotype annotations

expression_tcga.parquet
    TCGA tumor samples only

expression_gtex.parquet
    GTEx normal tissue samples only

Notes
-----
This script:
- Standardizes phenotype column names
- Transposes expression matrix to sample × gene format
- Attaches phenotype metadata
- Maps TCGA samples to cancer type abbreviations
- Separates TCGA and GTEx datasets
"""

import pandas as pd
import argparse
from pathlib import Path

# ----------------------------
# Argument parser
# ----------------------------
def get_args():
    parser = argparse.ArgumentParser(description="Merge TCGA GTEx expression and phenotype metadata.")
    parser.add_argument("--expr", required=True, help="Path to TcgaTargetGtex_rsem_gene_tpm.gz file")
    parser.add_argument("--meta", required=True, help="Path to TcgaTargetGTEX_phenotype.txt.gz file")
    parser.add_argument("--type_lookup", required=True, help="Path to tcga_cancer_type_lookup_orignal.csv")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--test_rows", type=int, default=None, help="Read only first N rows (for testing)")
    return parser.parse_args()

def main():
    args = get_args()
    expr_path = Path(args.expr)
    meta_path = Path(args.meta)
    type_lookup_path = Path(args.type_lookup)
    out_path = Path(args.out)

    out_path.mkdir(exist_ok=True, parents=True)
    
    # ------------------------------------------------
    # Load phenotype metadata
    # ------------------------------------------------
    pheno = pd.read_csv(meta_path, sep="\t",compression='gzip', encoding="latin-1")
    
    # Process phenotype data
    print(f"Loaded phenotype matrix: {pheno.shape[0]} samples × {pheno.shape[1]} features")
    # Standardize column names
    pheno = pheno.rename(columns={
    "sample": "sample_id",
    "primary disease or tissue": "disease_or_tissue",
    "_primary_site": "tissue",
    "_sample_type": "sample_type",
    "_gender": "sex",
    "_study": "dataset"})
    
    # Save cleaned phenotype table
    pheno.to_parquet(f"{out_path}/phenotype.parquet", index=False)
    
    # ------------------------------------------------
    # Load expression matrix
    # ------------------------------------------------
    print(f"Loading expression data from {args.expr} ...")
    expr = pd.read_csv(expr_path, sep="\t", index_col=0, nrows=args.test_rows).copy()
    
    print(f"Loaded matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

    # Transpose expression matrix
    # Original: gene × sample
    # Desired: sample × gene
    expr = expr.T
    expr.columns = expr.columns.str.split(".").str[0]
    expr.index.name = "sample"
    
    expr = expr.reset_index().rename(columns={"sample": "sample_id"})
    # ------------------------------------------------
    # Attach phenotype metadata
    # ------------------------------------------------
    expr = expr.merge(pheno, on="sample_id", how="left")
    
    # ------------------------------------------------
    # Split TCGA and GTEx datasets
    # ------------------------------------------------
    expr_tcga = expr.loc[expr.dataset == "TCGA"].copy()
    expr_gtex = expr.loc[expr.dataset == "GTEX"].copy() 
    
    # ------------------------------------------------
    # Load TCGA cancer type lookup table
    # ------------------------------------------------
    lookup = pd.read_csv(type_lookup_path)
    # # Normalize capitalization
    lookup["study_name"] = lookup["study_name"].str.title()
    
    # Fix naming mismatches between lookup and Xena dataset
    lookup.loc[lookup["abbreviation"]=="UCEC", "study_name"] = "Uterine Corpus Endometrioid Carcinoma"
    lookup.loc[lookup["abbreviation"]=="ACC", "study_name"] = "Adrenocortical Cancer"
    lookup.loc[lookup["abbreviation"]=="CESC", "study_name"] = "Cervical & Endocervical Cancer"
    lookup.loc[lookup["abbreviation"]=="KIRC", "study_name"] = "Kidney Clear Cell Carcinoma"
    lookup.loc[lookup["abbreviation"]=="KIRP", "study_name"] = "Kidney Papillary Cell Carcinoma"
    lookup.loc[lookup["abbreviation"]=="HNSC", "study_name"] = "Head & Neck Squamous Cell Carcinoma"
    lookup.loc[lookup["abbreviation"]=="TGCT", "study_name"] = "Testicular Germ Cell Tumor"
    lookup.loc[lookup["abbreviation"]=="PCPG", "study_name"] = "Pheochromocytoma & Paraganglioma"
    lookup.loc[lookup["abbreviation"]=="DLBC", "study_name"] = "Diffuse Large B-Cell Lymphoma"
    
    # ------------------------------------------------
    # Map TCGA samples to cancer type abbreviations
    # ------------------------------------------------
    expr_tcga = expr_tcga.merge(
    lookup,
    left_on="disease_or_tissue",
    right_on="study_name",
    how="left")
    expr_tcga = expr_tcga.drop(columns="study_name")
    
    # ------------------------------------------------
    # Save output datasets
    # ------------------------------------------------
    expr.to_parquet(f"{out_path}/expression_all.parquet", index=None, compression="zstd")
    expr_tcga.to_parquet(f"{out_path}/expression_tcga.parquet", index=None, compression="zstd")
    expr_gtex.to_parquet(f"{out_path}/expression_gtex.parquet", index=None, compression="zstd")

if __name__ == "__main__":
    main()