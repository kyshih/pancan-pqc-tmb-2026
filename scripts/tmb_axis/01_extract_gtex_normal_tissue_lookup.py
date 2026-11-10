#!/usr/bin/env python3
"""
Extract GTEx normal-tissue labels from UCSC Xena phenotype metadata.

Outputs a deduplicated lookup table of GTEx tissue labels for samples with
sample_type == "Normal Tissue".

Inputs
------
One of:
- --meta  TcgaTargetGTEX_phenotype.txt.gz (raw Xena phenotype file)
- --pheno phenotype.parquet (output from 01_process_tcga_gtex_toil_expression.py)

Outputs
-------
gtex_normal_tissue_lookup.csv
    Columns: tissue, disease_or_tissue, n_samples
"""

import argparse
from pathlib import Path

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(
        description="Extract GTEx normal tissue labels (tissue, disease_or_tissue)."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--meta",
        help="Path to TcgaTargetGTEX_phenotype.txt.gz",
    )
    input_group.add_argument(
        "--pheno",
        help="Path to processed phenotype.parquet",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output directory or file path (.csv or .parquet)",
    )
    return parser.parse_args()


def load_phenotype(meta_path=None, pheno_path=None):
    if pheno_path is not None:
        pheno = pd.read_parquet(pheno_path)
    else:
        pheno = pd.read_csv(
            meta_path,
            sep="\t",
            compression="gzip",
            encoding="latin-1",
        )
        pheno = pheno.rename(
            columns={
                "sample": "sample_id",
                "primary disease or tissue": "disease_or_tissue",
                "_primary_site": "tissue",
                "_sample_type": "sample_type",
                "_gender": "sex",
                "_study": "dataset",
            }
        )
    return pheno


def build_gtex_normal_tissue_lookup(pheno):
    gtex_normal = pheno[
        (pheno["dataset"] == "GTEX") & (pheno["sample_type"] == "Normal Tissue")
    ].copy()

    lookup = (
        gtex_normal.groupby(["tissue", "disease_or_tissue"], dropna=False)
        .size()
        .reset_index(name="n_samples")
        .sort_values(["tissue", "disease_or_tissue"])
        .reset_index(drop=True)
    )
    return lookup


def resolve_output_path(out_arg):
    out_path = Path(out_arg)
    if out_path.suffix.lower() in {".csv", ".parquet"}:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        return out_path

    out_path.mkdir(parents=True, exist_ok=True)
    return out_path / "gtex_normal_tissue_lookup.csv"


def main():
    args = get_args()
    pheno = load_phenotype(meta_path=args.meta, pheno_path=args.pheno)
    lookup = build_gtex_normal_tissue_lookup(pheno)

    print(
        f"Found {len(lookup)} GTEx normal-tissue label pairs "
        f"from {lookup['n_samples'].sum()} samples"
    )

    out_path = resolve_output_path(args.out)
    if out_path.suffix.lower() == ".parquet":
        lookup.to_parquet(out_path, index=False)
    else:
        lookup.to_csv(out_path, index=False)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
