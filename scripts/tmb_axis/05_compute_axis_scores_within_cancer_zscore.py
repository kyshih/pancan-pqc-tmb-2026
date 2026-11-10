import pandas as pd
import argparse
from pathlib import Path

# ----------------------------
# Argument parser
# ----------------------------
def get_args():
    parser = argparse.ArgumentParser(description="Compute sample-level TMB, PQC, and non-PQC axis scores.")
    parser.add_argument("--expr", required=True, help="Path to regression-ready TCGA expression table")
    parser.add_argument("--tmb_weights", required=True, help="Path to GLMM beta estimates")
    parser.add_argument("--pqc_genes", required=True, help="Path to PQC gene list")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--test_rows", type=int, default=None)
    return parser.parse_args()


def main():
    args = get_args()

    expr_path = Path(args.expr)
    tmb_weights_path = Path(args.tmb_weights)
    pqc_genes_path = Path(args.pqc_genes)
    out_path = Path(args.out)

    out_path.mkdir(exist_ok=True, parents=True)

    # ------------------------
    # Load data
    # ------------------------
    expr = pd.read_csv(expr_path)
    print(f"Loaded expression: {expr.shape}")

    tmb_weights = pd.read_csv(tmb_weights_path)
    print(f"Loaded TMB weights: {tmb_weights.shape}")

    pqc_genes = pd.read_csv(pqc_genes_path, sep="\t")
    print(f"Loaded PQC genes: {pqc_genes.shape}")

    # ------------------------
    # Fix PQC gene symbols
    # ------------------------
    pqc_genes.loc[pqc_genes["Hugo"] == "TOR4A", "Hugo"] = "C9orf167"
    pqc_genes.loc[pqc_genes["Hugo"] == "FAXC", "Hugo"] = "C6orf168"

    # ------------------------
    # Define gene sets
    # ------------------------
    axis_genes = set(tmb_weights["gene"])
    pqc_gene_set = set(pqc_genes["Hugo"])
    expr_genes = set(expr.columns)

    axis_all = list(axis_genes & expr_genes)
    pqc_axis_genes = list(axis_genes & pqc_gene_set & expr_genes)
    non_pqc_genes = list((axis_genes - pqc_gene_set) & expr_genes)

    print(f"{len(axis_all)} total axis genes")
    print(f"{len(pqc_axis_genes)} PQC axis genes")
    print(f"{len(non_pqc_genes)} non-PQC axis genes")

    # ------------------------
    # Build beta table
    # ------------------------
    tmb_weights = tmb_weights.drop_duplicates(subset="gene")
    beta_df = tmb_weights.set_index("gene")

    # 1. z-score ONCE
    X = expr[axis_all]
    X = X[beta_df.index]
    
    # Keep cancer type separately
    cancer_type = expr["cancer_type"]
    
    # ------------------------
    # Within-cancer z-score
    # ------------------------
    X_within_z = X.groupby(cancer_type).transform(
        lambda x: (x - x.mean()) / x.std()
    )

    # Handle edge cases (e.g. std = 0)
    X_within_z = X_within_z.fillna(0)

    # align betas to SAME X_within_z columns
    beta_all = beta_df.loc[X_within_z.columns, "estimate"]

    # PQC axis: keep only PQC genes
    beta_pqc = beta_all.copy()
    beta_pqc.loc[~beta_pqc.index.isin(pqc_axis_genes)] = 0
    
    # non-PQC axis: keep only non-PQC genes
    beta_non_pqc = beta_all.copy()
    beta_non_pqc.loc[~beta_non_pqc.index.isin(non_pqc_genes)] = 0

    print("PQC overlap:", len(set(pqc_axis_genes) & set(X_within_z.columns)))
    print("non-PQC overlap:", len(set(non_pqc_genes) & set(X_within_z.columns)))

    # compute scores
    tmb_axis = X_within_z @ beta_all
    pqc_axis = X_within_z @ beta_pqc
    non_pqc_axis = X_within_z @ beta_non_pqc

    # ------------------------
    # Save results
    # ------------------------
    axis_df = pd.DataFrame({
        "barcode": expr["barcode"],
        "cancer_type": expr["cancer_type"],
        "log_mut": expr["log_mut"],
        "TMB_axis_score": tmb_axis,
        "pqc_axis_score": pqc_axis,
        "non_pqc_axis_score": non_pqc_axis,
        "dataset": "tcga_xena"
    })

    out_file = out_path / "axis_sample_level.parquet"
    axis_df.to_parquet(out_file, index=False)

    print(f"Saved → {out_file}")


if __name__ == "__main__":
    main()