import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description="Compute TMB (non-PQC) and PQC axis scores for normal tissues."
    )

    parser.add_argument("--gtex_expr", required=True)
    parser.add_argument("--tcga_expr", required=True)

    parser.add_argument("--tmb_axis", required=True)
    parser.add_argument("--pqc_genes", required=True)

    parser.add_argument("--out_gtex", required=True)
    parser.add_argument("--out_tcga", required=True)

    return parser.parse_args()


def compute_axis(df, axis_genes, beta_df, pqc_gene_list):

    # --- Define gene sets ---
    axis_genes_set = set(axis_genes)
    pqc_genes_set = set(pqc_gene_list)
    df_genes = set(df.columns)

    # Intersections with data
    axis_all = list(axis_genes_set & df_genes)
    pqc_axis_genes = list(axis_genes_set & pqc_genes_set & df_genes)
    non_pqc_genes = list((axis_genes_set - pqc_genes_set) & df_genes)

    print(f"{len(axis_all)} total axis genes used")
    print(f"{len(pqc_axis_genes)} PQC axis genes")
    print(f"{len(non_pqc_genes)} non-PQC axis genes")

    # 1. z-score ONCE
    X = df[axis_all]
    X = X[beta_df.index]
    X_z = (X - X.mean()) / X.std()
    X_z = X_z.fillna(0)

    # align betas to SAME X_z columns
    beta_all = beta_df.loc[X_z.columns, "estimate"]
    
    # PQC axis: keep only PQC genes
    beta_pqc = beta_all.copy()
    beta_pqc.loc[~beta_pqc.index.isin(pqc_axis_genes)] = 0

    # non-PQC axis: keep only non-PQC genes
    beta_non_pqc = beta_all.copy()
    beta_non_pqc.loc[~beta_non_pqc.index.isin(non_pqc_genes)] = 0

    print("PQC overlap:", len(set(pqc_axis_genes) & set(X_z.columns)))
    print("non-PQC overlap:", len(set(non_pqc_genes) & set(X_z.columns)))
    
    # compute scores
    tmb_axis = X_z @ beta_all
    pqc_axis = X_z @ beta_pqc
    non_pqc_axis = X_z @ beta_non_pqc

    return tmb_axis, pqc_axis, non_pqc_axis

def main():

    args = get_args()

    print("Loading expression matrices...")
    df_gtex = pd.read_parquet(args.gtex_expr)
    df_tcga = pd.read_parquet(args.tcga_expr)

    # --- Filter normal tissues ---
    df_gtex_normal = df_gtex[df_gtex["sample_type"] == "Normal Tissue"].copy()
    df_tcga_normal = df_tcga[df_tcga["sample_type"] == "Solid Tissue Normal"].copy()

    print("Loading axis weights and gene sets...")
    tmb_axis = pd.read_csv(args.tmb_axis, sep="\t")
    pqc_genes = pd.read_csv(args.pqc_genes, sep="\t")

    tmb_axis = tmb_axis.drop_duplicates(subset="gene_id")
    axis_genes = tmb_axis["gene_id"]
    pqc_gene_list = pqc_genes["gene_id"]

    beta_df = tmb_axis.set_index("gene_id")

    # =========================
    # GTEx
    # =========================
    print("Computing GTEx axis scores...")

    tmb_score, pqc_score, non_pqc_score = compute_axis(
        df_gtex_normal,
        axis_genes,
        beta_df,
        pqc_gene_list
    )

    df_axis_gtex = pd.DataFrame({
        "barcode": df_gtex_normal["sample_id"],
        "tissue": df_gtex_normal["tissue"],
        "TMB_axis_score": tmb_score,
        "pqc_axis_score": pqc_score,
        "non_pqc_axis_score": non_pqc_score
    })

    df_axis_gtex.to_parquet(args.out_gtex, index=False)

    # =========================
    # TCGA
    # =========================
    print("Computing TCGA axis scores...")

    tmb_score, pqc_score, non_pqc_score = compute_axis(
        df_tcga_normal,
        axis_genes,
        beta_df,
        pqc_gene_list
    )

    df_axis_tcga = pd.DataFrame({
        "barcode": df_gtex_normal["sample_id"],
        "tissue": df_gtex_normal["tissue"],
        "TMB_axis_score": tmb_score,
        "pqc_axis_score": pqc_score,
        "non_pqc_axis_score": non_pqc_score
    })

    df_axis_tcga.to_parquet(args.out_tcga, index=False)

    print("Done.")


if __name__ == "__main__":
    main()