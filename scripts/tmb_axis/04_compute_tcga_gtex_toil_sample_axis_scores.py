import pandas as pd
import argparse
import os


def get_args():
    parser = argparse.ArgumentParser(description="Compute TMB and PQC axis scores for normal tissues.")

    parser.add_argument("--gtex_expr", required=True, help="GTEx expression parquet")
    parser.add_argument("--tcga_expr", required=True, help="TCGA expression parquet")

    parser.add_argument("--tmb_axis", required=True, help="TMB axis gene weights (tsv)")
    parser.add_argument("--pqc_genes", required=True, help="PQC gene list (tsv)")

    parser.add_argument("--out_gtex", required=True, help="Output parquet for GTEx normal tissue axis scores")
    parser.add_argument("--out_tcga", required=True, help="Output parquet for TCGA normal tissue axis scores")

    return parser.parse_args()


def compute_axis(df, axis_genes, beta_df, pqc_gene_list):

    common_genes = list(set(axis_genes) & set(df.columns))

    beta_sub = beta_df.loc[common_genes, ["estimate"]]

    X = df[common_genes]

    # z-score
    X_z = (X - X.mean()) / X.std()

    beta = beta_sub["estimate"].reindex(X_z.columns)

    tmb_axis_score = X_z @ beta

    # PQC
    common_pqc = df.columns.intersection(pqc_gene_list)
    print(f"{len(common_pqc)} PQC genes overlap")

    X_pqc = df[common_pqc]

    X_pqc_z = (X_pqc - X_pqc.mean()) / X_pqc.std()

    pqc_axis = X_pqc_z.mean(axis=1)

    return tmb_axis_score, pqc_axis


def main():

    args = get_args()

    print("Loading expression matrices...")

    df_gtex = pd.read_parquet(args.gtex_expr)
    df_tcga = pd.read_parquet(args.tcga_expr)

    df_gtex_normal = df_gtex[df_gtex["sample_type"] == "Normal Tissue"].copy()
    df_tcga_normal = df_tcga[df_tcga["sample_type"] == "Solid Tissue Normal"].copy()

    print("Loading axis weights and gene sets...")

    tmb_axis = pd.read_csv(args.tmb_axis, sep="\t")
    pqc_genes = pd.read_csv(args.pqc_genes, sep="\t")

    pqc_gene_list = pqc_genes["gene_id"]
    axis_genes = tmb_axis["gene_id"]

    beta_df = tmb_axis.set_index("gene_id")

    print("Computing GTEx axis scores...")

    tmb_score, pqc_score = compute_axis(
        df_gtex_normal,
        axis_genes,
        beta_df,
        pqc_gene_list
    )

    df_axis_gtex = pd.DataFrame({
        "barcode": df_gtex_normal["sample_id"],
        "tissue": df_gtex_normal["tissue"],
        "TMB_axis_score": tmb_score,
        "pqc_axis_score": pqc_score
    })

    df_axis_gtex.to_parquet(args.out_gtex, index=False)

    print("Computing TCGA axis scores...")

    tmb_score, pqc_score = compute_axis(
        df_tcga_normal,
        axis_genes,
        beta_df,
        pqc_gene_list
    )

    df_axis_tcga = pd.DataFrame({
        "barcode": df_tcga_normal["sample_id"],
        "tissue": df_tcga_normal["tissue"],
        "TMB_axis_score": tmb_score,
        "pqc_axis_score": pqc_score
    })

    df_axis_tcga.to_parquet(args.out_tcga, index=False)

    print("Done.")


if __name__ == "__main__":
    main()