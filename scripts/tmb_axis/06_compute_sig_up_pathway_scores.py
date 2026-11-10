#!/usr/bin/env python3
"""
Compute beta-weighted pathway/complex scores and TMB correlations for
CORUM complexes that contain at least one significantly upregulated
TMB-axis gene (posBeta, FDR < 0.05).

Mirrors local/exploratory/29_calculate_pqc_complex_scores.ipynb and
local/figures/fig1d_pathway_corr_bs_types.ipynb.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection

METADATA_COLS = {
    "barcode",
    "cancer_type",
    "log_mut",
    "purity",
    "sample_type_id",
    "sample_type",
}


def get_args():
    parser = argparse.ArgumentParser(
        description="Score sig-up CORUM pathways and correlate with TMB."
    )
    parser.add_argument("--expr", required=True, help="Regression-ready expression CSV.")
    parser.add_argument("--tmb_weights", required=True, help="posBeta gene list with estimate column.")
    parser.add_argument("--corum_annot", required=True, help="CORUM annotation TSV (complex, Hugo).")
    parser.add_argument("--output_scores", required=True, help="Output CSV for per-sample pathway scores.")
    parser.add_argument("--output_corr", required=True, help="Output CSV for pathway-vs-TMB correlations.")
    parser.add_argument("--n_boot", type=int, default=1000, help="Bootstrap iterations for correlation.")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--test_rows", type=int, default=None)
    return parser.parse_args()


def load_corum_annotation(path):
    annot = pd.read_csv(path, sep="\t")
    gene_col = "Hugo" if "Hugo" in annot.columns else "gene"
    annot = annot.rename(columns={gene_col: "gene"})
    annot.loc[annot["gene"] == "TOR4A", "gene"] = "C9orf167"
    annot.loc[annot["gene"] == "FAXC", "gene"] = "C6orf168"
    return annot


def complexes_with_sig_up_genes(annot, sig_genes):
    sig_by_complex = (
        annot.assign(is_sig=annot["gene"].isin(sig_genes))
        .groupby("complex")["is_sig"]
        .any()
    )
    keep = sig_by_complex[sig_by_complex].index.tolist()
    print(f"{len(keep)} / {annot['complex'].nunique()} CORUM complexes have >=1 sig-up gene")
    return keep


def compute_pathway_scores(expr, beta_df, complex_to_genes):
    axis_genes = set(beta_df.index)
    expr_genes = set(expr.columns)
    axis_all = sorted(axis_genes & expr_genes)
    if not axis_all:
        raise ValueError("No axis genes overlap expression matrix.")

    X = expr[axis_all]
    X = X[beta_df.index]
    beta_all = beta_df.loc[X.columns, "estimate"]

    scores = {}
    for pathway, genes in complex_to_genes.items():
        genes_in_data = sorted(set(genes) & set(X.columns))
        if not genes_in_data:
            print(f"  {pathway}: 0 genes overlap — skipped")
            continue

        beta_pathway = beta_all.copy()
        beta_pathway.loc[~beta_pathway.index.isin(genes_in_data)] = 0
        scores[pathway] = X.to_numpy() @ beta_pathway.to_numpy()
        print(f"  {pathway}: {len(genes_in_data)} genes")

    return pd.DataFrame(scores, index=expr.index)


def bootstrap_global_corr(df_long, n_boot=1000, seed=2026):
    """Cancer-type bootstrap, matching fig1d_pathway_corr_bs_types.ipynb."""
    np.random.seed(seed)
    rows = []

    for complex_name, df_c in df_long.groupby("complex"):
        df_c = df_c.dropna(subset=["score", "log_mut"])
        if len(df_c) < 3:
            continue

        r_obs = pearsonr(df_c["score"], df_c["log_mut"])[0]
        cancer_types = df_c["cancer_type"].unique()
        n_cancer = len(cancer_types)
        groups = {ct: df_c[df_c["cancer_type"] == ct] for ct in cancer_types}

        boot_rs = []
        for _ in range(n_boot):
            sampled_types = np.random.choice(cancer_types, size=n_cancer, replace=True)
            df_boot = pd.concat([groups[ct] for ct in sampled_types])
            try:
                boot_rs.append(pearsonr(df_boot["score"], df_boot["log_mut"])[0])
            except Exception:
                continue

        if not boot_rs:
            continue

        boot_rs = np.array(boot_rs)
        n = len(boot_rs)
        p_boot = 2 * min(
            (np.sum(boot_rs <= 0) + 1) / (n + 1),
            (np.sum(boot_rs > 0) + 1) / (n + 1),
        )

        rows.append(
            {
                "complex": complex_name,
                "r": r_obs,
                "ci_low": np.percentile(boot_rs, 2.5),
                "ci_high": np.percentile(boot_rs, 97.5),
                "p_boot": p_boot,
                "n_samples": len(df_c),
                "n_sig_up_genes_in_complex": df_c["n_sig_up_genes"].iloc[0],
                "n_genes_scored": df_c["n_genes_scored"].iloc[0],
            }
        )

    corr_df = pd.DataFrame(rows)
    if not corr_df.empty:
        _, corr_df["fdr"] = fdrcorrection(corr_df["p_boot"])
        corr_df = corr_df.sort_values("r", ascending=False)
    return corr_df


def main():
    args = get_args()

    expr = pd.read_csv(args.expr, nrows=args.test_rows)
    tmb_weights = pd.read_csv(args.tmb_weights).drop_duplicates(subset="gene")
    annot = load_corum_annotation(args.corum_annot)

    sig_genes = set(tmb_weights["gene"])
    sig_complexes = complexes_with_sig_up_genes(annot, sig_genes)

    complex_to_genes = (
        annot[annot["complex"].isin(sig_complexes)]
        .groupby("complex")["gene"]
        .apply(lambda s: sorted(set(s)))
        .to_dict()
    )

    beta_df = tmb_weights.set_index("gene")
    gene_cols = [c for c in expr.columns if c not in METADATA_COLS]

    print("Computing pathway scores...")
    score_mat = compute_pathway_scores(expr[gene_cols], beta_df, complex_to_genes)

    scores_df = score_mat.copy()
    for col in ["barcode", "cancer_type", "log_mut", "purity"]:
        if col in expr.columns:
            scores_df[col] = expr[col].values

    out_scores = Path(args.output_scores)
    out_scores.parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(out_scores, index=False)
    print(f"Saved pathway scores → {out_scores}")

    meta = (
        annot[annot["complex"].isin(sig_complexes)]
        .assign(is_sig=lambda d: d["gene"].isin(sig_genes))
        .groupby("complex")
        .agg(
            n_sig_up_genes=("gene", lambda s: s[s.isin(sig_genes)].nunique()),
            n_genes_scored=("gene", "nunique"),
        )
        .reset_index()
    )

    pathway_cols = [c for c in scores_df.columns if c not in METADATA_COLS]
    id_cols = [c for c in ["barcode", "cancer_type", "log_mut", "purity"] if c in scores_df.columns]
    df_long = scores_df.melt(
        id_vars=id_cols,
        value_vars=pathway_cols,
        var_name="complex",
        value_name="score",
    ).merge(meta, on="complex", how="left")

    print(f"Computing bootstrap correlations (n_boot={args.n_boot})...")
    corr_df = bootstrap_global_corr(df_long, n_boot=args.n_boot, seed=args.seed)

    out_corr = Path(args.output_corr)
    corr_df.to_csv(out_corr, index=False)
    print(f"Saved pathway–TMB correlations → {out_corr}")


if __name__ == "__main__":
    main()
