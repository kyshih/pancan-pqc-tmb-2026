#!/usr/bin/env python3
"""
Compute LOCO TMB axis scores on held-out cancer samples and per-cancer
correlation with log_mut (TMB).

Mirrors local/exploratory/21_loco_analysis.ipynb.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

METADATA_COLS = {"barcode", "cancer_type", "log_mut", "purity"}


def get_args():
    parser = argparse.ArgumentParser(
        description="Compute LOCO axis sample scores and per-cancer TMB correlations."
    )
    parser.add_argument("--loco_combined", required=True)
    parser.add_argument("--expr", required=True, help="Regression-ready expression CSV.")
    parser.add_argument(
        "--input_beta",
        default=None,
        help="Optional posBeta gene list CSV; only these genes are used for scoring.",
    )
    parser.add_argument("--output_scores", required=True)
    parser.add_argument("--output_corr", required=True)
    return parser.parse_args()


def corr_with_p(df, x_col, y_col):
    if len(df) < 3:
        return pd.Series({"pearson_r_tmb_vs_tmb_axis": np.nan, "pval_tmb_vs_tmb_axis": np.nan})
    r, p = pearsonr(df[x_col], df[y_col])
    return pd.Series({"pearson_r_tmb_vs_tmb_axis": r, "pval_tmb_vs_tmb_axis": p})


def main():
    args = get_args()

    beta_df = pd.read_csv(args.loco_combined)
    expr = pd.read_csv(args.expr)

    if args.input_beta:
        allowed_genes = set(pd.read_csv(args.input_beta)["gene"].unique())
        beta_df = beta_df[beta_df["gene"].isin(allowed_genes)].copy()
        print(f"Filtered to {len(allowed_genes)} genes from posBeta list")

    expr_gene_cols = [c for c in expr.columns if c not in METADATA_COLS]
    common_genes = set(expr_gene_cols) & set(beta_df["gene"])
    print(f"{len(common_genes)} genes overlap for TMB axis")

    score_frames = []
    for cancer in beta_df["held_out_cancer"].unique():
        expr_c = expr.loc[expr["cancer_type"] == cancer].copy()
        if expr_c.empty:
            print(f"Skipping {cancer} — no expression samples")
            continue

        beta_c = (
            beta_df.loc[beta_df["held_out_cancer"] == cancer, ["gene", "estimate"]]
            .drop_duplicates(subset="gene")
            .set_index("gene")
        )

        genes = beta_c.index.intersection(expr_gene_cols)
        X_c = expr_c[list(genes)]
        beta_vec = beta_c.loc[genes, "estimate"]
        tmb_axis_c = X_c.to_numpy() @ beta_vec.to_numpy()

        score_frames.append(
            pd.DataFrame(
                {
                    "barcode": expr_c["barcode"].values,
                    "cancer_type": cancer,
                    "log_mut": expr_c["log_mut"].values,
                    "purity": expr_c["purity"].values,
                    "TMB_axis_score": tmb_axis_c,
                }
            )
        )

    scores_df = pd.concat(score_frames, ignore_index=True)

    out_scores = Path(args.output_scores)
    out_scores.parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(out_scores, index=False)
    print(f"Saved LOCO axis scores → {out_scores} ({len(scores_df)} samples)")

    corr_rows = []
    for cancer, group in scores_df.groupby("cancer_type"):
        c = corr_with_p(group, "log_mut", "TMB_axis_score")
        pval = c["pval_tmb_vs_tmb_axis"]
        corr_rows.append(
            {
                "cancer_type": cancer,
                "pearson_r_tmb_vs_tmb_axis": c["pearson_r_tmb_vs_tmb_axis"],
                "pval_tmb_vs_tmb_axis": pval,
                "minus_log10_pval": -np.log10(max(pval, 1e-300)) if pd.notna(pval) else np.nan,
            }
        )

    corr_df = pd.DataFrame(corr_rows).sort_values(
        "pearson_r_tmb_vs_tmb_axis", ascending=False
    )

    out_corr = Path(args.output_corr)
    corr_df.to_csv(out_corr, index=False)
    print(f"Saved per-cancer correlations → {out_corr}")


if __name__ == "__main__":
    main()
