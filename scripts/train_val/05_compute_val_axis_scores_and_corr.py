#!/usr/bin/env python3
"""
Compute TMB, PQC, and non-PQC axis scores on held-out validation samples using
train-fit betas, then calculate Pearson correlation with TMB (log_mut) per cancer
type and overall.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

SCORE_COLS = ["TMB_axis_score", "pqc_axis_score", "non_pqc_axis_score"]


def get_args():
    parser = argparse.ArgumentParser(
        description="Score validation samples and compute axis-vs-TMB correlations."
    )
    parser.add_argument("--expr", required=True, help="Regression-ready expression table.")
    parser.add_argument("--split", required=True, help="Train/val split CSV.")
    parser.add_argument("--beta", required=True, help="Train-fit beta TSV (gene, estimate).")
    parser.add_argument(
        "--pqc_genes",
        required=True,
        help="Path to PQC gene annotation TSV (Hugo column).",
    )
    parser.add_argument(
        "--output_scores",
        required=True,
        help="Output CSV for validation sample axis scores.",
    )
    parser.add_argument(
        "--output_corr_all",
        required=True,
        help="Output CSV for per-cancer correlations, all score types (fig2-style long format).",
    )
    parser.add_argument(
        "--output_corr_pancan",
        required=True,
        help="Output CSV for pan-cancer (cross-cancer) correlations, all score types.",
    )
    return parser.parse_args()


def load_pqc_gene_set(pqc_genes_path):
    pqc_genes = pd.read_csv(pqc_genes_path, sep="\t")
    pqc_genes.loc[pqc_genes["Hugo"] == "TOR4A", "Hugo"] = "C9orf167"
    pqc_genes.loc[pqc_genes["Hugo"] == "FAXC", "Hugo"] = "C6orf168"
    return set(pqc_genes["Hugo"])


def corr_with_p(df, x_col, y_col):
    if len(df) < 3:
        return pd.Series({"pearson_r": np.nan, "pval": np.nan, "n_samples": len(df)})
    r, p = pearsonr(df[x_col], df[y_col])
    return pd.Series({"pearson_r": r, "pval": p, "n_samples": len(df)})


def compute_axis_scores(expr_val, beta_df, pqc_gene_set):
    axis_genes = set(beta_df["gene"])
    expr_genes = set(expr_val.columns)

    common_genes = sorted(axis_genes & expr_genes)
    pqc_axis_genes = sorted(axis_genes & pqc_gene_set & expr_genes)
    non_pqc_genes = sorted((axis_genes - pqc_gene_set) & expr_genes)

    print(f"{len(common_genes)} total axis genes")
    print(f"{len(pqc_axis_genes)} PQC axis genes")
    print(f"{len(non_pqc_genes)} non-PQC axis genes")

    beta_all = beta_df.set_index("gene").loc[common_genes, "estimate"]
    beta_pqc = beta_all.copy()
    beta_pqc.loc[~beta_pqc.index.isin(pqc_axis_genes)] = 0
    beta_non_pqc = beta_all.copy()
    beta_non_pqc.loc[~beta_non_pqc.index.isin(non_pqc_genes)] = 0

    X_val = expr_val[common_genes]
    return {
        "TMB_axis_score": X_val @ beta_all,
        "pqc_axis_score": X_val @ beta_pqc,
        "non_pqc_axis_score": X_val @ beta_non_pqc,
    }


def _corr_by_cancer(group, score_col, **_kwargs):
    return corr_with_p(group, score_col, "log_mut")


def build_corr_table(scores_df, score_col):
    grouped = scores_df.groupby("cancer_type", group_keys=False)
    try:
        corr_by_cancer = grouped.apply(
            _corr_by_cancer,
            score_col=score_col,
            include_groups=False,
        )
    except TypeError:
        # pandas < 2.2
        corr_by_cancer = grouped.apply(
            lambda g: corr_with_p(g, score_col, "log_mut"),
        )
    corr_by_cancer = corr_by_cancer.reset_index()
    corr_by_cancer["minus_log10_pval"] = -np.log10(corr_by_cancer["pval"].clip(lower=1e-300))

    overall = corr_with_p(scores_df, score_col, "log_mut")
    overall_row = pd.DataFrame(
        [
            {
                "cancer_type": "ALL",
                "pearson_r": overall["pearson_r"],
                "pval": overall["pval"],
                "n_samples": overall["n_samples"],
                "minus_log10_pval": -np.log10(max(overall["pval"], 1e-300)),
            }
        ]
    )

    return pd.concat([corr_by_cancer, overall_row], ignore_index=True)


def build_pancan_corr_table(scores_df):
    rows = []
    for score_col in SCORE_COLS:
        overall = corr_with_p(scores_df, score_col, "log_mut")
        rows.append(
            {
                "score": score_col,
                "pearson_r": overall["pearson_r"],
                "pval": overall["pval"],
                "n_samples": overall["n_samples"],
                "minus_log10_pval": -np.log10(max(overall["pval"], 1e-300)),
            }
        )
    return pd.DataFrame(rows)


def main():
    args = get_args()

    expr = pd.read_csv(args.expr)
    split_df = pd.read_csv(args.split)
    beta_df = pd.read_csv(args.beta, sep="\t").drop_duplicates(subset="gene")
    pqc_gene_set = load_pqc_gene_set(args.pqc_genes)

    val_barcodes = split_df.loc[split_df["split"] == "val", "barcode"]
    expr_val = expr[expr["barcode"].isin(val_barcodes)].copy()

    axis_scores = compute_axis_scores(expr_val, beta_df, pqc_gene_set)

    scores_df = pd.DataFrame(
        {
            "barcode": expr_val["barcode"].values,
            "cancer_type": expr_val["cancer_type"].values,
            "log_mut": expr_val["log_mut"].values,
            "purity": expr_val["purity"].values,
            **{col: values.values for col, values in axis_scores.items()},
        }
    )

    out_scores = Path(args.output_scores)
    out_scores.parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(out_scores, index=False)
    print(f"Saved validation axis scores → {out_scores}")

    corr_all_rows = []
    for score_col in SCORE_COLS:
        corr_df = build_corr_table(scores_df, score_col)
        corr_df["score"] = score_col
        corr_all_rows.append(corr_df)

    corr_all = pd.concat(corr_all_rows, ignore_index=True)
    corr_all = corr_all[
        ["cancer_type", "score", "pearson_r", "pval", "n_samples", "minus_log10_pval"]
    ]
    corr_all = corr_all.sort_values(["score", "pearson_r"], ascending=[True, False])

    out_corr_all = Path(args.output_corr_all)
    corr_all.to_csv(out_corr_all, index=False)
    print(f"Saved per-cancer validation correlations → {out_corr_all}")

    pancan_corr = build_pancan_corr_table(scores_df)
    out_corr_pancan = Path(args.output_corr_pancan)
    pancan_corr.to_csv(out_corr_pancan, index=False)
    print(f"Saved pan-cancer validation correlations → {out_corr_pancan}")

    tmb_per_cancer = corr_all.loc[
        (corr_all["score"] == "TMB_axis_score") & (corr_all["cancer_type"] != "ALL")
    ].sort_values("pearson_r", ascending=False)
    print("\nTMB axis — top per-cancer correlations:")
    print(tmb_per_cancer.head(10).to_string(index=False))

    print("\nPan-cancer validation r by score type:")
    for _, row in pancan_corr.iterrows():
        print(f"  {row['score']}: r={row['pearson_r']:.4f}, p={row['pval']:.2e}")


if __name__ == "__main__":
    main()
