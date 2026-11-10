#!/usr/bin/env python3
"""
Sign test for per-cancer correlation direction: full data vs validation.

Tests whether Pearson r signs agree across cancer types for each score type.
"""

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest

SCORE_LABELS = {
    "TMB_axis_score": "TMB axis",
    "pqc_axis_score": "PQC axis",
    "non_pqc_axis_score": "Non-PQC axis",
}


def get_args():
    parser = argparse.ArgumentParser(
        description="Sign test for per-cancer r sign concordance (full data vs validation)."
    )
    parser.add_argument(
        "--val_per_cancer",
        required=True,
        help="Validation per-cancer correlation CSV (corr_val_axis_scores_vs_tmb.csv).",
    )
    parser.add_argument(
        "--full_per_cancer",
        required=True,
        help="Full-data per-cancer correlation CSV (e.g. score_corr_per_type_summary.csv).",
    )
    parser.add_argument(
        "--output_merged",
        required=True,
        help="Output merged per-cancer comparison CSV.",
    )
    parser.add_argument(
        "--output_sign_test",
        required=True,
        help="Output sign-test summary CSV.",
    )
    return parser.parse_args()


def load_validation_per_cancer(path):
    df = pd.read_csv(path)
    df = df.loc[df["cancer_type"] != "ALL", ["cancer_type", "score", "pearson_r", "n_samples"]]
    return df.rename(columns={"pearson_r": "r_validation", "n_samples": "n_validation"})


def load_full_per_cancer(path):
    df = pd.read_csv(path)
    if "r" in df.columns:
        r_col = "r"
    elif "pearson_r" in df.columns:
        r_col = "pearson_r"
    else:
        raise ValueError("Full-data table must contain 'r' or 'pearson_r'.")

    n_col = "n" if "n" in df.columns else "n_samples"
    df = df[["cancer_type", "score", r_col, n_col]].copy()
    return df.rename(columns={r_col: "r_full_data", n_col: "n_full_data"})


def run_sign_test(merged):
    rows = []
    for score, sub in merged.groupby("score"):
        sub = sub.copy()
        sub["sign_validation"] = sub["r_validation"] > 0
        sub["sign_full_data"] = sub["r_full_data"] > 0
        sub["same_sign"] = sub["sign_validation"] == sub["sign_full_data"]
        sub["both_positive"] = sub["sign_validation"] & sub["sign_full_data"]

        n = len(sub)
        n_same = int(sub["same_sign"].sum())
        n_both_pos = int(sub["both_positive"].sum())
        n_full_pos = int(sub["sign_full_data"].sum())
        n_val_pos = int(sub["sign_validation"].sum())

        full_pos = sub.loc[sub["sign_full_data"]]
        n_full_pos_val_pos = int(full_pos["sign_validation"].sum()) if len(full_pos) else 0

        p_concordance = binomtest(n_same, n, 0.5, alternative="greater").pvalue
        p_both_positive = binomtest(n_both_pos, n, 0.5, alternative="greater").pvalue

        rows.append(
            {
                "score": score,
                "score_label": SCORE_LABELS.get(score, score),
                "n_cancers": n,
                "n_same_sign": n_same,
                "fraction_same_sign": n_same / n,
                "p_sign_concordance": p_concordance,
                "n_both_positive": n_both_pos,
                "fraction_both_positive": n_both_pos / n,
                "p_both_positive_gt_half": p_both_positive,
                "n_full_data_positive": n_full_pos,
                "n_validation_positive": n_val_pos,
                "n_full_pos_and_val_pos": n_full_pos_val_pos,
                "n_full_data_positive": n_full_pos,
                "fraction_val_pos_given_full_pos": (
                    n_full_pos_val_pos / n_full_pos if n_full_pos else float("nan")
                ),
            }
        )
    return pd.DataFrame(rows)


def main():
    args = get_args()

    val_df = load_validation_per_cancer(args.val_per_cancer)
    full_df = load_full_per_cancer(args.full_per_cancer)

    merged = val_df.merge(full_df, on=["cancer_type", "score"], how="inner")
    merged["same_sign"] = (merged["r_validation"] > 0) == (merged["r_full_data"] > 0)
    merged["both_positive"] = (merged["r_validation"] > 0) & (merged["r_full_data"] > 0)
    merged["score_label"] = merged["score"].map(SCORE_LABELS)
    merged = merged.sort_values(["score", "cancer_type"])

    sign_df = run_sign_test(merged)

    out_merged = Path(args.output_merged)
    out_sign = Path(args.output_sign_test)
    out_merged.parent.mkdir(parents=True, exist_ok=True)
    out_sign.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(out_merged, index=False)
    sign_df.to_csv(out_sign, index=False)

    print(f"Saved merged per-cancer table → {out_merged}")
    print(f"Saved sign-test summary → {out_sign}")
    print()
    print(sign_df[
        [
            "score_label",
            "n_same_sign",
            "n_cancers",
            "p_sign_concordance",
            "n_both_positive",
            "n_full_pos_and_val_pos",
            "n_full_data_positive",
        ]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
