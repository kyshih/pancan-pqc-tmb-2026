#!/usr/bin/env python3
"""
Build a pan-cancer Pearson r comparison table: full-data vs train/validation.

Full-data correlations are computed from precomputed sample scores
(axis_sample_level.parquet), matching fig2c.
"""

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr

SCORE_COLS = ["TMB_axis_score", "pqc_axis_score", "non_pqc_axis_score"]
SCORE_LABELS = {
    "TMB_axis_score": "TMB axis",
    "pqc_axis_score": "PQC axis",
    "non_pqc_axis_score": "Non-PQC axis",
}


def get_args():
    parser = argparse.ArgumentParser(
        description="Compare pan-cancer axis-vs-TMB r between full data and validation."
    )
    parser.add_argument(
        "--val_pancan",
        required=True,
        help="Validation pan-cancer correlation CSV (corr_val_pancan_axis_scores_vs_tmb.csv).",
    )
    parser.add_argument(
        "--full_scores",
        required=True,
        help="Precomputed full-data sample scores (axis_sample_level.parquet or CSV).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output comparison table CSV.",
    )
    return parser.parse_args()


def pancan_corr_from_scores(scores_df):
    rows = []
    for score_col in SCORE_COLS:
        r, p = pearsonr(scores_df["log_mut"], scores_df[score_col])
        rows.append(
            {
                "score": score_col,
                "pearson_r": r,
                "pval": p,
                "n_samples": len(scores_df),
            }
        )
    return pd.DataFrame(rows)


def load_full_scores(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Full-data scores file not found: {path}")

    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)

    missing = [col for col in ["log_mut"] + SCORE_COLS if col not in df.columns]
    if missing:
        raise ValueError(f"Full-data scores file missing columns: {missing}")
    return df[["log_mut"] + SCORE_COLS].copy()


def main():
    args = get_args()

    val_df = pd.read_csv(args.val_pancan)
    val_df = val_df.rename(
        columns={
            "pearson_r": "validation_r",
            "pval": "validation_pval",
            "n_samples": "validation_n",
        }
    )

    print(f"Using precomputed full-data scores → {args.full_scores}")
    full_scores_df = load_full_scores(args.full_scores)

    full_df = pancan_corr_from_scores(full_scores_df).rename(
        columns={
            "pearson_r": "full_data_r",
            "pval": "full_data_pval",
            "n_samples": "full_data_n",
        }
    )

    table = val_df.merge(full_df, on="score", how="inner")
    table["score_label"] = table["score"].map(SCORE_LABELS)
    table["delta_r"] = table["validation_r"] - table["full_data_r"]
    table = table[
        [
            "score",
            "score_label",
            "full_data_r",
            "full_data_pval",
            "full_data_n",
            "validation_r",
            "validation_pval",
            "validation_n",
            "delta_r",
        ]
    ]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path, index=False)
    print(f"Saved comparison table → {out_path}")
    print()
    print(table.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
