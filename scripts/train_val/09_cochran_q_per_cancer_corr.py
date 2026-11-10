#!/usr/bin/env python3
"""
Meta-analytic Cochran's Q for heterogeneity of per-cancer correlations.

For each score type, tests whether Pearson r with TMB differs across cancer types
(H0: one common true correlation; Ha: heterogeneity across cancers).

Uses Fisher z transformation with sample-size-based weights.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

SCORE_COLS = ["TMB_axis_score", "pqc_axis_score", "non_pqc_axis_score"]
SCORE_LABELS = {
    "TMB_axis_score": "TMB axis",
    "pqc_axis_score": "PQC axis",
    "non_pqc_axis_score": "Non-PQC axis",
}
EXCLUDE_CANCERS = ["THYM", "UVM"]


def get_args():
    parser = argparse.ArgumentParser(
        description="Cochran's Q heterogeneity test for per-cancer correlations."
    )
    parser.add_argument("--val_per_cancer", required=True)
    parser.add_argument("--full_per_cancer", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def load_long_per_cancer(path):
    df = pd.read_csv(path)
    r_col = "pearson_r" if "pearson_r" in df.columns else "r"
    n_col = "n_samples" if "n_samples" in df.columns else "n"
    if r_col not in df.columns or n_col not in df.columns:
        raise ValueError(f"{path} must contain correlation and sample-size columns.")

    df = df.loc[
        df["cancer_type"] != "ALL",
        ["cancer_type", "score", r_col, n_col],
    ].copy()
    return df.rename(columns={r_col: "pearson_r", n_col: "n_samples"})


def cochran_q_heterogeneity(r, n):
    """Cochran's Q on Fisher-z transformed per-cancer correlations."""
    r = np.clip(np.asarray(r, dtype=float), -0.999999, 0.999999)
    n = np.asarray(n, dtype=float)
    valid = n > 3
    r = r[valid]
    n = n[valid]

    if len(r) < 2:
        return None

    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    weights = 1.0 / (se**2)
    z_bar = np.sum(weights * z) / np.sum(weights)
    q_stat = np.sum(weights * (z - z_bar) ** 2)
    df = len(r) - 1
    p_value = chi2.sf(q_stat, df)
    i2 = max(0.0, (q_stat - df) / q_stat * 100) if q_stat > 0 else 0.0
    pooled_r = np.tanh(z_bar)

    return {
        "n_cancers": len(r),
        "cochran_q": q_stat,
        "p_value": p_value,
        "df": df,
        "i2_percent": i2,
        "pooled_r": pooled_r,
        "median_r": float(np.median(r)),
        "min_r": float(np.min(r)),
        "max_r": float(np.max(r)),
    }


def run_for_subset(long_df, cancers, dataset, subset_name):
    rows = []
    sub = long_df.loc[long_df["cancer_type"].isin(cancers)]
    print(f"\n{dataset} | {subset_name} | up to {len(cancers)} cancers")

    for score in SCORE_COLS:
        score_df = sub.loc[sub["score"] == score].copy()
        result = cochran_q_heterogeneity(
            score_df["pearson_r"].values,
            score_df["n_samples"].values,
        )
        if result is None:
            continue

        row = {
            "dataset": dataset,
            "cancer_subset": subset_name,
            "score": score,
            "score_label": SCORE_LABELS[score],
            **result,
        }
        rows.append(row)
        print(
            f"  {SCORE_LABELS[score]}: Q={result['cochran_q']:.2f}, "
            f"p={result['p_value']:.2e}, I²={result['i2_percent']:.1f}%, "
            f"pooled r={result['pooled_r']:.3f}, "
            f"median r={result['median_r']:.3f} "
            f"[{result['min_r']:.3f}, {result['max_r']:.3f}]"
        )
    return rows


def main():
    args = get_args()

    full_long = load_long_per_cancer(args.full_per_cancer)
    val_long = load_long_per_cancer(args.val_per_cancer)

    all_cancers = sorted(set(full_long["cancer_type"]) & set(val_long["cancer_type"]))
    filtered_cancers = [c for c in all_cancers if c not in EXCLUDE_CANCERS]

    subsets = [
        ("all_cancers", all_cancers),
        ("exclude_THYM_UVM", filtered_cancers),
    ]

    rows = []
    for subset_name, cancers in subsets:
        rows.extend(run_for_subset(full_long, cancers, "full_data", subset_name))
        rows.extend(run_for_subset(val_long, cancers, "validation", subset_name))

    out = pd.DataFrame(rows)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
