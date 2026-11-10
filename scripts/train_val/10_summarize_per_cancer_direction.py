#!/usr/bin/env python3
"""
Summarize cross-cancer TMB-axis correlations:

- Pan-cancer r: Pearson r across all samples pooled (matches fig2c).
- Per-cancer median (IQR): distribution of cancer-type-specific r values.
- Binomial sign test: fraction of cancer types with r > 0.
"""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

SCORE_COLS = ["TMB_axis_score", "pqc_axis_score", "non_pqc_axis_score"]
SCORE_LABELS = {
    "TMB_axis_score": "TMB axis",
    "pqc_axis_score": "PQC axis",
    "non_pqc_axis_score": "Non-PQC axis",
}
PALETTE = {
    "TMB_axis_score": "#0072B2",
    "pqc_axis_score": "#D55E00",
    "non_pqc_axis_score": "#009E73",
}


def get_args():
    parser = argparse.ArgumentParser(description="Pan-cancer + per-cancer direction summary.")
    parser.add_argument("--val_per_cancer", required=True)
    parser.add_argument("--full_per_cancer", required=True)
    parser.add_argument(
        "--pancan_full",
        required=True,
        help="Pan-cancer full-data r table (pancan_corr_full_vs_validation.csv).",
    )
    parser.add_argument(
        "--pancan_val",
        required=True,
        help="Pan-cancer validation r table (same file; uses validation_r column).",
    )
    parser.add_argument("--output_table", required=True)
    parser.add_argument("--output_figure", required=True)
    return parser.parse_args()


def load_long_per_cancer(path):
    df = pd.read_csv(path)
    r_col = "pearson_r" if "pearson_r" in df.columns else "r"
    df = df.loc[df["cancer_type"] != "ALL", ["cancer_type", "score", r_col]].copy()
    return df.rename(columns={r_col: "pearson_r"})


def load_pancan(path, dataset):
    df = pd.read_csv(path)
    if dataset == "full_data":
        return df.set_index("score")[["full_data_r", "full_data_n"]].rename(
            columns={"full_data_r": "pancan_r", "full_data_n": "pancan_n"}
        )
    return df.set_index("score")[["validation_r", "validation_n"]].rename(
        columns={"validation_r": "pancan_r", "validation_n": "pancan_n"}
    )


def summarize_per_cancer(values):
    values = pd.Series(values, dtype=float)
    n_pos = int((values > 0).sum())
    n = len(values)
    bt = binomtest(n_pos, n, 0.5, alternative="greater")
    q25 = values.quantile(0.25)
    q75 = values.quantile(0.75)
    return {
        "n_cancers": n,
        "n_positive": n_pos,
        "fraction_positive": n_pos / n,
        "binomial_p": bt.pvalue,
        "median_per_cancer_r": values.median(),
        "q25_per_cancer_r": q25,
        "q75_per_cancer_r": q75,
    }


def build_summary(long_df, pancan_df, dataset):
    rows = []
    for score in SCORE_COLS:
        vals = long_df.loc[long_df["score"] == score, "pearson_r"]
        row = {
            "dataset": dataset,
            "score": score,
            "score_label": SCORE_LABELS[score],
            "pancan_r": pancan_df.loc[score, "pancan_r"],
            "pancan_n": int(pancan_df.loc[score, "pancan_n"]),
            **summarize_per_cancer(vals),
        }
        rows.append(row)
    return rows


def format_p(p):
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def plot_summary(summary_df, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    offsets = {"pancan": -0.15, "per_cancer": 0.15}

    for ax, dataset, title in zip(
        axes,
        ["full_data", "validation"],
        ["Full data", "Validation"],
    ):
        sub = summary_df.loc[summary_df["dataset"] == dataset].set_index("score")
        sub = sub.reindex(SCORE_COLS)
        x = np.arange(len(SCORE_COLS))

        for i, score in enumerate(SCORE_COLS):
            row = sub.loc[score]
            color = PALETTE[score]

            ax.scatter(
                i + offsets["pancan"],
                row["pancan_r"],
                marker="s",
                s=55,
                color=color,
                edgecolor="black",
                linewidth=0.6,
                zorder=4,
                label="Pan-cancer r" if i == 0 else None,
            )
            ax.errorbar(
                i + offsets["per_cancer"],
                row["median_per_cancer_r"],
                yerr=[
                    [row["median_per_cancer_r"] - row["q25_per_cancer_r"]],
                    [row["q75_per_cancer_r"] - row["median_per_cancer_r"]],
                ],
                fmt="o",
                color=color,
                ecolor=color,
                elinewidth=2,
                capsize=4,
                capthick=1.5,
                markersize=7,
                zorder=3,
                label="Per-cancer median (IQR)" if i == 0 else None,
            )
            ax.text(
                i + offsets["per_cancer"],
                row["q75_per_cancer_r"] + 0.035,
                f"{int(row['n_positive'])}/{int(row['n_cancers'])} +",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )
            ax.text(
                i + offsets["per_cancer"],
                row["q75_per_cancer_r"] + 0.075,
                format_p(row["binomial_p"]),
                ha="center",
                va="bottom",
                fontsize=7.5,
            )

        ax.axhline(0, color="0.5", linewidth=0.8, linestyle="--", zorder=1)
        ax.set_xticks(x)
        ax.set_xticklabels([SCORE_LABELS[s] for s in SCORE_COLS])
        ax.set_title(title)
        ax.set_ylim(-0.15, 0.75)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, loc="lower right", fontsize=8)

    axes[0].set_ylabel("Pearson r with TMB")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = get_args()

    full_long = load_long_per_cancer(args.full_per_cancer)
    val_long = load_long_per_cancer(args.val_per_cancer)
    pancan_df = pd.read_csv(args.pancan_full)

    rows = build_summary(full_long, load_pancan(args.pancan_full, "full_data"), "full_data")
    rows += build_summary(val_long, load_pancan(args.pancan_val, "validation"), "validation")
    summary = pd.DataFrame(rows)

    table_path = Path(args.output_table)
    table_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(table_path, index=False)

    plot_summary(summary, args.output_figure)

    print(f"Saved table → {table_path}")
    print(f"Saved figure → {args.output_figure}")
    print()
    print(
        summary[
            [
                "dataset",
                "score_label",
                "pancan_r",
                "median_per_cancer_r",
                "q25_per_cancer_r",
                "q75_per_cancer_r",
                "n_positive",
                "n_cancers",
                "binomial_p",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.4g}")
    )


if __name__ == "__main__":
    main()
