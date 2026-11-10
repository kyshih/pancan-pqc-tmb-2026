#!/usr/bin/env python3
"""
Write a markdown summary report for a train/validation run.
"""

import argparse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCORE_LABELS = {
    "TMB_axis_score": "TMB axis",
    "pqc_axis_score": "PQC axis",
    "non_pqc_axis_score": "Non-PQC axis",
}


def get_args():
    parser = argparse.ArgumentParser(description="Generate train/val markdown report.")
    parser.add_argument("--mode", required=True, choices=["fixed_genes", "train_selected"])
    parser.add_argument("--output_dir", required=True, help="Run output directory.")
    parser.add_argument("--split", required=True, help="Train/val split CSV.")
    parser.add_argument("--beta", required=True, help="Beta TSV used for scoring.")
    parser.add_argument(
        "--corr_all",
        default=None,
        help="Per-cancer correlation CSV for all score types.",
    )
    parser.add_argument(
        "--corr_pancan",
        default=None,
        help="Pan-cancer correlation CSV for all score types.",
    )
    parser.add_argument(
        "--compare",
        default=None,
        help="Optional full-data vs validation pan-cancer comparison CSV.",
    )
    parser.add_argument(
        "--sign_test",
        default=None,
        help="Optional per-cancer sign-test summary CSV.",
    )
    parser.add_argument(
        "--cochran_q",
        default=None,
        help="Optional Cochran's Q summary CSV.",
    )
    parser.add_argument(
        "--direction_summary",
        default=None,
        help="Optional per-cancer median/IQR + sign test summary CSV.",
    )
    parser.add_argument("--scores", required=True, help="Validation axis scores CSV.")
    parser.add_argument(
        "--axis_genes",
        default=None,
        help="Optional axis gene list TSV (train_selected mode).",
    )
    parser.add_argument(
        "--regression_with_fdr",
        default=None,
        help="Optional regression results with fdr_bh (train_selected mode).",
    )
    parser.add_argument(
        "--reference_corr",
        default=None,
        help="Optional reference correlation CSV (e.g. full-data or LOCO) for comparison.",
    )
    parser.add_argument(
        "--report",
        required=True,
        help="Output markdown report path.",
    )
    return parser.parse_args()


def fmt_r(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.4f}"


def score_rows(corr_all_df, score_col):
    rows = corr_all_df.loc[corr_all_df["score"] == score_col]
    overall = rows.loc[rows["cancer_type"] == "ALL"].iloc[0]
    per_cancer = rows.loc[rows["cancer_type"] != "ALL"].copy()
    return overall, per_cancer


def main():
    args = get_args()
    out_dir = Path(args.output_dir)

    split_df = pd.read_csv(args.split)
    beta_df = pd.read_csv(args.beta, sep="\t")
    scores_df = pd.read_csv(args.scores)

    corr_all_path = (
        Path(args.corr_all)
        if args.corr_all
        else out_dir / "corr_val_axis_scores_vs_tmb.csv"
    )
    corr_pancan_path = (
        Path(args.corr_pancan)
        if args.corr_pancan
        else out_dir / "corr_val_pancan_axis_scores_vs_tmb.csv"
    )

    if not corr_all_path.exists():
        raise FileNotFoundError(f"Missing per-cancer correlation file: {corr_all_path}")
    if not corr_pancan_path.exists():
        raise FileNotFoundError(f"Missing pan-cancer correlation file: {corr_pancan_path}")

    corr_all_df = pd.read_csv(corr_all_path)
    corr_pancan_df = pd.read_csv(corr_pancan_path)

    compare_path = (
        Path(args.compare)
        if args.compare
        else out_dir / "pancan_corr_full_vs_validation.csv"
    )
    compare_df = pd.read_csv(compare_path) if compare_path.exists() else None

    sign_test_path = (
        Path(args.sign_test)
        if args.sign_test
        else out_dir / "sign_test_per_cancer_corr_conservation.csv"
    )
    sign_test_df = pd.read_csv(sign_test_path) if sign_test_path.exists() else None

    cochran_path = (
        Path(args.cochran_q) if args.cochran_q else out_dir / "cochran_q_per_cancer_corr.csv"
    )
    cochran_df = pd.read_csv(cochran_path) if cochran_path.exists() else None

    direction_path = (
        Path(args.direction_summary)
        if args.direction_summary
        else out_dir / "per_cancer_direction_summary.csv"
    )
    direction_df = pd.read_csv(direction_path) if direction_path.exists() else None

    tmb_overall, tmb_per_cancer = score_rows(corr_all_df, "TMB_axis_score")
    median_r = tmb_per_cancer["pearson_r"].median()

    n_train = (split_df["split"] == "train").sum()
    n_val = (split_df["split"] == "val").sum()
    n_cancers = split_df["cancer_type"].nunique()

    mode_title = {
        "fixed_genes": "Fixed gene list (refit betas only)",
        "train_selected": "Train-only gene selection (no leakage)",
    }[args.mode]

    lines = [
        "# Train/Validation TMB Axis Report",
        "",
        f"- **Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"- **Mode:** `{args.mode}` — {mode_title}",
        f"- **Output directory:** `{out_dir}`",
        "",
        "## What this run tests",
        "",
    ]

    if args.mode == "fixed_genes":
        lines += [
            "Uses the **full-data posBeta gene list** (FDR selected on all samples),",
            "but refits GLMM weights on 70% training samples only.",
            "Scores and correlates on held-out 30% validation samples.",
            "",
            "> **Note:** Gene selection still used validation samples indirectly.",
            "> For a fully held-out validation, run with `--select_genes_on_train`.",
            "",
        ]
    else:
        lines += [
            "Refits GLMM on **all genes** using training samples only,",
            "applies BH FDR on train p-values, selects posBeta genes (`estimate > 0`, `fdr_bh < 0.05`),",
            "then scores and correlates on held-out validation samples.",
            "",
            "> **Note:** This is the principled no-leakage validation mode.",
            "",
        ]

    lines += [
        "## Split summary",
        "",
        f"| Metric | Value |",
        f"|--------|------:|",
        f"| Cancer types | {n_cancers} |",
        f"| Train samples | {n_train} |",
        f"| Validation samples | {n_val} |",
        f"| Train fraction | {n_train / (n_train + n_val):.1%} |",
        f"| Axis genes used for scoring | {len(beta_df)} |",
        f"| Validation samples scored | {len(scores_df)} |",
        "",
        "## Validation correlation (TMB axis score vs log10 TMB)",
        "",
        f"| Metric | Value |",
        f"|--------|------:|",
        f"| Overall Pearson r (ALL) | {fmt_r(tmb_overall['pearson_r'])} |",
        f"| Overall p-value | {tmb_overall['pval']:.2e} |",
        f"| Median per-cancer r | {fmt_r(median_r)} |",
        f"| Cancers with r > 0 | {(tmb_per_cancer['pearson_r'] > 0).sum()} / {len(tmb_per_cancer)} |",
        "",
        "### Top 10 cancer types by validation r",
        "",
        "| cancer_type | r | p-value | n_val |",
        "|-------------|---:|--------:|------:|",
    ]

    top = tmb_per_cancer.sort_values("pearson_r", ascending=False).head(10)
    for _, row in top.iterrows():
        lines.append(
            f"| {row['cancer_type']} | {fmt_r(row['pearson_r'])} "
            f"| {row['pval']:.2e} | {int(row['n_samples'])} |"
        )

    lines += [
        "",
        "## Pan-cancer validation correlation (all cancers pooled)",
        "",
        "| score | overall r | p-value | n_val |",
        "|-------|----------:|--------:|------:|",
    ]
    for score_col, label in SCORE_LABELS.items():
        row = corr_pancan_df.loc[corr_pancan_df["score"] == score_col].iloc[0]
        lines.append(
            f"| {label} | {fmt_r(row['pearson_r'])} "
            f"| {row['pval']:.2e} | {int(row['n_samples'])} |"
        )
    lines.append(f"\nPan-cancer table: `{corr_pancan_path}`")

    if compare_df is not None:
        lines += [
            "",
            "## Pan-cancer r: full data vs validation",
            "",
            "| score | full-data r | validation r | delta (val - full) |",
            "|-------|------------:|-------------:|-------------------:|",
        ]
        for _, row in compare_df.iterrows():
            lines.append(
                f"| {row['score_label']} | {fmt_r(row['full_data_r'])} "
                f"| {fmt_r(row['validation_r'])} | {fmt_r(row['delta_r'])} |"
            )
        lines.append(f"\nComparison table: `{compare_path}`")

    if sign_test_df is not None:
        lines += [
            "",
            "## Per-cancer sign concordance (full data vs validation)",
            "",
            "Binomial sign test: H₀ = probability of same sign is 0.5; "
            "Hₐ = same sign more often than expected by chance.",
            "",
            "| score | same sign (of 32) | p-value | both r > 0 (of 32) | val r > 0 given full r > 0 |",
            "|-------|------------------:|--------:|-------------------:|----------------------------:|",
        ]
        for _, row in sign_test_df.iterrows():
            conditional = (
                f"{int(row['n_full_pos_and_val_pos'])}/{int(row['n_full_data_positive'])}"
                if row["n_full_data_positive"] > 0
                else "NA"
            )
            lines.append(
                f"| {row['score_label']} | "
                f"{int(row['n_same_sign'])}/{int(row['n_cancers'])} | "
                f"{row['p_sign_concordance']:.2e} | "
                f"{int(row['n_both_positive'])}/{int(row['n_cancers'])} | "
                f"{conditional} |"
            )
        lines += [
            "",
            "Notes:",
            "- **Same sign (32):** full-data and validation Pearson r have the same direction.",
            "- **Val r > 0 given full r > 0:** conditional on cancers with positive full-data r only "
            "(excludes cancers already ≤ 0 in full data).",
            "",
        ]
        lines.append(f"\nMerged per-cancer table: `{out_dir / 'per_cancer_corr_full_vs_validation.csv'}`")
        lines.append(f"Sign-test summary: `{sign_test_path}`")

    if cochran_df is not None:
        lines += [
            "",
            "## Cochran's Q heterogeneity (per-cancer r across cancer types)",
            "",
            "Meta-analytic Cochran's Q on Fisher-z transformed per-cancer Pearson r. "
            "Tests whether the TMB–score association differs across cancer types "
            "(H₀: one common correlation; significant Q = heterogeneous response).",
            "",
            "| dataset | cancers | score | n | Q | p-value | I² | pooled r | median r |",
            "|---------|---------|-------|--:|--:|--------:|---:|---------:|---------:|",
        ]
        for _, row in cochran_df.iterrows():
            subset = row["cancer_subset"].replace("_", " ")
            lines.append(
                f"| {row['dataset']} | {subset} | {row['score_label']} | "
                f"{int(row['n_cancers'])} | {row['cochran_q']:.2f} | {row['p_value']:.2e} | "
                f"{row['i2_percent']:.1f}% | {fmt_r(row['pooled_r'])} | {fmt_r(row['median_r'])} |"
            )
        lines.append(f"\nCochran's Q table: `{cochran_path}`")

    if direction_df is not None:
        lines += [
            "",
            "## Cross-cancer correlation summary",
            "",
            "**Pan-cancer r** = Pearson r across all samples pooled (fig2c). "
            "**Per-cancer median (IQR)** = distribution of cancer-type-specific r. "
            "**Sign test** = binomial test of r > 0 across cancer types (one-sided).",
            "",
            "| dataset | score | pan-cancer r | median per-cancer r | IQR | n positive | binomial p |",
            "|---------|-------|-------------:|--------------------:|-----|----------:|-----------:|",
        ]
        for _, row in direction_df.iterrows():
            iqr = f"{row['q25_per_cancer_r']:.3f}–{row['q75_per_cancer_r']:.3f}"
            lines.append(
                f"| {row['dataset']} | {row['score_label']} | "
                f"{fmt_r(row['pancan_r'])} | {fmt_r(row['median_per_cancer_r'])} | {iqr} | "
                f"{int(row['n_positive'])}/{int(row['n_cancers'])} | {row['binomial_p']:.2e} |"
            )
        fig_path = (
            out_dir.parent.parent.parent
            / "results/figures/supplementary/per_cancer_r_median_iqr_signtest.pdf"
        )
        lines += [
            "",
            f"Summary table: `{direction_path}`",
            f"Figure: `{fig_path}`",
        ]

    lines += [
        "",
        "## Per-cancer validation correlation by score type (fig2-style)",
        "",
        "| score | overall r | p-value | median per-cancer r |",
        "|-------|----------:|--------:|--------------------:|",
    ]
    for score_col, label in SCORE_LABELS.items():
        overall, per_cancer = score_rows(corr_all_df, score_col)
        lines.append(
            f"| {label} | {fmt_r(overall['pearson_r'])} "
            f"| {overall['pval']:.2e} "
            f"| {fmt_r(per_cancer['pearson_r'].median())} |"
        )
    lines += [
        "",
        f"Full per-cancer table: `{corr_all_path}`",
    ]

    if args.mode == "train_selected" and args.axis_genes:
        axis_genes = pd.read_csv(args.axis_genes, sep="\t")
        lines += [
            "",
            "## Train-only gene selection",
            "",
            f"- Selected axis genes: **{len(axis_genes)}**",
            f"- Selection rule: `estimate > 0` and `fdr_bh < 0.05` on training GLMM",
            f"- Gene list: `{args.axis_genes}`",
        ]

    if args.reference_corr and Path(args.reference_corr).exists():
        ref = pd.read_csv(args.reference_corr)
        ref_overall = ref.loc[ref["cancer_type"] == "ALL"]
        if len(ref_overall):
            ref_r = ref_overall.iloc[0]["pearson_r_tmb_vs_tmb_axis"]
            val_r = tmb_overall["pearson_r"]
            lines += [
                "",
                "## Reference comparison",
                "",
                f"- Reference file: `{args.reference_corr}`",
                f"- Reference overall r: **{fmt_r(ref_r)}**",
                f"- This run overall r: **{fmt_r(val_r)}**",
                f"- Delta (this - reference): **{fmt_r(val_r - ref_r)}**",
            ]

    lines += [
        "",
        "## Output files",
        "",
        f"- Split: `{args.split}`",
        f"- Betas: `{args.beta}`",
        f"- Validation scores: `{args.scores}`",
        f"- Validation per-cancer correlations: `{corr_all_path}`",
        f"- Validation pan-cancer correlations: `{corr_pancan_path}`",
    ]
    if compare_df is not None:
        lines.append(f"- Full vs validation pan-cancer comparison: `{compare_path}`")
    if sign_test_df is not None:
        lines.append(f"- Per-cancer sign concordance test: `{sign_test_path}`")
    if cochran_df is not None:
        lines.append(f"- Cochran's Q heterogeneity: `{cochran_path}`")
    if direction_df is not None:
        lines.append(f"- Per-cancer median/IQR + sign test: `{direction_path}`")

    if args.regression_with_fdr:
        lines.append(f"- Regression + FDR: `{args.regression_with_fdr}`")
    if args.axis_genes:
        lines.append(f"- Axis gene list: `{args.axis_genes}`")

    lines += [
        "",
        "## How to re-run",
        "",
        "```bash",
        "# Fast mode (fixed posBeta genes, refit betas on train):",
        "bash workflows/train_val/01_run_train_val_pipeline.sh",
        "",
        "# Slow mode (all genes, BH selection on train only):",
        "bash workflows/train_val/01_run_train_val_pipeline.sh --select_genes_on_train",
        "```",
        "",
        "See `workflows/train_val/README.md` for full documentation.",
    ]

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n")
    print(f"Saved report → {report_path}")


if __name__ == "__main__":
    main()
