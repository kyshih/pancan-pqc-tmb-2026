#!/usr/bin/env python3
"""Grouped bar chart: pan-cancer r, full data vs validation (fig2 colors)."""

import os
import socket
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

if "scg" in socket.gethostname():
    BASE_DIR = Path("/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected")
else:
    BASE_DIR = Path("~/Documents/research/scg/protein-master/v2_pancan_batch_corrected").expanduser()

TABLE = BASE_DIR / "analysis/03_train_val/fixed_genes/pancan_corr_full_vs_validation.csv"
OUT_DIR = BASE_DIR / "results/figures/supplementary/validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "TMB axis": "#0072B2",
    "PQC axis": "#D55E00",
    "Non-PQC axis": "#009E73",
}

df = pd.read_csv(TABLE)
df = df.set_index("score_label")

labels = ["TMB axis", "PQC axis", "Non-PQC axis"]
x = np.arange(len(labels))
width = 0.34

fig, ax = plt.subplots(figsize=(5.5, 4.2))
full_vals = [df.loc[l, "full_data_r"] for l in labels]
val_vals = [df.loc[l, "validation_r"] for l in labels]

ax.bar(
    x - width / 2,
    full_vals,
    width,
    label="Full data",
    color="#BDBDBD",
    edgecolor="black",
    linewidth=0.6,
)
bars_val = ax.bar(
    x + width / 2,
    val_vals,
    width,
    label="Validation",
    color=[PALETTE[l] for l in labels],
    edgecolor="black",
    linewidth=0.6,
)

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Pan-cancer Pearson r")
ax.set_ylim(0, 0.65)
ax.legend(frameon=False, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, val in zip(bars_val, val_vals):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.015,
        f"{val:.2f}",
        ha="center",
        va="bottom",
        fontsize=9,
    )

fig.tight_layout()
out_pdf = OUT_DIR / "val_pancan_corr_full_vs_validation.pdf"
out_png = OUT_DIR / "val_pancan_corr_full_vs_validation.png"
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved → {out_pdf}")
print(f"Saved → {out_png}")

# LaTeX table for copy-paste
tex_path = BASE_DIR / "results/tables/val_pancan_corr_full_vs_validation.tex"
tex_path.parent.mkdir(parents=True, exist_ok=True)
lines = [
    r"\begin{table}[ht]",
    r"\centering",
    r"\caption{Pan-cancer Pearson correlation between axis scores and $\log_{10}$(TMB). "
    r"Full data ($n=9104$): betas fit on all samples. "
    r"Validation ($n=2746$): betas fit on 70\% training samples; "
    r"correlations on held-out 30\% per cancer type.}",
    r"\label{tab:val_pancan_corr}",
    r"\begin{tabular}{lrr}",
    r"\toprule",
    r"Score & Full data & Validation \\",
    r"\midrule",
]
for _, row in pd.read_csv(TABLE).iterrows():
    lines.append(
        f"{row['score_label']} & {row['full_data_r']:.3f} & {row['validation_r']:.3f} \\\\"
    )
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
tex_path.write_text("\n".join(lines))
print(f"Saved → {tex_path}")
