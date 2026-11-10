#!/usr/bin/env python3
"""
Create a stratified 70/30 train/validation split within each cancer type.

Outputs a CSV with columns: barcode, cancer_type, split (train | val).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(
        description="Create stratified train/validation split within each cancer type."
    )
    parser.add_argument(
        "--expr",
        required=True,
        help="Path to regression-ready expression table (wide format).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV path for split assignments.",
    )
    parser.add_argument(
        "--train_frac",
        type=float,
        default=0.7,
        help="Fraction of samples assigned to training within each cancer type.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible splits.",
    )
    return parser.parse_args()


def main():
    args = get_args()
    rng = np.random.default_rng(args.seed)

    expr = pd.read_csv(args.expr, usecols=["barcode", "cancer_type"])
    expr = expr.drop_duplicates(subset="barcode")

    records = []
    for cancer_type, group in expr.groupby("cancer_type", sort=True):
        barcodes = group["barcode"].to_numpy()
        rng.shuffle(barcodes)

        n_train = int(np.floor(len(barcodes) * args.train_frac))
        n_train = max(n_train, 1)
        if n_train >= len(barcodes):
            n_train = max(len(barcodes) - 1, 1)

        train_barcodes = set(barcodes[:n_train])
        for barcode in barcodes:
            records.append(
                {
                    "barcode": barcode,
                    "cancer_type": cancer_type,
                    "split": "train" if barcode in train_barcodes else "val",
                }
            )

    split_df = pd.DataFrame(records)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    split_df.to_csv(out_path, index=False)

    summary = (
        split_df.groupby(["cancer_type", "split"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={"train": "n_train", "val": "n_val"})
    )
    summary["n_total"] = summary.sum(axis=1)
    summary["train_frac"] = summary["n_train"] / summary["n_total"]

    print(f"Saved split assignments → {out_path}")
    print(f"Total samples: {len(split_df)}")
    print(
        f"Train: {(split_df['split'] == 'train').sum()} | "
        f"Val: {(split_df['split'] == 'val').sum()}"
    )
    print("\nPer-cancer summary:")
    print(summary.to_string())


if __name__ == "__main__":
    main()
