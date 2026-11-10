#!/usr/bin/env python3
"""
Compute beta-weighted KEGG pathway scores and TMB correlations for all
pathways shown in fig1b_kegg_ora_all_tmb_genes.Rmd (filtered kegg_ora_tmb_up).

Mirrors local/exploratory/32_calculate_non_pqc_pathway_scores.ipynb and
local/figures/fig1d_pathway_corr_bs_types.ipynb.
"""

import argparse
import re
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

EXCLUDE_SUBCATEGORIES = {
    "Global and overview maps",
    "Transport and catabolism",
    "Endocrine system",
}
EXCLUDE_DESCRIPTIONS = {
    "Oocyte meiosis",
    "Progesterone-mediated oocyte maturation",
}
METABOLISM_SUBCATEGORIES = {
    "Energy metabolism",
    "Carbohydrate metabolism",
    "Metabolism of terpenoids and polyketides",
    "Glycan biosynthesis and metabolism",
    "Amino acid metabolism",
}

# Keep column names aligned with complex_axis_score_non_pqc.csv where applicable.
LEGACY_SLUG = {
    "Ribosome biogenesis in eukaryotes": "ribosome_biogenesis_in_eukaryotes",
    "Spliceosome": "spliceosome",
    "DNA replication": "dna_replication",
    "Cell cycle": "cell_cycle",
    "Oxidative phosphorylation": "oxidative_phosphorylation",
    "Antigen processing and presentation": "antigen_processing_and_presentation",
}


def get_args():
    parser = argparse.ArgumentParser(
        description="Score fig1b KEGG ORA pathways and correlate with TMB."
    )
    parser.add_argument("--expr", required=True)
    parser.add_argument("--tmb_weights", required=True)
    parser.add_argument("--kegg_ora", required=True)
    parser.add_argument("--output_scores", required=True)
    parser.add_argument(
        "--output_corr_cancer_type",
        required=True,
        help="Correlations with cancer-type bootstrap CIs (fig1d style).",
    )
    parser.add_argument(
        "--output_corr_within_cancer",
        required=True,
        help="Correlations with within-cancer-type sample bootstrap CIs (fig1c style).",
    )
    parser.add_argument("--n_boot", type=int, default=1000)
    parser.add_argument(
        "--min_samples_per_cancer",
        type=int,
        default=5,
        help="Min samples per cancer type for within-cancer bootstrap.",
    )
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--test_rows", type=int, default=None)
    return parser.parse_args()


def slugify(description: str) -> str:
    if description in LEGACY_SLUG:
        return LEGACY_SLUG[description]
    slug = re.sub(r"[^a-z0-9]+", "_", description.lower()).strip("_")
    return slug


def assign_group(subcategory: str) -> str:
    if subcategory == "Replication and repair":
        return "DNA replication and repair"
    if subcategory in METABOLISM_SUBCATEGORIES:
        return "Metabolism"
    if subcategory == "Immune system":
        return "Immune"
    return subcategory


def filter_fig1b_pathways(df: pd.DataFrame) -> pd.DataFrame:
    out = df[
        (df["category"] != "Human Diseases")
        & (~df["subcategory"].isin(EXCLUDE_SUBCATEGORIES))
        & (~df["Description"].isin(EXCLUDE_DESCRIPTIONS))
    ].copy()
    out["pathway_slug"] = out["Description"].map(slugify)
    out["group"] = out["subcategory"].map(assign_group)
    return out.sort_values("neglogP", ascending=False).reset_index(drop=True)


def pathway_gene_dict(df: pd.DataFrame) -> dict:
    pathways = {}
    meta = {}
    for _, row in df.iterrows():
        slug = row["pathway_slug"]
        genes = sorted(set(row["geneID"].split("/")))
        pathways[slug] = genes
        meta[slug] = {
            "kegg_id": row["ID"],
            "description": row["Description"],
            "subcategory": row["subcategory"],
            "group": row["group"],
            "n_genes_in_pathway": len(genes),
        }
    return pathways, meta


def compute_pathway_scores(expr, beta_df, pathway_to_genes):
    axis_genes = set(beta_df.index)
    expr_genes = set(expr.columns)
    axis_all = sorted(axis_genes & expr_genes)
    if not axis_all:
        raise ValueError("No axis genes overlap expression matrix.")

    X = expr[axis_all]
    X = X[beta_df.index]
    beta_all = beta_df.loc[X.columns, "estimate"]

    scores = {}
    n_scored = {}
    for pathway, genes in pathway_to_genes.items():
        genes_in_data = sorted(set(genes) & set(X.columns))
        n_scored[pathway] = len(genes_in_data)
        if not genes_in_data:
            print(f"  {pathway}: 0 genes overlap — skipped")
            continue

        beta_pathway = beta_all.copy()
        beta_pathway.loc[~beta_pathway.index.isin(genes_in_data)] = 0
        scores[pathway] = X.to_numpy() @ beta_pathway.to_numpy()
        print(f"  {pathway}: {len(genes_in_data)} genes")

    return pd.DataFrame(scores, index=expr.index), n_scored


def _empirical_p_boot(boot_rs: np.ndarray) -> float:
    n = len(boot_rs)
    return 2 * min(
        (np.sum(boot_rs <= 0) + 1) / (n + 1),
        (np.sum(boot_rs > 0) + 1) / (n + 1),
    )


def _pathway_corr_row(pathway, df_p, boot_rs):
    boot_rs = np.array(boot_rs)
    return {
        "pathway": pathway,
        "description": df_p["description"].iloc[0],
        "group": df_p["group"].iloc[0],
        "kegg_id": df_p["kegg_id"].iloc[0],
        "r": pearsonr(df_p["score"], df_p["log_mut"])[0],
        "ci_low": np.percentile(boot_rs, 2.5),
        "ci_high": np.percentile(boot_rs, 97.5),
        "p_boot": _empirical_p_boot(boot_rs),
        "n_samples": len(df_p),
        "n_genes_scored": df_p["n_genes_scored"].iloc[0],
        "n_genes_in_pathway": df_p["n_genes_in_pathway"].iloc[0],
    }


def bootstrap_corr_by_cancer_type(df_long, n_boot=1000, seed=2026):
    """Resample cancer types with replacement, then pool samples (fig1d)."""
    np.random.seed(seed)
    rows = []

    for pathway, df_p in df_long.groupby("pathway"):
        df_p = df_p.dropna(subset=["score", "log_mut"])
        if len(df_p) < 3:
            continue

        cancer_types = df_p["cancer_type"].unique()
        n_cancer = len(cancer_types)
        groups = {ct: df_p[df_p["cancer_type"] == ct] for ct in cancer_types}

        boot_rs = []
        for _ in range(n_boot):
            sampled_types = np.random.choice(cancer_types, size=n_cancer, replace=True)
            df_boot = pd.concat([groups[ct] for ct in sampled_types])
            try:
                boot_rs.append(pearsonr(df_boot["score"], df_boot["log_mut"])[0])
            except Exception:
                continue

        if boot_rs:
            row = _pathway_corr_row(pathway, df_p, boot_rs)
            row["bootstrap_method"] = "cancer_type"
            rows.append(row)

    corr_df = pd.DataFrame(rows)
    if not corr_df.empty:
        _, corr_df["fdr"] = fdrcorrection(corr_df["p_boot"])
        corr_df = corr_df.sort_values("r", ascending=False)
    return corr_df


def bootstrap_corr_within_cancer(df_long, n_boot=1000, seed=2026, min_samples_per_cancer=5):
    """Resample samples within each cancer type, then pool (fig1c)."""
    np.random.seed(seed)
    rows = []

    for pathway, df_p in df_long.groupby("pathway"):
        df_p = df_p.dropna(subset=["score", "log_mut"])
        if len(df_p) < 3:
            continue

        boot_rs = []
        for _ in range(n_boot):
            resampled = []
            for _, df_sub in df_p.groupby("cancer_type"):
                if len(df_sub) < min_samples_per_cancer:
                    continue
                resampled.append(df_sub.sample(n=len(df_sub), replace=True))
            if not resampled:
                continue
            df_boot = pd.concat(resampled, axis=0)
            try:
                boot_rs.append(pearsonr(df_boot["score"], df_boot["log_mut"])[0])
            except Exception:
                continue

        if boot_rs:
            row = _pathway_corr_row(pathway, df_p, boot_rs)
            row["bootstrap_method"] = "within_cancer"
            rows.append(row)

    corr_df = pd.DataFrame(rows)
    if not corr_df.empty:
        _, corr_df["fdr"] = fdrcorrection(corr_df["p_boot"])
        corr_df = corr_df.sort_values("r", ascending=False)
    return corr_df


def main():
    args = get_args()

    expr = pd.read_csv(args.expr, nrows=args.test_rows)
    tmb_weights = pd.read_csv(args.tmb_weights).drop_duplicates(subset="gene")
    kegg = filter_fig1b_pathways(pd.read_csv(args.kegg_ora))

    print(f"{len(kegg)} pathways after fig1b filters")
    pathway_to_genes, meta = pathway_gene_dict(kegg)

    beta_df = tmb_weights.set_index("gene")
    gene_cols = [c for c in expr.columns if c not in METADATA_COLS]

    print("Computing pathway scores...")
    score_mat, n_scored = compute_pathway_scores(expr[gene_cols], beta_df, pathway_to_genes)

    scores_df = score_mat.copy()
    for col in ["barcode", "cancer_type", "log_mut", "purity"]:
        if col in expr.columns:
            scores_df[col] = expr[col].values

    out_scores = Path(args.output_scores)
    out_scores.parent.mkdir(parents=True, exist_ok=True)
    scores_df.to_csv(out_scores, index=False)
    print(f"Saved pathway scores → {out_scores}")

    meta_df = pd.DataFrame.from_dict(meta, orient="index").reset_index(names="pathway")
    meta_df["n_genes_scored"] = meta_df["pathway"].map(n_scored)

    pathway_cols = [c for c in scores_df.columns if c not in METADATA_COLS]
    id_cols = [c for c in ["barcode", "cancer_type", "log_mut", "purity"] if c in scores_df.columns]
    df_long = scores_df.melt(
        id_vars=id_cols,
        value_vars=pathway_cols,
        var_name="pathway",
        value_name="score",
    ).merge(meta_df, on="pathway", how="left")

    print(f"Computing cancer-type bootstrap correlations (n_boot={args.n_boot})...")
    corr_cancer_type = bootstrap_corr_by_cancer_type(
        df_long, n_boot=args.n_boot, seed=args.seed
    )
    out_ct = Path(args.output_corr_cancer_type)
    corr_cancer_type.to_csv(out_ct, index=False)
    print(f"Saved cancer-type bootstrap correlations → {out_ct}")

    print(f"Computing within-cancer bootstrap correlations (n_boot={args.n_boot})...")
    corr_within = bootstrap_corr_within_cancer(
        df_long,
        n_boot=args.n_boot,
        seed=args.seed,
        min_samples_per_cancer=args.min_samples_per_cancer,
    )
    out_within = Path(args.output_corr_within_cancer)
    corr_within.to_csv(out_within, index=False)
    print(f"Saved within-cancer bootstrap correlations → {out_within}")


if __name__ == "__main__":
    main()
