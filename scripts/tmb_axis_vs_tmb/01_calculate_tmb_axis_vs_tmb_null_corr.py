import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

def get_args():
    parser = argparse.ArgumentParser(description="Null TMB-axis vs TMB correlation analysis")
    parser.add_argument("--expr",        required=True, help="Path to regression-ready expression table (samples × genes + log_mut column)")
    parser.add_argument("--tmb_weights", required=True, help="Path to GLMM beta estimates (columns: gene, estimate)")
    parser.add_argument("--out",         required=True, help="Output directory")
    return parser.parse_args()

def main():
    args = get_args()
    out  = Path(args.out)

    # ── Load inputs ───────────────────────────────────────────────────────────
    expr        = pd.read_csv(args.expr)
    tmb_weights = pd.read_csv(args.tmb_weights, sep="\t")
    
    tmb_weights = tmb_weights.drop_duplicates(subset="gene")
    axis_genes = tmb_weights[(tmb_weights["estimate"]>0)&(tmb_weights["fdr_bh"]<0.05)]
    print(f"{len(axis_genes)} tmb axis genes.")
    
    axis_genes = set(axis_genes["gene"])
    expr_genes = set(expr.columns)
    META_COLS  = {"log_mut", "cancer_type", "purity", "barcode"}
    
    axis_genes = list(axis_genes & expr_genes)
    
    log10_tmb = expr["log_mut"]
    betas     = tmb_weights.set_index("gene")["estimate"]
    
    # ── Setup ─────────────────────────────────────────────────────────────────
    nonsig_genes = [g for g in expr.columns if g not in axis_genes and g not in META_COLS]
    nonsig_betas = betas.reindex(nonsig_genes).dropna()   # exclude genes absent from GLMM
    nonsig_genes = list(nonsig_betas.index)
    expr_nonsig  = expr[nonsig_genes]
    print(f"{len(nonsig_genes)} genes in null pool ({len(expr.columns) - len(axis_genes) - len(META_COLS) - len(nonsig_genes)} unmodeled genes excluded)")
    tmb_values   = log10_tmb.values   # fixed order matches expr rows

    gene_set_sizes = [5, 10, 50, 100, 500, 1000, 2000, 4903]
    n_iter         = 2
    rng            = np.random.default_rng(seed=2026)

    # ── Null sampling ─────────────────────────────────────────────────────────
    for size in gene_set_sizes:
        if size > len(nonsig_genes):
            warnings.warn(
                f"Gene-set size {size} exceeds nonsig pool ({len(nonsig_genes)} genes); skipping."
            )
            continue

        size_dir = out / "random_gene_sets" / f"n_{size}"
        size_dir.mkdir(parents=True, exist_ok=True)

        corr_records = []
        gene_records = []

        for i in range(n_iter):
            sampled = rng.choice(nonsig_genes, size=size, replace=False)
            
            w     = nonsig_betas.reindex(sampled).values         # (size,)
            X     = expr_nonsig.reindex(columns=sampled).values  # (n_samples, size)
            score = X @ w                                         # (n_samples,)

            r, _ = stats.pearsonr(score, tmb_values)
            corr_records.append({"iter": i, "r": r})
            for g in sampled:
                gene_records.append({"iter": i, "gene": g})

        pd.DataFrame(corr_records, columns=["iter", "r"]).to_csv(
            size_dir / "correlations.tsv", sep="\t", index=False
        )
        pd.DataFrame(gene_records, columns=["iter", "gene"]).to_csv(
            size_dir / "sampled_gene_sets.tsv", sep="\t", index=False
        )
        print(f"n={size}: wrote {size_dir}")

if __name__ == "__main__":
    main()
