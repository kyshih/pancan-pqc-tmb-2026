"""
Compute sample-level transcriptional axis scores

This script computes two transcriptional module scores for each sample:

1. TMB Axis Score
2. PQC Axis Score

These scores summarize the activity of biologically meaningful gene
expression programs using linear projections or module scoring.

------------------------------------------------------------
TMB Axis
------------------------------------------------------------

The TMB axis represents the transcriptional program associated with
increasing tumor mutational burden (TMB).

Axis construction:
    A generalized linear mixed model (GLMM) was fit across TCGA tumors:

        expression_gene ~ log10(TMB) + tumor_purity + (1 | cancer_type)

Genes significantly associated with TMB (FDR < 0.05, β > 0) were retained.
The GLMM coefficients define the direction of the TMB transcriptional axis.

Gene expression is first standardized across samples:

        z_i = (expression_i - mean_i) / sd_i
        
Sample scoring:
    For each sample, the TMB axis score is computed as a weighted sum:

        TMB_axis_score = Σ (β_i * expression_i)

This can be interpreted as projecting each sample's expression profile
onto the TMB-associated transcriptional axis.


------------------------------------------------------------
PQC Axis
------------------------------------------------------------

The PQC (protein quality control) axis summarizes the expression of
proteostasis-related pathways, including ribosome, proteasome,
UPR, and heat shock response genes.

Gene expression is first standardized across samples:

        z_i = (expression_i - mean_i) / sd_i

The PQC axis score is then computed as the average standardized
expression across PQC genes:

        PQC_axis_score = mean(z_1, z_2, ..., z_n)

This score reflects the overall expression of the PQC transcriptional program.
"""

#TODO: compute axis doesnt work because expr reindex

import pandas as pd
import argparse
from pathlib import Path
import numpy as np
    
# ----------------------------
# Argument parser
# ----------------------------
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Compute sample-level TMB axis and PQC axis scores.")
    parser.add_argument("--expr", required=True, help="Path to regression-ready TCGA expression table (tcga_regression_ready_imputed.csv)")
    parser.add_argument("--tmb_weights", required=True, help="Path to GLMM beta estimates for TMB axis genes (tcga_tmb_genes_posBeta_fdr0.05.csv)")
    parser.add_argument("--pqc_genes", required=True, help="Path to PQC gene list")
    parser.add_argument("--out", required=True, help="Output file path for sample axis scores")
    parser.add_argument("--test_rows", type=int, default=None, help="Read only first N rows (for testing)")
    
    return parser.parse_args()

def main():
    args = get_args()
    
    expr_path = Path(args.expr)
    tmb_weights_path = Path(args.tmb_weights)
    pqc_genes_path = Path(args.pqc_genes)
    out_path = Path(args.out)

    # create output directory
    out_path.mkdir(exist_ok=True, parents=True)
    
    # ------------------------
    # Load data
    # ------------------------
    # load expression data
    expr = pd.read_csv(expr_path)
    print(f"Loaded tcga cancer sample expression matrix: {expr.shape[0]} samples × {expr.shape[1]} features")

    # load GLMM beta estimate for tmb axis genes
    tmb_weights = pd.read_csv(tmb_weights_path)
    print(f"Loaded GLMM beta estimates: {tmb_weights.shape[0]} genes")
    
    # load pqc genes
    pqc_genes = pd.read_csv(pqc_genes_path, sep="\t")
    print(f"Loaded PQC gene list: {pqc_genes.shape[0]} genes")
    
    #------------------------
    # calculate TMB and pqc axis score
    #------------------------
    beta = tmb_weights.set_index("gene")["estimate"]

    pqc_genes.loc[pqc_genes["Hugo"]=="TOR4A", "Hugo"] = "C9orf167" # TOR4A = C9orf167 in tcga
    pqc_genes.loc[pqc_genes["Hugo"]=="FAXC", "Hugo"] = "C6orf168" # FAXC = C6orf168 in tcga
    
    expr["TMB_axis_score"] = compute_axis(expr, beta.index, beta)

    expr["pqc_axis_score"] = compute_axis(expr, pqc_genes["Hugo"])
    
    # ------------------------
    # Save results
    # ------------------------
    expr_axis = expr[["barcode", "cancer_type", "log_mut", "purity", "TMB_axis_score", "pqc_axis_score"]].copy()
    out_file = out_path / "axis_sample_level.parquet"
    
    expr_axis.to_parquet(out_file, index=False)
    print(f"Saved axis scores → {out_file}")

def compute_axis(
    expr: pd.DataFrame,
    genes,
    weights = None,
    zscore = True,
):
    """
    Compute a gene program axis score from expression data.

    Parameters
    ----------
    expr : pd.DataFrame
        Expression matrix (samples x genes).
    genes : list-like
        Gene symbols defining the axis.
    weights : pd.Series, optional
        Gene weights (indexed by gene). If None, computes simple mean.
    zscore : bool
        Whether to z-score genes across samples before computing axis.

    Returns
    -------
    pd.Series
        Axis score per sample.
    """

    # ensure genes exist in expression matrix
    genes = expr.columns.intersection(genes) # preserves the order of the left index.
    print(f"{len(genes)} genes used for axis")
    
    if len(genes) == 0:
        raise ValueError("No genes overlap with expression matrix.")

    X = expr[genes]

    # z-score per gene if requested
    if zscore:
        X = (X - X.mean()) / X.std()

    # weighted axis
    if weights is not None:
        weights = weights.reindex(genes)
        score = X @ weights

    # unweighted axis
    else:
        score = X.mean(axis=1)

    return score

if __name__ == "__main__":
    main()
    

    