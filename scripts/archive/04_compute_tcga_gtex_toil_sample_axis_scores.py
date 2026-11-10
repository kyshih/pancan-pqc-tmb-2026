# this code works. I just refactored

import pandas as pd
import socket, os

hostname = socket.gethostname()
if "scg" in hostname:
    BASE_DIR = "/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
else:
    BASE_DIR = os.path.expanduser("~/Documents/research/scg/protein-master/v2_pancan_batch_corrected")

WORKING_DIR=f"{BASE_DIR}/data/tcga_gtex_toil"

os.chdir(WORKING_DIR)
print(f"📂 Working directory set to: {WORKING_DIR}")

def main():
    df_gtex = pd.read_parquet(f"{WORKING_DIR}/processed/expression_gtex.parquet")
    df_tcga = pd.read_parquet(f"{WORKING_DIR}/processed/expression_tcga.parquet")
    df_gtex_normal_tissue = df_gtex[(df_gtex["sample_type"]=="Normal Tissue")].copy()
    df_tcga_normal_tissue = df_tcga[(df_tcga["sample_type"]=="Solid Tissue Normal")].copy()

    tmb_axis = pd.read_csv(f"{BASE_DIR}/analysis/01_axis_analysis/axes/tcga_tmb_genes_posBeta.tsv", sep="\t")
    pqc_genes = pd.read_csv(f"{BASE_DIR}/data/corum/processed/corum_pqc_genes.tsv", sep="\t")

    pqc_gene_list = pqc_genes["gene_id"]
    axis_genes = tmb_axis["gene_id"]

    common_genes = list(set(axis_genes) & set(df_gtex_normal_tissue.columns))
    
    # 2. Subset beta table
    beta_df = tmb_axis.set_index("gene_id").loc[common_genes, ["estimate"]]

    # 3. Expression matrix
    X = df_gtex_normal_tissue[common_genes]

    # 4. Z-score genes
    X_z = (X - X.mean()) / X.std()

    # 5. Align explicitly
    beta = beta_df["estimate"].reindex(X_z.columns)

    # 6. Compute score
    tmb_axis_score = X_z @ beta

    common_pqc_genes = df_gtex_normal_tissue.columns.intersection(pqc_gene_list)
    print(f"{len(common_pqc_genes)} PQC genes overlap")

    # expression matrix for PQC genes
    X_pqc = df_gtex_normal_tissue[common_pqc_genes]
    # z-score per gene
    X_pqc_z = (X_pqc - X_pqc.mean()) / X_pqc.std()
    # PQC axis = mean z-score
    pqc_axis = X_pqc_z.mean(axis=1)

    df_axis_gtex_normal_tissue = pd.DataFrame({
    "barcode": df_gtex_normal_tissue["sample_id"],
    "tissue": df_gtex_normal_tissue["tissue"],
    "TMB_axis_score": tmb_axis_score,
    "pqc_axis_score": pqc_axis})
    
    df_axis_gtex_normal_tissue.to_parquet(f"{BASE_DIR}/exploratory/axes/derived/tcga_gtex_toil/gtex_normal_tissue_axis_sample_level.parquet", index=False)

    common_genes = list(set(axis_genes) & set(df_tcga_normal_tissue.columns))

    # 2. Subset beta table
    beta_df = tmb_axis.set_index("gene_id").loc[common_genes, ["estimate"]]

    # 3. Expression matrix
    X = df_tcga_normal_tissue[common_genes]

    # 4. Z-score genes
    X_z = (X - X.mean()) / X.std()

    # 5. Align explicitly
    beta = beta_df["estimate"].reindex(X_z.columns)

    # 6. Compute score
    tmb_axis_score = X_z @ beta

    common_pqc_genes = df_tcga_normal_tissue.columns.intersection(pqc_gene_list)
    print(f"{len(common_pqc_genes)} PQC genes overlap")

    # expression matrix for PQC genes
    X_pqc = df_tcga_normal_tissue[common_pqc_genes]
    # z-score per gene
    X_pqc_z = (X_pqc - X_pqc.mean()) / X_pqc.std()
    # PQC axis = mean z-score
    pqc_axis = X_pqc_z.mean(axis=1)
    
    df_axis_tcga_normal_tissue = pd.DataFrame({
    "barcode": df_tcga_normal_tissue["sample_id"],
    "tissue": df_tcga_normal_tissue["tissue"],
    "TMB_axis_score": tmb_axis_score,
    "pqc_axis_score": pqc_axis})
    
    df_axis_tcga_normal_tissue.to_parquet(f"{BASE_DIR}/exploratory/axes/derived/tcga_gtex_toil/tcga_normal_tissue_axis_sample_level.parquet", index=False)

if __name__ == "__main__":
    main()