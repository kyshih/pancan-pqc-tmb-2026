import pandas as pd
import socket, os

hostname = socket.gethostname()
if "scg" in hostname:
    BASE_DIR = "/labs/mwinslow/Karen/protein-master/v2_pancan_batch_corrected"
else:
    BASE_DIR = os.path.expanduser("~/Documents/research/scg/protein-master/v2_pancan_batch_corrected")

# WORKING_DIR=f"{BASE_DIR}/analysis/01_axis_analysis/axes"
WORKING_DIR=f"{BASE_DIR}/analysis/01_axis_analysis/bh/axes"

os.chdir(WORKING_DIR)
print(f"📂 Working directory set to: {WORKING_DIR}")

def main():
    gencode = pd.read_csv(f"{BASE_DIR}/data/gencode/processed/gencode_v23_gene_id_to_symbol_protein_coding.tsv", sep="\t")

    beta = pd.read_csv(f"{BASE_DIR}/analysis/00_reproduction/tilk_et_al_2023/gene_sets/bh/tcga_tmb_genes_posBeta_fdr0.05.csv")
    # beta = pd.read_csv(f"{BASE_DIR}/analysis/00_reproduction/tilk_et_al_2023/gene_sets/tcga_tmb_genes_negBeta_fdr0.05.csv")
    beta = beta.copy()
    print(f"Loaded GLMM beta estimates: {beta.shape[0]} genes")
    
    pqc_genes = pd.read_csv(f"{BASE_DIR}/data/corum/pqc_annotation.tsv", sep="\t")
    pqc_genes = pqc_genes.copy()
    
    # merge to get ENSG
    mapped_beta = beta.merge(gencode, left_on="gene", right_on="gene_name", how="left").drop(columns=["gene_name"])
    mapped_beta["source_dataset"] = "TCGA"
    dup = mapped_beta.groupby("gene")["gene_id"].nunique()
    multi_map = dup[dup > 1]

    total_genes = len(beta)
    mapped_genes = mapped_beta["gene_id"].notna().sum()
    missing_genes = mapped_beta["gene_id"].isna().sum()

    print(f"Total genes from GLMM: {total_genes}")
    print(f"Mapped to ENSG: {mapped_genes}")
    print(f"Missing genes: {missing_genes}")
    print(f"Multiple ENSG mappings resolved: {len(multi_map)}")
    
    missing = mapped_beta[mapped_beta["gene_id"].isna()]
    # missing.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/axes/qc/missing_tcga_tmb_genes.tsv", index=False)
    missing.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/bh/axes/qc/missing_tcga_tmb_genes.tsv", index=False)
    
    multi_map = mapped_beta[mapped_beta["gene"].isin(multi_map.index)]
    # multi_map.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/axes/qc/multi_mapped_tcga_tmb_genes.tsv", index=False)
    multi_map.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/bh/axes/qc/multi_mapped_tcga_tmb_genes.tsv", index=False)
    
    mapped_beta = mapped_beta[mapped_beta["gene_id"].notna()]
    mapped_beta = mapped_beta.drop_duplicates("gene")
    mapped_beta_slim = mapped_beta[["gene", "gene_id", "estimate", "source_dataset"]]
    print(f"Total genes after removing multi-mapped genes and missing genes: {len(mapped_beta)}")
    
    # pqc_genes = pqc_genes.merge(gencode, left_on="Hugo", right_on="gene_name", how="left").drop(columns=["gene_name"]).rename(columns={"Hugo":"gene"})
    # print("pqc genes mapped to ENSG:", pqc_genes["gene_id"].notna().sum())
    # pqc_genes["dataset"] = "corum"
    # pqc_genes_slim = pqc_genes[["Group", "subgroup", "gene", "gene_id", "dataset"]]
    
    pqc_genes_slim.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/axes/corum_pqc_genes.tsv", sep="\t", index=False)
    mapped_beta_slim.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/bh/axes/tcga_tmb_genes_posBeta.tsv", sep="\t", index=False)
    # mapped_beta_slim.to_csv(f"{BASE_DIR}/analysis/01_axis_analysis/axes/tcga_tmb_genes_negBeta.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()