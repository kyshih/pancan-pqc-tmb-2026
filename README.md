# pancan-pqc-tmb

Analysis code for the pan-cancer TCGA analysis of protein quality control (PQC) and tumor mutational burden (TMB) (Shih et al., 2026).

- `scripts/` — pipeline: expression/metadata merging, TMB calculation, per-gene regressions with permutation-based FDR (`shell/run_all.sh`), TMB-axis scoring, leave-one-cancer-out (`loco/`), and train/validation analyses (`train_val/`).
- `local/figures/` — notebooks and R Markdown that generate the main and supplementary figures.

Input data (TCGA / GTEx via UCSC Xena) are not included. Some scripts contain absolute paths that need to be edited to point at your local data.
