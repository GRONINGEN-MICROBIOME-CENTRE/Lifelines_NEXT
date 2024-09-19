#!/bin/bash
#SBATCH --job-name=combine_coverage
#SBATCH --error=./err/10.cco/PD_combine_cov.err
#SBATCH --output=./out/10.cco/PD_combine_cov.out
#SBATCH --mem=128gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# --- LOAD MODULES ---
module purge
module load R

Rscript combine_coverage.R /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/mapping/VLP_to_w_neg_der95_NCP/coverage_NCP/ *_NCP.coverage.txt
