#!/bin/bash
#SBATCH --job-name=combine_coverage_decontam_99
#SBATCH --error=./err/13.grpkm/PD_get_RPKM_decontam_99.err
#SBATCH --output=./out/13.grpkm/PD_get_RPKM_decontam_99.out
#SBATCH --mem=128gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# --- LOAD MODULES ---
module purge
module load R

Rscript 02.get_RPKM_decontam_99.R
