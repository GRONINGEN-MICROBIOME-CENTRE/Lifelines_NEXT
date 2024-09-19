#!/bin/bash
#SBATCH --job-name=combine_coverage
#SBATCH --error=./err/11.grpkm/PD_get_RPKM.err
#SBATCH --output=./out/11.grpkm/PD_get_RPKM.out
#SBATCH --mem=128gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# --- LOAD MODULES ---
module purge
module load R

Rscript 01a.get_RPKM.R
