#!/bin/bash
#SBATCH --job-name=dec_AT_generation
#SBATCH --output=./out/13.mdct/decontamination_w_mgs_nc.out
#SBATCH --mem=120gb
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

SAMPLE_DIR=$1 # path to count table

echo ${SAMPLE_DIR}

module load R/4.3.2-gfbf-2023a

# Calculating the alignment rate
Rscript new_algn_rate.R ${SAMPLE_DIR}

# Generating decontaminated table
Rscript decontam_final.R

module purge
