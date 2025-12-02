#!/bin/bash
#SBATCH --job-name=AT_generation
#SBATCH --output=./out/12.mrmc/at_generation_w_mgs_nc.out
#SBATCH --mem=256gb
#SBATCH --time=23:59:30
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

COV_DIR=$1
ETOF=$2
MAP_DIR=$(dirname "${COV_DIR}")

echo "BEDtools coverage source: $(realpath "${COV_DIR}")"
echo "ETOF source: $(realpath "${ETOF}")"
echo "Mapping directory: $(realpath "${MAP_DIR}")"

module load R/4.3.2-gfbf-2023a

# Combine coverage
Rscript combine_coverage.R ${COV_DIR} *.coverage.txt

# Create RPKM table for decontamination
Rscript RPKM_table2.R ${MAP_DIR} ${ETOF}

mkdir -p ../VIR_DB/VLP_MGS_decontamination/keep_vOTUs_per_sample
# Create decontamination guides
Rscript Filtering_vOTUs_and_samples_to_decontam2.R

module purge
