#!/bin/bash
#SBATCH --job-name=Bin_dereplication
#SBATCH --output=Bin_dereplication.out
#SBATCH --mem=60gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

bins_dir=$1 #directory with the bins output of metaWRAP
output_dir=$2 #output directory

find $bins_dir -type f -name *_bin*fa > bins_paths.txt

echo -e '\n-------------------- DEREPLICATING BINS --------------------'

# Clean environment, load modules 
module purge; ml dRep; module list

dRep dereplicate \
                -p ${SLURM_CPUS_PER_TASK} \
                $output_dir \
                -g bins_paths.txt \
                -comp 75 \
                -con 10 \
                --S_algorithm ANImf \
                --S_ani 0.98 \
                --cov_thresh 0.25

rm bins_paths.txt
