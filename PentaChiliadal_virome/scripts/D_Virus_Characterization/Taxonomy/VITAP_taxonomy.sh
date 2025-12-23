#!/bin/bash
#SBATCH --job-name=VITAP_tax
#SBATCH --output=VITAP_tax.out
#SBATCH --mem=120gb
#SBATCH --time=3-0
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_genomes=$1 #path to FASTA file with viral sequences
output_dir=$2 #path to output directory

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/VITAP_env; conda list 

# Get VITAP taxonomy
VITAP assignment \
    -i $viral_genomes \
    -d /scratch/hb-llnext/databases/VITAP_DB_MSL37_RefSeq209_IMGVR/ \
    -o $output_dir \
    -p ${SLURM_CPUS_PER_TASK}


