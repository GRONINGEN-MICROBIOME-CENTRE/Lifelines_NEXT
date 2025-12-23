#!/bin/bash
#SBATCH --job-name=geNomad_ME
#SBATCH --output=geNomad_ME.out
#SBATCH --mem=100gb
#SBATCH --time=4-0
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with contigs
output_dir=$2 #path to output directory

echo -e '\n---- RUNNING geNomad ----'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/geNomad_conda/; conda list 

# Run geNomad
genomad end-to-end \
        --enable-score-calibration \
        --cleanup \
        $contig_file \
        $output_dir \
        --threads ${SLURM_CPUS_PER_TASK} \
        /scratch/hb-llnext/databases/geNomad_db/

conda deactivate

