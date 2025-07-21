#!/bin/bash
#SBATCH --job-name=VITAP_tax
#SBATCH --output=./out/14.fan/VITAP_tax.out
#SBATCH --mem=160gb
#SBATCH --time=2-0
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_genomes=$1 #path to FASTA file with viral sequences
output_dir=$2 #path to output directory

echo "Viral genomes file: $(realpath "$viral_genomes")"
echo "Output directory: $(realpath "$output_dir")"

# Clean environment, load modules and activate conda environment
module purge 
module load Anaconda3
module load prodigal
module list

conda activate /scratch/p282752/tools/conda_envs/vitap_env
conda list 

# Get VITAP taxonomy
/scratch/p282752/tools/conda_envs/VITAP/scripts/VITAP assignment \
    -i $viral_genomes \
    -d /scratch/p282752/tools/conda_envs/VITAP/DB_hybrid_MSL37_RefSeq209_IMGVR/ \
    -o $output_dir \
    -p ${SLURM_CPUS_PER_TASK}

conda deactivate
module purge
