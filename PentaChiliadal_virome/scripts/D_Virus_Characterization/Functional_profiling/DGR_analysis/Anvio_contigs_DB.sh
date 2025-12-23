#!/bin/bash
#SBATCH --job-name=Anvio_contigs_DB 
#SBATCH --output=Anvio_contigs_DB.out
#SBATCH --mem=32gb
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

vOTUs=$1 #path to FASTA file with vOTU genomes
external_gene_calls=$2 #path to table with external gene calls

echo -e '\n---- RUNNING anvio ----'

# Clean environment, load modules 
ml purge; ml Anaconda3; ml list
conda activate /scratch/hb-llnext/conda_envs/anvio-8; conda list

# Run anvi'o
anvi-gen-contigs-database \
    -f $vOTUs \
    -o contigs.db \
    --external-gene-calls $external_gene_calls \
    -n "vOTUs" \
    -T ${SLURM_CPUS_PER_TASK} 
