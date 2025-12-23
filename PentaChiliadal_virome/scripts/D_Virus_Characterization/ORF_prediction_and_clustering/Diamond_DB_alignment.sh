#!/bin/bash
#SBATCH --job-name=Diamond_DB_alignment
#SBATCH --output=Diamond_DB_alignment.out
#SBATCH --mem=40gb
#SBATCH --time=02:20:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

protein_file=$1 #path to FASTA file with predicted viral proteins 

# Load Diamond
module purge; ml DIAMOND; ml list

#First, create a DIAMOND database:
diamond makedb \
	--in $protein_file \
  --db viral_proteins \
	--threads ${SLURM_CPUS_PER_TASK}

#Next, use diamond blastp to perform all-vs-all comparison of sequences:
diamond blastp \
	--query $protein_file \
	--db viral_proteins \
	--out blastp.tsv \
	--outfmt 6 \
	--evalue 1e-5 \
	--max-target-seqs 10000 \
	--query-cover 50 \
	--subject-cover 50
