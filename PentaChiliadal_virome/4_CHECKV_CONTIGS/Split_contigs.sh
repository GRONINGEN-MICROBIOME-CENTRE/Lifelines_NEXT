#!/bin/bash
#SBATCH --job-name=split_contig_file
#SBATCH --output=split_contig.out
#SBATCH --mem=8gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with the predicted viral contigs
mkdir -p SPLIT_CONTIGS

# Load BBMap and split the contig file in 50 files
module purge; ml BBMap; module list
partition.sh \
	in=$contig_file \
	out=SPLIT_CONTIGS/Viral_sequences_%.fa \
	ways=200

# Generate list of files
cd SPLIT_CONTIGS
ls *fa  > list.txt
