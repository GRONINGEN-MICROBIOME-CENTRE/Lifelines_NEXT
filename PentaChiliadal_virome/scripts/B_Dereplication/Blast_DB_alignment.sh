#!/bin/bash
#SBATCH --job-name=Blast_DB_alignment
#SBATCH --output=Blast_DB_alignment.out
#SBATCH --mem=150gb
#SBATCH --time=6-0
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

contig_file=$1 #path to FASTA file with predicted viral contigs

mkdir -p BLASTN_RESULTS

# Load BLAST
module purge; ml BLAST+/2.13.0-gompi-2022a; ml list

#First, create a blast+ database:
makeblastdb \
	-in $contig_file \
	-dbtype nucl \
	-out BLASTN_RESULTS/viruses

#Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn \
	-query $contig_file \
	-db BLASTN_RESULTS/viruses \
	-outfmt '6 std qlen slen' \
	-max_target_seqs 10000 \
	-out BLASTN_RESULTS/viruses_blast.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}
