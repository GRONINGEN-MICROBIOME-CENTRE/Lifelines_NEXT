#!/bin/bash

#SBATCH --job-name=DSP_DBA
#SBATCH --output=./out/DSP/DSP_DBA_%A.out
#SBATCH --error=./err/DSP/DSP_DBA_%A.err
#SBATCH --time=03:59:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

VIRAL_PROTEINS=$1 #path to FASTA file with viral sequences
OUTPUT_DIR=$2 #path to output directory

echo "Viral protein file: $(realpath "$VIRAL_PROTEINS")"
echo "Output directory: $(realpath "$OUTPUT_DIR")"

cd /home1/p309176/dbAPIS

module purge
module load HMMER/3.4-gompi-2023a

module list

hmmscan \
	--cpu 8 \
	--domtblout "${OUTPUT_DIR}/hmmscan.out" \
	--noali dbAPIS.hmm \
	"$VIRAL_PROTEINS"


module purge
module load DIAMOND/2.1.8-GCC-12.2.0

module list

diamond blastp \
	--db APIS_db \
	-q "$VIRAL_PROTEINS" \
	-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
	-o "${OUTPUT_DIR}/diamond.out" \
	--max-target-seqs 10000

module purge

bash /home1/p309176/dbAPIS/parse_annotation_result.sh "${OUTPUT_DIR}/hmmscan.out" "${OUTPUT_DIR}/diamond.out"


#module load HMMER/3.4-gompi-2023a

# run hmmscan for your amino acid sequences
#hmmscan --domtblout hmmscan.out --noali dbAPIS.hmm your_sequence.faa
# "--domtblout" option produces the space-separated domain hits table. There is one line for each domain. "--noali" option is used to omit the alignment section from output and reduce the output volume. More hmmscan information please see http://eddylab.org/software/hmmer/Userguide.pdf

#module purge

#module load DIAMOND/2.1.8-GCC-12.2.0


# run diamond for your amino acid sequences
#diamond blastp --db APIS_db -q your_sequence.faa -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen -o diamond.out --max-target-seqs 10000
# "-f 6" option generates tabular-separated format (a BLAST output format using the option -outfmt 6), which composed of the customized fields. "--max-target-seqs" means maximum number of target sequences to report alignments for. More diamond details please see https://github.com/bbuchfink/diamond/wiki/1.-Tutorial


# run script to parse annotation output files
#bash parse_annotation_result.sh hmmscan.out diamond.out
# This will generate parsed output files of hmmscan and diamond respectively

