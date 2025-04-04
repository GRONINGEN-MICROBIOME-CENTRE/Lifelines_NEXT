#!/bin/bash
#SBATCH --job-name=PostDiscovery
#SBATCH --output=./out/11.lfs/DB_launch_%A_%a.out
#SBATCH --mem=8gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

contig_file=$1 #path to FASTA file with the predicted viral contigs

DB_DIR=$(dirname ${contig_file})

SPLIT_FACTOR=$2

mkdir -p ${DB_DIR}/fragmented
cd ${DB_DIR}/fragmented

echo "Splitting the fasta in ${SPLIT_FACTOR} pieces"
/scratch/p282752/ANALYSIS_CHILIADAL/scripts/fastasplitn ${contig_file} ${SPLIT_FACTOR}
ls *.fa | sed 's/\.fa//g'> split_list

cd /scratch/p282752/ANALYSIS_CHILIADAL/scripts
sbatch --array=1-${SPLIT_FACTOR} 11.pd_lifestyle.sh ${DB_DIR}/fragmented/split_list

module purge
