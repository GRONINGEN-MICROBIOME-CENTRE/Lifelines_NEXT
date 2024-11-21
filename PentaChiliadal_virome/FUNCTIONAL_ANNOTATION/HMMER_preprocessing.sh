#!/bin/bash
#SBATCH --job-name=HMMER_preprocessing 
#SBATCH --output=HMMER_preprocessing.out
#SBATCH --mem=4gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=himem

protein_file=$1 #path to FASTA file with the predicted viral proteins
HMM_DB=$2 #directory with HMM profiles of public DBs

echo -e '\n-------------------- SPLITTING PROTEIN FASTA file ------------------'

mkdir -p SPLIT_PROTEINS

# Load SeqKit and split the contig file in 1000 files
module purge; ml SeqKit; module list
seqkit split2 $protein_file -p 1000 -O SPLIT_PROTEINS -f

# Generate list of files
cd SPLIT_PROTEINS
ls *faa  > list.txt 

echo -e '\n-------------------- PROCESSING HMMs from DBs ------------------'

# Clean environment and load modules 
module purge; ml HMMER; module list

# Convert the individual hmms to the most recent version using hmmconvert
for i in $HMM_DB/*hmm; do
    DB_name=$(basename "$i")
    base="${DB_name%.*}"
    hmmconvert $i > $HMM_DB/${base}_new.hmm
done

# Create an HMM database flatfile combining NCBIFam and Pfam
cat $HMM_DB/PGAP_NCBIFam_new.hmm $HMM_DB/Pfam-A_new.hmm > $HMM_DB/Pfam-A_NCBIFam.hmm
