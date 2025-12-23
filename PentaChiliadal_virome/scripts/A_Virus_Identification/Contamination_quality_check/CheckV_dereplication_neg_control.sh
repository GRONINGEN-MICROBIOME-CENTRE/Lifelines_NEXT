#!/bin/bash
#SBATCH --job-name=CheckV_der_neg_control
#SBATCH --output=CheckV_der_neg_control.out
#SBATCH --mem=20gb
#SBATCH --time=16:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #path to FASTA file with all the viral sequences (own viral database + public DBs)
blast_db=$2 #path to blast+ database generated from all viral sequences

mkdir -p DEREPLICATION_NEG_CONTROL_RESULTS

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 CheckV seqtk; module list

# Run dereplication for contigs with CheckV scripts

#Using the blast+ database as input, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/hb-llnext/PentaChiliadal/5_DEREPLICATION_vOTUs/anicalc.py \
    -i $blast_db  \
    -o DEREPLICATION_NEG_CONTROL_RESULTS/viruses_ani.tsv

#Finally, perform UCLUST-like clustering using the following -parameters (99% ANI + 85% AF):
python /scratch/hb-llnext/PentaChiliadal/5_DEREPLICATION_vOTUs/aniclust.py \
    --fna $contig_file \
    --ani DEREPLICATION_NEG_CONTROL_RESULTS/viruses_ani.tsv \
    --out DEREPLICATION_NEG_CONTROL_RESULTS/viral_clusters.tsv \
    --min_ani 99 \
    --min_tcov 85 \
    --min_qcov 0
