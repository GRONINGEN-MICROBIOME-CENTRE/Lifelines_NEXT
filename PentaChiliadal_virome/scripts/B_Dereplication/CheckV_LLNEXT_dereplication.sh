#!/bin/bash
#SBATCH --job-name=CheckV_LLNEXT_dereplication
#SBATCH --output=CheckV_LLNEXT_dereplication.out
#SBATCH --mem=20gb
#SBATCH --time=16:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

contig_file=$1 #file with all the renamed viral genomes after removal of negative control sequences
blast_db=$2 #path to blast+ database generated from all viral sequences

mkdir -p LLNEXT_DEREPLICATION

# Load SeqKit
module purge; ml SeqKit Python/3.10.8-GCCcore-12.2.0 CheckV seqtk; ml list

# Extract only LLNEXT genomes from the FASTA file with the renamed viral sequences
seqkit grep -r -p "LLNEXT_" $contig_file > LLNEXT_viral_genomes.fa

# Exclude all comparisons with negative controls and contigs clustering with them from the blast results
# Alsp, extract only comparisons between LLNEXT genomes 
grep -v -f Excluded_sequences.txt $blast_db | grep -Ev "MGV|GPD|Gulyaeva|Benler|Shah|IMG_VR|ELGV|RefSeq|Guerin|Yutin" > BLASTN_RESULTS/LLNEXT_viruses_blast_no_neg.tsv


echo -e '\n-------------------- DEDUPLICATING LLNEXT VIRAL GENOMES --------------------'

# Run deduplication for contigs with CheckV scripts

# Using the blast+ database as input, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/hb-llnext/PentaChiliadal/5_DEREPLICATION_vOTUs/anicalc.py \
    -i BLASTN_RESULTS/LLNEXT_viruses_blast_no_neg.tsv  \
    -o LLNEXT_DEREPLICATION/LLNEXT_viruses_ani.tsv

# Perform UCLUST-like clustering (99% ANI + 95% AF):
python /scratch/hb-llnext/PentaChiliadal/5_DEREPLICATION_vOTUs/aniclust.py \
    --fna LLNEXT_viral_genomes.fa \
    --ani LLNEXT_DEREPLICATION/LLNEXT_viruses_ani.tsv \
    --out LLNEXT_DEREPLICATION/LLNEXT_viral_genome_clusters.tsv \
    --min_ani 99 \
    --min_tcov 95 \
    --min_qcov 0

# Extract the representatives viral genomes (deduplicated genomes)
cat LLNEXT_DEREPLICATION/LLNEXT_viral_genome_clusters.tsv  | cut -f1 > LLNEXT_dedup_viral_genomes.txt
seqtk subseq LLNEXT_viral_genomes.fa LLNEXT_dedup_viral_genomes.txt > LLNEXT_dedup_viral_genomes.fa

echo -e '\n-------------------- DEREPLICATING LLNEXT VIRAL GENOMES INTO vOTUs --------------------'

# Perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
# Use deduplicated genomes as input
python /scratch/hb-llnext/PentaChiliadal/5_DEREPLICATION_vOTUs/aniclust.py \
    --fna LLNEXT_dedup_viral_genomes.fa \
    --ani LLNEXT_DEREPLICATION/LLNEXT_viruses_ani.tsv \
    --out LLNEXT_DEREPLICATION/LLNEXT_vOTU_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0

# Generate a TXT and a FASTA file with the representative sequences of each LLNEXT vOTU (only LLNEXT)
cat LLNEXT_DEREPLICATION/LLNEXT_vOTU_clusters.tsv | cut -f1 > rep_seqs_LLNEXT_vOTUs.txt
seqtk subseq LLNEXT_viral_genomes.fa rep_seqs_LLNEXT_vOTUs.txt > rep_seqs_LLNEXT_vOTUs.fa

# Set permissions 
chmod 440 rep_seqs_LLNEXT_vOTUs.txt rep_seqs_LLNEXT_vOTUs.fa
