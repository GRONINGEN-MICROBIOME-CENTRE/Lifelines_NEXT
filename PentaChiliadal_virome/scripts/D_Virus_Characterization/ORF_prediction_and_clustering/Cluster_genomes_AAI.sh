#!/bin/bash
#SBATCH --job-name=Cluster_genomes_AAI
#SBATCH --output=Cluster_genomes_AAI.out
#SBATCH --mem=150gb
#SBATCH --time=02:29:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

protein_file=$1 #path to FASTA file with predicted viral proteins 
diamond_db=$2 #path to diamond database generated from all viral sequences

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 CheckV; module list

# Run clustering of viruses with MGV scripts

# Compute AAI from BLAST results (Amino acid identity is computed based on the average BLAST percent identity between all genes shared between each pair of genomes (E-value <1e-5))
python amino_acid_identity.py --in_faa $protein_file --in_blast $diamond_db --out_tsv viral_aai.tsv

# Filter edges and prepare MCL input (keeping edges between genomes with >=20% AAI and genomes with either 8 shared genes or at least 20% of shared genes (relative to both genomes)
python filter_aai.py --in_aai viral_aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv genus_edges.tsv
python filter_aai.py --in_aai viral_aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv family_edges.tsv

# Load MCL conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/MCL_env; conda list

#Finally, perform MCL-based clustering (In the output each row indicates the members belonging to each cluster (including singletons))
mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.txt
mcl family_edges.tsv -te 8 -I 1.2 --abc -o family_clusters.txt
