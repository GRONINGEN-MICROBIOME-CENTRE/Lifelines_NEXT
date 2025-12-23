#!/bin/bash
#SBATCH --job-name=Cluster_proteins_MMseqs2
#SBATCH --output=Cluster_proteins_MMseqs2.out
#SBATCH --mem=60gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_proteins_file=$1 #FASTA file with predicted viral proteins
viral_proteins="$(basename "${viral_proteins_file}")" #extract filename
viral_proteins="${viral_proteins%.*}" #extract filename without the extension

# Create empty directories
mkdir -p viral_proteins_db cluster_db msa_db profile_db consensus_db 
mkdir -p profile_consensus_search_db profile_cluster_db merged_cluster_db unaligned_db aligned_db ${viral_proteins}_aligned_unpacked_db_10seq

# Clean environment, load modules
mamba activate /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/MMseqs2; mamba list

# Before clustering, convert the FASTA DB into the MMseqs2 database (DB) format
mmseqs createdb \
    $viral_proteins_file \
    viral_proteins_db/${viral_proteins}_db

# Perform the clustering of the proteins
# By default a 3 step clustering is performed (increase in speed and allows large cluster sizes)
# --cluster-reassign can reassign some sequences to their correct clusters as in each iteration of the cascade clustering the representative seq might change
# -s sensitivity - 1 fast to 7.5 more sensitive
# -c coverage
# --cov-mode 0 coverage of both query and target
# --cluster-steps 3 #cascaded clustering steps - default 3 
# --cluster-mode 0 #set-cover - greedy
# --cluster-steps 3 #cascaded clustering steps - default 3 
# --cluster-reassign 1 #set to 1 -yes- to correct possible errors in cascaded clustering
mmseqs cluster \
    viral_proteins_db/${viral_proteins}_db \
    cluster_db/${viral_proteins}_cluster_db \
    tmp \
    -s 7.0 \
    -c 0.8 \
    --cov-mode 0 \
    -e 1e-3 \
    --cluster-mode 0 \
    --cluster-steps 3 \
    --cluster-reassign 1 \
    --kmer-per-seq 50

rm -rf tmp

# Generate a MSA from each cluster
# result2msa produces an MSA using a centre star alignment without insertions in the query
# --msa-format-mode 2 #2: aligned FASTA DB
mmseqs result2msa \
    viral_proteins_db/${viral_proteins}_db \
    viral_proteins_db/${viral_proteins}_db \
    cluster_db/${viral_proteins}_cluster_db \
    msa_db/${viral_proteins}_msa_db \
    --msa-format-mode 2

# Compute profiles from the MSA
# A) Match-mode 1 + match-ratio 0.5 turns all columns with at least 50% residues to profile columns
# B) Match-mode 0: All columns of the first sequence except gaps ‘-’ will be turned into profile columns
# *B option can be used for centre star MSAs but is risky for large MSAs (1st seq might not be representative of the MSA)
# --msa-type 2 #2: FASTA
# --match-mode 1  #1: columns that have a residue in --match-ratio of all sequences are kept
mmseqs msa2profile \
    --msa-type 2 \
    --match-mode 1 \
    --match-ratio 0.5 \
    msa_db/${viral_proteins}_msa_db \
    profile_db/${viral_proteins}_profile_db

# Extract the consensus sequence of the profile into a normal MMseqs2 DB
mmseqs profile2consensus \
    profile_db/${viral_proteins}_profile_db \
    consensus_db/${viral_proteins}_consensus_db

# Do a consensus vs. profile search (which is much more sensitive than the protein vs. protein search that was performed in the initial clustering)  
# -a 1 #include additional information about the alignment      
mmseqs search \
    profile_db/${viral_proteins}_profile_db \
    consensus_db/${viral_proteins}_consensus_db \
    profile_consensus_search_db/${viral_proteins}_profile_consensus_search_db \
    tmp \
    --cov-mode 0 \
    -c 0.9 \
    -s 8.0 \
    -e 1e-4 \
    --add-self-matches 1 \
    -a 1 \
    
rm -rf tmp

# Cluster the profiles (and, consequently, the clusters they represent) to get larger and more diverse clusters
mmseqs clust \
    profile_db/${viral_proteins}_profile_db \
    profile_consensus_search_db/${viral_proteins}_profile_consensus_search_db \
    profile_cluster_db/${viral_proteins}_profile_cluster_db
    
# Merge two clusterings (cluster_db and profile_cluster_db) into one result database
mmseqs mergeclusters \
    viral_proteins_db/${viral_proteins}_db \
    merged_cluster_db/${viral_proteins}_merged_cluster_db \
    cluster_db/${viral_proteins}_cluster_db \
    profile_cluster_db/${viral_proteins}_profile_cluster_db
       
# Generate a TSV formatted output file with the final clusters
mmseqs createtsv \
    viral_proteins_db/${viral_proteins}_db \
    viral_proteins_db/${viral_proteins}_db \
    merged_cluster_db/${viral_proteins}_merged_cluster_db \
    ${viral_proteins}_clusters.tsv

# Build MSAs using star center alignment of the final clusters (2nd clustering) 
mmseqs result2msa \
    viral_proteins_db/${viral_proteins}_db \
    viral_proteins_db/${viral_proteins}_db \
    merged_cluster_db/${viral_proteins}_merged_cluster_db \
    aligned_db/${viral_proteins}_aligned_db \
    --msa-format-mode 2

# Finally, generate a FASTA file from the MMSeqs2 formatted DB
mmseqs unpackdb \
    --unpack-suffix ".faa" \
    --unpack-name-mode 0 \
    aligned_db/${viral_proteins}_aligned_db \
    ${viral_proteins}_aligned_unpacked_db

