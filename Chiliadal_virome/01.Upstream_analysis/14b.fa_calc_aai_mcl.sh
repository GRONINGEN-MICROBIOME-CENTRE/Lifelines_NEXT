#!/bin/bash
#SBATCH --job-name=aai_calculation
#SBATCH --output=./out/14.fan/aai_debug.out
#SBATCH --mem=64gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

INPUT=$1
WORK_DIR=$2

echo "Viral proteomes file: $(realpath "${INPUT}")"
echo "Work directory: $(realpath "${WORK_DIR}")"

INPUT="${INPUT%.faa}"

echo ${INPUT}

# --- LOAD MODULES --- 
module purge
Python/3.12.3-GCCcore-13.3.0
module load SciPy-bundle/2024.05-gfbf-2024a
module load Biopython/1.84-foss-2024a
module list

# --- CALCULATE AAI ---
# using python3-updated https://github.com/snayfach/MGV/blob/master/aai_cluster/amino_acid_identity.py

python amino_acid_identity_p3.py \
	--in_faa ${INPUT}.faa \
	--in_blast ${WORK_DIR}/blastp.tsv \
	--out_tsv ${WORK_DIR}/aai.tsv

# --- FILTER AAI ---
# using https://github.com/snayfach/MGV/blob/master/aai_cluster/filter_aai.py

# GENUS LEVEL:
python filter_aai.py \
	--in_aai ${WORK_DIR}/aai.tsv \
	--min_percent_shared 20 \
	--min_num_shared 16 \
	--min_aai 40 \
	--out_tsv ${WORK_DIR}/genus_edges.tsv

# FAMILY LEVEL:
python filter_aai.py \
	--in_aai ${WORK_DIR}/aai.tsv \
	--min_percent_shared 10 \
	--min_num_shared 8 \
	--min_aai 20 \
	--out_tsv ${WORK_DIR}/family_edges.tsv

# --- Perform MCL-based clustering ---

# GENUS CLUSTERS:
/home1/p282752/local/bin/mcl ${WORK_DIR}/genus_edges.tsv \
	-te 8 \
	-I 2.0 \
	--abc \
	-o ${WORK_DIR}/genus_clusters.txt

# FAMILY CLUSTERS:
/home1/p282752/local/bin/mcl ${WORK_DIR}/family_edges.tsv \
	-te 8 \
	-I 1.2 \
	--abc \
	-o ${WORK_DIR}/family_clusters.txt

# --- CLEAN ENV --- 
module purge

# --- LOAD MODULES ---
module load R/4.3.2-gfbf-2023a

# Reformatting clustering info:

# GENUS:
tr '\t' ',' < ${WORK_DIR}/genus_clusters.txt > ${WORK_DIR}/genus_clusters_commas.txt
awk -v prefix="Genus_" '{printf "%s%05d\t%s\n", prefix, NR, $0}' ${WORK_DIR}/genus_clusters_commas.txt > ${WORK_DIR}/genus_clusters_labeled.tsv

Rscript dereplication_stat.R ${WORK_DIR}/genus_clusters_labeled.tsv

rm ${WORK_DIR}/genus_clusters_commas.txt
rm ${WORK_DIR}/genus_clusters_labeled.tsv

# FAMILY:
tr '\t' ',' < ${WORK_DIR}/family_clusters.txt > ${WORK_DIR}/family_clusters_commas.txt
awk -v prefix="Family_" '{printf "%s%04d\t%s\n", prefix, NR, $0}' ${WORK_DIR}/family_clusters_commas.txt > ${WORK_DIR}/family_clusters_labeled.tsv

Rscript dereplication_stat.R ${WORK_DIR}/family_clusters_labeled.tsv

rm ${WORK_DIR}/family_clusters_commas.txt
rm ${WORK_DIR}/family_clusters_labeled.tsv

# --- CLEAN ENV ---
module purge
