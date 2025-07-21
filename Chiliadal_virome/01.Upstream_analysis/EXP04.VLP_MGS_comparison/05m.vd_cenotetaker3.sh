#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --output=./out/05.mct3/VD_PentaChiliadal_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- SWITCH TO TMP --- 
mkdir -p ${TMPDIR}/${SAMPLE_ID}/CT3

cd ${TMPDIR}/${SAMPLE_ID}/CT3 # since Mike has designed the logger checker this way & I do not want to rewrite his scripts

echo "$(pwd)"
# --- LOAD MODULES --- 
module purge
module load Anaconda3
conda activate /scratch/p282752/tools/conda_envs/ct3.2.1_env

# --- RUNNING Cenote-Taker3 ---
echo "> Running Cenote-Taker3"

cenotetaker3 \
	-c ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	-r CenoteTaker3 \
	-p F \
	-t ${SLURM_CPUS_PER_TASK}

# -p is FALSE since we did it same way for VLP-enriched data and will be prunning prophages at the later stage	

# --- MOVING TO /SCRATCH ---
echo "> Moving results to /scratch"
rm -r ./CenoteTaker3/ct_processing # intermediate
rm -r ./CenoteTaker3/sequin_and_genome_maps # we will re-create it afterwards for selected contigs

mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/CenoteTaker3
#cp ./CenoteTaker3/CenoteTaker3_cenotetaker.log ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/CenoteTaker3/ # decided to skip saving the log as .out is identical
cp ./CenoteTaker3/CenoteTaker3_virus_summary.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/CenoteTaker3/${SAMPLE_ID}_CenoteTaker3_virus_summary.tsv
cp ./CenoteTaker3/final_genes_to_contigs_annotation_summary.tsv ${SAMPLE_DIR}/05_CT3_all_genes_summary/${SAMPLE_ID}_final_genes_to_contigs_annotation_summary.tsv

conda list
conda deactivate

module list

module purge
