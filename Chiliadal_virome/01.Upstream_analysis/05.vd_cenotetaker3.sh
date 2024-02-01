#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --error=./err/05.ct3/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/05.ct3/VD_Chiliadal_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

##### Attempted to use, ended up not including Cenote-Taker 3 to the virus discovery pipeline for the time being as I did not find an option to skip generation of visualization files that filled up my scratch limit after running 1 sample. 


SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- SWITCH TO TMP --- 
mkdir -p ${TMPDIR}/${SAMPLE_ID}/CT3

cd ${TMPDIR}/${SAMPLE_ID}/CT3 # since Mike has designed the logger checker this way & I do not want to rewrite his scripts

echo "$(pwd)"
# --- LOAD MODULES --- 
module purge
module load Anaconda3
conda activate /scratch/p282752/tools/conda_envs/ct3_env

# --- RUNNING Cenote-Taker3 ---
echo "> Running Cenote-Taker3"

cenotetaker3 \
	-c /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	-r CenoteTaker3 \
	-p F \
	-t ${SLURM_CPUS_PER_TASK}

# -p is FALSE since we have VLP-enriched data and will be prunning prophages at the later stage	

# --- MOVING TO /SCRATCH ---
echo "> Moving results to /scratch"
rm -r ./CenoteTaker3/ct_processing # intermediate
rm -r ./CenoteTaker3/sequin_and_genome_maps # we will re-create it afterwards for selected contigs

mkdir -p /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/virome_discovery/CenoteTaker3
cp ./CenoteTaker3/CenoteTaker3_cenotetaker.log /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/virome_discovery/CenoteTaker3/
cp ./CenoteTaker3/CenoteTaker3_virus_summary.tsv /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/virome_discovery/CenoteTaker3/
cp ./CenoteTaker3/final_genes_to_contigs_annotation_summary.tsv /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/virome_discovery/CenoteTaker3/

conda list
conda deactivate

module list

module purge
