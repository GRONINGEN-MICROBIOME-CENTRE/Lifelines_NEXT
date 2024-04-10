#!/bin/bash
#SBATCH --job-name=PostDiscovery
#SBATCH --error=/scratch/p282752/ANALYSIS_CHILIADAL/scripts/err/08.sti/VD_Chiliadal_%A_%a.err
#SBATCH --output=/scratch/p282752/ANALYSIS_CHILIADAL/scripts/out/08.sti/VD_Chiliadal_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=25:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

##############################################
# Use-cases: in-house negative controls only #
##############################################

FRAG_LIST=$1

echo ${FRAG_LIST}

FRAG_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${FRAG_LIST} | cut -d "_" -f1)

echo "FRAG_ID=${FRAG_ID}"

SAMPLE_ID=$2

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- RENAME CONTIGS ---
awk 'NR>1' /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/yasli/${FRAG_ID}_Extended_TOF | awk '{print $15"\t"$17}' | while IFS=$'\t' read -r old_id new_id; do
        sed "s/>$old_id\b/>$new_id/" -i /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/yasli/${FRAG_ID}.fa
done

rm /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/01_sc_assembly/yasli/${FRAG_ID}_Extended_TOF
module purge
