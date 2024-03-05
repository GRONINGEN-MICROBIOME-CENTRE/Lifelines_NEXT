#!/bin/bash
#SBATCH --job-name=PostDiscovery
#SBATCH --error=./err/08.chv/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/08.chv/VD_Chiliadal_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

mkdir -p ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/CheckV_pruning

# PREPARING TMP
mkdir -p ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning

# --- LOAD MODULES --- 
module purge
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# --- RUNNING CHECKV --- 
checkv end_to_end \
    ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
    ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning \
    -t ${SLURM_CPUS_PER_TASK} \
    -d /scratch/hb-llnext/databases/checkv-db-v1.5

# --- CONCATENATING VIRUSES WITH PRUNED PROVIRUSES ---
sed -i 's/\ /_/g' ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning/proviruses.fna 
sed -i 's/\//_/g' ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning/proviruses.fna 

cat ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning/*viruses.fna > ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral.fasta

# --- SAVING BACTERIAL CONTAMINATION REPORT ---
mv ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning/contamination.tsv ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/CheckV_pruning/


# --- REMOVING OLD CHECKV RUN RESULTS ---
rm -r ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning
mkdir ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning

# --- RERUNNING CHECKV --- 
checkv end_to_end \
    ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral.fasta \
    ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning \
    -t ${SLURM_CPUS_PER_TASK} \
    -d /scratch/hb-llnext/databases/checkv-db-v1.5

mv ${TMPDIR}/${SAMPLE_ID}/CheckV_pruning/quality_summary.tsv ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/CheckV_pruning/

# --- CREATE NEW VIRUS CONTIGS METADATA AND IDs ---
module load R
Rscript New_contigs_ID_and_metadata.R /scratch/p282752/ANALYSIS_CHILIADAL/SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/ ${SAMPLE_ID}

# --- RENAME CONTIGS ---
cp ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral.fasta ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral_renamed.fasta
awk 'NR>1' ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/Extended_TOF  | awk '{print $15"\t"$17}' | while IFS=$'\t' read -r old_id new_id; do
        sed "s/>$old_id\b/>$new_id/" -i ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral_renamed.fasta
done

module purge
