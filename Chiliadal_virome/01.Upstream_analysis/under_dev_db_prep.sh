#!/bin/bash
#SBATCH --job-name=VirusDiscovery
#SBATCH --output=./out/09.dbs/DB_prep_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=2

DB_LIST=$1

echo "DB list: ${DB_LIST}"

DB=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DB_LIST})

echo "DB=${DB}"

DB_DIR=$(dirname ${DB})

DB_ID=$(basename ${DB_DIR})
echo "DB_ID=${DB_ID}"
echo "DB_DIR=${DB_DIR}"

mkdir -p ${DB_DIR}/Prophage_pruning

# --- GETTING contig lengths ---

# --- LOAD MODULES ---
module load bioawk
module list

bioawk -c fastx '{ print $name, length($seq) }' \
	${DB} \
	> ${DB_DIR}/Prophage_pruning/POST_CBR_LENGTH

module purge

# PREPARING TMP
mkdir -p ${TMPDIR}/${DB_ID}/Prophage_pruning

# --- LOAD MODULES ---
module load ARAGORN/1.2.41-foss-2021b
module load Python/3.9.5-GCCcore-10.3.0
source /scratch/p282752/tools/python_envs/geNomad/bin/activate
module list

# --- RUNNING geNomad for prepruning ---
echo "> Running geNomad"

genomad \
        end-to-end \
        --enable-score-calibration \
        --cleanup \
        ${DB} \
        ${TMPDIR}/${DB_ID}/Prophage_pruning \
	/scratch/p282752/databases/genomad_db

genomad --version

deactivate

# 1. Get putative virus contigs that have DTR or ITR:
grep -E "(DTR|ITR)" \
	${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_summary/${DB_ID}_virus_summary.tsv | \
	awk -F '\t' '{print $1}' \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs

# 2. Get the provirus sequences:
## We have to take it from the sample virus summary, because some proviruses might be discarded by the post-classification filtering process
grep "Provirus" \
	${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_summary/${DB_ID}_virus_summary.tsv | \
	awk -F '\t' '{print $1}' | \
	sort | \
	uniq \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus_vIDs

awk -F '\|' '{print $1}' \
	${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus_vIDs | \
	sort | \
	uniq \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus_source_vIDs

# 3. Get the rest of putative virus contigs (including those not recognized by geNomad as viral)
grep '>' ${DB} | \
	sed 's/>//g' | \
	sort \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_IDs

cat ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs \
       ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus_source_vIDs | \
       sort | \
       uniq \
       > ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR_provirus_vIDs

comm -23 \
	<(sort ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_IDs) \
	<(sort ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR_provirus_vIDs) \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_rest_vIDs

# --- LOAD MODULES ---
module purge
module load seqtk/1.3-GCC-11.3.0
module list

# 4. Create fastas
seqtk \
        subseq \
        -l60 \
	${DB} \
	${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR.fasta

seqtk \
        subseq \
        -l60 \
	${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_find_proviruses/${DB_ID}_provirus.fna \
	${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus_vIDs \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus.fasta

seqtk \
        subseq \
        -l60 \
        ${DB} \
	${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_rest_vIDs \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_rest.fasta

cat ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_provirus.fasta \
	${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_rest.fasta \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/POST_GND.fasta

N_VIRUSES_PRE=$(grep '>' ${DB} | wc -l)
N_VIRUSES_POST=$(grep '>' ${TMPDIR}/${DB_ID}/Prophage_pruning/POST_GND.fasta | sed 's/|.*//' | sort | uniq | wc -l)
N_VIRUSES_COMPLETE=$(cat ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs | wc -l)

if [ ${N_VIRUSES_PRE} -eq $((N_VIRUSES_POST + N_VIRUSES_COMPLETE)) ]; then
    echo "All putative viruses were processed"
else
    echo "N viruses in input and intermediate output is different"
fi

# --- RUNNING CheckV for pruning ---

# --- LOAD MODULES --- 
module purge
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# 1. Running CheckV
checkv end_to_end \
	${TMPDIR}/${DB_ID}/Prophage_pruning/POST_GND.fasta \
	${TMPDIR}/${DB_ID}/Prophage_pruning/CHV \
	-t ${SLURM_CPUS_PER_TASK} \
    	-d /scratch/hb-llnext/databases/checkv-db-v1.5

# --- CONCATENATING VIRUSES WITH PRUNED PROVIRUSES ---
sed -i 's/\ /_/g' ${TMPDIR}/${DB_ID}/Prophage_pruning/CHV/proviruses.fna
sed -i 's/\//_/g' ${TMPDIR}/${DB_ID}/Prophage_pruning/CHV/proviruses.fna

# --- RUNNING CheckV for quality assessment ---

# 1. Concatenating complete genomes & CheckV-pruned seqeunces & CheckV-untouched sequences
cat ${TMPDIR}/${DB_ID}/Prophage_pruning/geNomad_DTR_ITR.fasta \
	${TMPDIR}/${DB_ID}/Prophage_pruning/CHV/proviruses.fna \
	${TMPDIR}/${DB_ID}/Prophage_pruning/CHV/viruses.fna \
	> ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_pruned.fasta

# 2. Running CheckV
checkv end_to_end \
        ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_pruned.fasta \
        ${TMPDIR}/${DB_ID}/Prophage_pruning/CHV_new \
        -t ${SLURM_CPUS_PER_TASK} \
        -d /scratch/hb-llnext/databases/checkv-db-v1.5

# --- COPYING the output ---
#cp -r ${TMPDIR}/${DB_ID}/Prophage_pruning/* ${DB_DIR}/Prophage_pruning/
cp ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_pruned.fasta ${DB_DIR}/Prophage_pruning
cp ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_summary/${DB_ID}_plasmid_summary.tsv ${DB_DIR}/Prophage_pruning
cp ${TMPDIR}/${DB_ID}/Prophage_pruning/${DB_ID}_summary/${DB_ID}_virus_summary.tsv ${DB_DIR}/Prophage_pruning
cp ${TMPDIR}/${DB_ID}/Prophage_pruning/CHV/contamination.tsv ${DB_DIR}/Prophage_pruning
cp ${TMPDIR}/${DB_ID}/Prophage_pruning/CHV_new/quality_summary.tsv ${DB_DIR}/Prophage_pruning

# --- CREATE NEW VIRUS CONTIGS METADATA AND IDs ---
module load R
Rscript New_DBcontigs_ID_and_metadata.R ${DB_DIR} ${DB_ID}

# --- RENAME CONTIGS ---
cp ${DB_DIR}/Prophage_pruning/${DB_ID}_pruned.fasta ${DB_DIR}/Prophage_pruning/${DB_ID}_pruned_renamed.fasta
awk 'NR>1' ${DB_DIR}/Prophage_pruning/Extended_TOF | awk -F '\t' '{print $35"\t"$1}' | while IFS=$'\t' read -r old_id new_id; do
	sed "s/>$old_id\b/>$new_id/" -i ${DB_DIR}/Prophage_pruning/${DB_ID}_pruned_renamed.fasta
done

# --- CLEANING output ---
rm ${DB_DIR}/Prophage_pruning/${DB_ID}_plasmid_summary.tsv
rm ${DB_DIR}/Prophage_pruning/${DB_ID}_virus_summary.tsv
rm ${DB_DIR}/Prophage_pruning/contamination.tsv
rm ${DB_DIR}/Prophage_pruning/POST_CBR_LENGTH
rm ${DB_DIR}/Prophage_pruning/quality_summary.tsv
mv ${DB_DIR}/Prophage_pruning/*.fasta ${DB_DIR}
mv ${DB_DIR}/Prophage_pruning/Extended_TOF ${DB_DIR}
rm -r ${DB_DIR}/Prophage_pruning

module purge
