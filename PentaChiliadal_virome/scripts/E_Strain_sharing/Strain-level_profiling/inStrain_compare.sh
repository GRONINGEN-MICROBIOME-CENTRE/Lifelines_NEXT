#!/bin/bash
#SBATCH --job-name=inStrain_compare
#SBATCH --output=inStrain_compare_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

IS_profiles=$1 #path to directory with the inStrain profile results
vOTU_list=$2 #path to file with the list of vOTU names to be compared
metadata=$3  # Metadata with at least: sampleID[TAB]mother_or_infant
results_dir=$4 # Directory to store the result files
vOTU_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${vOTU_list})

echo '-------------------- RUNNING inStrain COMPARE ON vOTU '${vOTU_ID}' --------------------'

# Activate environment
ml purge; module load Python/3.10.8-GCCcore-12.2.0; ml list
source /scratch/hb-llnext/python_venvs/inStrain/bin/activate

# Temp directory for this vOTU
TMP_VOTU_DIR="${TMPDIR}/${vOTU_ID}"
mkdir -p "${TMP_VOTU_DIR}"

# First, select the inStrain profiles (samples) to be included in the comparison (breadth ≥ 0.75)
grep -rH -w ${vOTU_ID} $IS_profiles --include="*inStrain_scaffold_info.tsv" \
  | awk -F'\t' '$4 > 0.74 {print $1}' \
  | cut -d '/' -f1-8 \
  | sort | uniq > "${TMP_VOTU_DIR}/list_profiles_raw.txt"
  
profile_count=$(wc -l < "${TMP_VOTU_DIR}/list_profiles_raw.txt")

# Check if the number of profiles with the vOTU present is <2
if (( profile_count < 2 )); then
  echo "${vOTU_ID}" >> "Skipped_vOTUs_few_profiles.txt"
  echo "Only ${profile_count} profile(s) found for ${vOTU_ID}, skipping comparisons."
  exit 0
fi

# Get sample names
sed 's|.*/||;s/_inStrain$//' "${TMP_VOTU_DIR}/list_profiles_raw.txt" > "${TMP_VOTU_DIR}/samples_pass.txt"

# Filter metadata for these samples
grep -wf "${TMP_VOTU_DIR}/samples_pass.txt" "${metadata}" > "${TMP_VOTU_DIR}/meta_filtered.tsv"

# Separate mother/infant
awk '$2 == "mother" {print $1}' "${TMP_VOTU_DIR}/meta_filtered.tsv" > "${TMP_VOTU_DIR}/mothers.txt"
awk '$2 == "infant" {print $1}' "${TMP_VOTU_DIR}/meta_filtered.tsv" > "${TMP_VOTU_DIR}/infants.txt"

# Check if the number of maternal or infant profiles with the vOTU present is = 0
mother_count=$(wc -l < "${TMP_VOTU_DIR}/mothers.txt")
infant_count=$(wc -l < "${TMP_VOTU_DIR}/infants.txt")

if (( mother_count == 0 )); then
  echo "${vOTU_ID}" >> "Skipped_vOTUs_no_mothers.txt"
  echo "No maternal samples found for ${vOTU_ID}, skipping comparisons."
  exit 0
fi

if (( infant_count == 0 )); then
  echo "${vOTU_ID}" >> "Skipped_vOTUs_no_infants.txt"
  echo "No infant samples found for ${vOTU_ID}, skipping comparisons."
  exit 0
fi

# Output folder
mkdir -p ${TMP_VOTU_DIR}/${vOTU_ID}_output
echo "${vOTU_ID}" > ${TMP_VOTU_DIR}/${vOTU_ID}.txt

# Pairwise mother-infant comparisons
while read -r infant; do
  while read -r mother; do
    if [[ "$infant" != "$mother" ]]; then
      echo "Comparing $infant vs $mother for $vOTU_ID"
      outdir="${TMP_VOTU_DIR}/${vOTU_ID}_output/${infant}_${mother}"
      mkdir -p "$outdir"

      inStrain compare \
        -i ${IS_profiles}/${infant}_inStrain ${IS_profiles}/${mother}_inStrain \
        -o "$outdir" \
        --scaffolds ${TMP_VOTU_DIR}/${vOTU_ID}.txt \
        --min_cov 5 \
        --min_freq 0.05 \
        --fdr 1e-06 \
        -p ${SLURM_CPUS_PER_TASK}
    fi
  done < "${TMP_VOTU_DIR}/mothers.txt"
done < "${TMP_VOTU_DIR}/infants.txt"

# Move results to /scratch (add header from the first file)
comp_files=$(find "${TMP_VOTU_DIR}/${vOTU_ID}_output/" -type f -name '*comparisonsTable.tsv')
if [[ -n "$comp_files" ]]; then
    mkdir -p "${results_dir}"
    out_file="${results_dir}/${vOTU_ID}_comparisons.tsv"
    first=1
    for f in $comp_files; do
        if [[ $first -eq 1 ]]; then
            cat "$f" > "$out_file" # keep header
            first=0
        else
            tail -n +2 "$f" >> "$out_file"
        fi
    done
    echo "Results written to $out_file"
else
    echo "No comparison tables found for ${vOTU_ID}"
fi

# Clean TMP dir
rm -rf ${TMP_VOTU_DIR}

echo "-------------------- inStrain compare DONE for ${vOTU_ID} --------------------"
