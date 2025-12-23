#!/bin/bash
#SBATCH --job-name=Anvio_gen_var
#SBATCH --output=Anvio_gen_var_%A_%a.out
#SBATCH --mem=90gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

vOTU_DB=$1 # path to the vOTU DB
PROFILE=$2 # path to the merged profile
metadata=$3 # path to sample metadata file
abundance_table=$4 # path to abundance table
vOTU_persistence_mother=$5 # path to file with vOTU persistence per mother
DGR_info_table=$6 # path to file with DGR information
list_mothers=$7 # list of mothers

MOTHER_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $list_mothers)

echo '-------------------- WORKING WITH' ${MOTHER_ID} '--------------------'

# Clean environment, load modules
ml purge; ml Anaconda3; ml list
conda activate /scratch/hb-llnext/conda_envs/anvio-8; conda list

mkdir -p 2_VARIABILITY 2_VARIABILITY/${MOTHER_ID} OUTPUT_files
mkdir -p ${TMPDIR}/${MOTHER_ID}

# Select persistent DGR+ vOTUs in the mother
persistent_vOTUs=$(awk -v id="$MOTHER_ID" '$1 == id {print $3}' "$vOTU_persistence_mother" | tr ',' '\n')
num_vOTUs=$(echo "$persistent_vOTUs" | wc -l)
echo
echo "Mother $MOTHER_ID has $num_vOTUs persistent vOTUs"
echo

# Loop over persistent vOTUs for this mother
for votu in $persistent_vOTUs; do
    echo "---- Processing vOTU: $votu ----"

    # A: Select mother-specific samples with abundance > 0 (via R script)
    samples_file=${TMPDIR}/${MOTHER_ID}/samples.txt
    ./select_samples.R "$metadata" "$abundance_table" "$MOTHER_ID" "$votu" "$samples_file"

    # B: Extract DGR target gene(s) for this vOTU
    awk -v id="$votu" '$1 == id {print $12}' "$DGR_info_table" > ${TMPDIR}/${MOTHER_ID}/target_genes.txt
    echo "Target genes:"
    cat ${TMPDIR}/${MOTHER_ID}/target_genes.txt

    # C: Run anvi-gen-variability-profile
    outdir="2_VARIABILITY/${MOTHER_ID}/${votu}"
    mkdir -p $outdir

    # Nucleotide level
    anvi-gen-variability-profile -p $PROFILE \
                                 -c $vOTU_DB \
                                 --genes-of-interest ${TMPDIR}/${MOTHER_ID}/target_genes.txt \
                                 --samples-of-interest $samples_file \
                                 --engine NT \
                                 --quince-mode \
                                 -o $outdir/${votu}_NT_variability.txt
    # Amino acid level
    anvi-gen-variability-profile -p $PROFILE \
                                 -c $vOTU_DB \
                                 --genes-of-interest ${TMPDIR}/${MOTHER_ID}/target_genes.txt \
                                 --samples-of-interest $samples_file \
                                 --engine AA \
                                 --kiefl-mode \
                                 --include-additional-data \
                                 -o $outdir/${votu}_AA_variability.txt

done

