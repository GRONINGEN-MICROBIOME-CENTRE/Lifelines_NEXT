#!/bin/bash

#SBATCH --job-name=PP_4
#SBATCH --output=./out/PP/PP_4.out
#SBATCH --error=./err/PP/PP_4.err
#SBATCH --time=00:39:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

OUTPUT_DIR=$1 #path to results (output) directory

echo "Results (output directory): $(realpath "$OUTPUT_DIR")"

module purge

cat ${OUTPUT_DIR}/PP_PNT/batch_*.faa > "${OUTPUT_DIR}/High_quality_vOTUs_PNT.faa"
cat ${OUTPUT_DIR}/PP_PNT/batch_*.gff > "${OUTPUT_DIR}/High_quality_vOTUs_PNT.gff"

# cat ${OUTPUT_DIR}/PP_PGV/batch_*.faa > "${OUTPUT_DIR}/High_quality_vOTUs_PGV.faa"
# cat ${OUTPUT_DIR}/PP_PGV/batch_*.gff > "${OUTPUT_DIR}/High_quality_vOTUs_PGV.gff"

# sed 's/[#+]$//' -i "${OUTPUT_DIR}/High_quality_vOTUs_PNT.faa"

module load Python/3.13.1-GCCcore-14.2.0

module list

python clean_faa.py "${OUTPUT_DIR}/High_quality_vOTUs_PNT.faa" "${OUTPUT_DIR}/High_quality_vOTUs_PNT_ed.faa"
python clean_gff.py "${OUTPUT_DIR}/High_quality_vOTUs_PNT.gff" "${OUTPUT_DIR}/High_quality_vOTUs_PNT_ed.gff"
