#!/bin/bash

#SBATCH --job-name=HP_4
#SBATCH --output=./out/HP/HP_4_%A.out
#SBATCH --error=./err/HP/HP_4_%A.err
#SBATCH --time=00:39:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

OUTPUT_DIR=$1 #path to results (output) directory

echo "Results (output directory): $(realpath "$OUTPUT_DIR")"

module purge

for f in ${OUTPUT_DIR}/host_predictions/Date_and_version_batch*; do [ -f "$f" ] && echo "=== $f ===" && cat "$f" && echo; done > ${OUTPUT_DIR}/host_predictions/Date_and_version_default_db.log
awk '(NR == 1) || (FNR > 1)' ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90_batch* > ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90_default_db.csv
awk '(NR == 1) || (FNR > 1)' ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90_batch* > ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90_default_db.csv
awk 'FNR <= 10 && NR <= FNR { print; next } FNR > 10 { print }' ${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool_batch* > "${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool_default_db.csv"

for f in ${OUTPUT_DIR}/host_predictions/Date_and_version_enriched_db_batch*; do [ -f "$f" ] && echo "=== $f ===" && cat "$f" && echo; done > ${OUTPUT_DIR}/host_predictions/Date_and_version_enriched_db.log
awk '(NR == 1) || (FNR > 1)' ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90_enriched_db_batch* > ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90_enriched_db.csv
awk '(NR == 1) || (FNR > 1)' ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90_enriched_db_batch* > ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90_enriched_db.csv
awk 'FNR <= 10 && NR <= FNR { print; next } FNR > 10 { print }' ${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool_enriched_db_batch* > "${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool_enriched_db.csv"

rm ${OUTPUT_DIR}/host_predictions/Date_and_version_*batch*
rm ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90_*batch*
rm ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90_*batch*
rm ${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool_*batch*

for f in ${OUTPUT_DIR}/host_predictions/Date_and_version*; do [ -f "$f" ] && echo "=== $f ===" && cat "$f" && echo; done > ${OUTPUT_DIR}/host_predictions/Date_and_version.log
awk '(NR == 1) || (FNR > 1)' ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90* > ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90.csv
awk '(NR == 1) || (FNR > 1)' ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90* > ${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90.csv
awk 'FNR <= 10 && NR <= FNR { print; next } FNR > 10 { print }' ${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool* > "${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool.csv"
