#!/bin/bash

for SAMPLE in $@; do
 	sbatch --output ./out/${SAMPLE}_03.out --error ./err/${SAMPLE}_03.err --job-name vQC_${SAMPLE} 03.virus_enrichment_index.sh ${SAMPLE}
done
