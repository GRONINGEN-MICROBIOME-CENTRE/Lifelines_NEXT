#!/bin/bash

for SAMPLE in $@; do
 	# Running the regular assembly script:
	sbatch --output ./out/02.rAs/${SAMPLE}_02_se_s1test.out --error ./err/02.rAs/${SAMPLE}_02_se_s1test.err --job-name RAs_se_s1test${SAMPLE} 02.sc_assembly_se_s1test.sh ${SAMPLE}
	# Running the time-extended assembly script:
	#sbatch --output ./out/${SAMPLE}_02.out --error ./err/${SAMPLE}_02.err --job-name RAs_${SAMPLE} 02a.sc_assembly.sh ${SAMPLE}
done
