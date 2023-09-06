#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=12gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge

echo -e "Sample\tReads\tReads_HQ\tSSU rRNA alignment rate\tLSU rRNA alignment rate\tBacterial_Markers alignment rate\ttotal enrichmnet score" >> 08_virus_enrichment.txt

for SAMPLE_ID in $@; do 
	tail -1 ../${SAMPLE_ID}/viromeqc_stat.txt >> tmp
	less tmp | paste -s >> tmp2
	paste -d "\n" tmp2 >> 08_virus_enrichment.txt
	rm tmp*
done

module list

module purge

