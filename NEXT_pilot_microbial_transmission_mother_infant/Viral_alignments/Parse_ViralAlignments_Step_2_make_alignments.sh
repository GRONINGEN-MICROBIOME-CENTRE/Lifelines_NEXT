#!/bin/bash
#SBATCH --mem=5gb
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name="kalign"

### Bash script to make alignments of consensus sequences
### To the database of non-redundant viral contigs 
## Date: 21 Apr 2023
## By: Alex Kurilshikov
### Notes:
## requires kalign program
## sometimes it fails and I have to run kalign manually on my local computer

cd Consensus3
mkdir all_aln
for i in `ls`; do cd ${i}
        echo ${i}
        ../../kalign-1.04/kalign -i merged.fa -o ${i}.aln.fa
		cp ${i}.aln.fa ../all_aln/
        cd ../
done
cd ../Consensus2

for i in `ls`; do cd ${i}
        echo ${i}
        ../../kalign-1.04/kalign -i merged.fa -o ${i}.aln.fa
        cp ${i}.aln.fa ../all_aln/
		cd ../
done
cd ../
echo "Finished"