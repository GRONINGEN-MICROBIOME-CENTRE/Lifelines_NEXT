# Strainphlan 4.0 analysis Lifelines NEXT 

Adapted from Biobakery (StrainPhlAn 4.0). 
https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4

Authors: Trishla Sinha
Description: The script shows how strain profiling was performed using Strainphlan 4 for all maternal and infant samples post QC, for all species.   
Languages: Bash and R.   

## Step 1: Reconstruct all species strains

https://github.com/GRONINGEN-MICROBIOME-CENTRE/gmc-mgs-pipeline/blob/main/GMH_pipe.py 

## Step 2: Profile the clades present in the samples (profileClades.sh)

```
#!/bin/bash

#SBATCH --mem=24gb
#SBATCH --time=0-07:59:59
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --job-name=SP4pr
#SBATCH --error=__SP4_profile.err
#SBATCH --output=__SP4_profile.out

# NOTES:
# script profiles all clades in the dataset in given folder ($1)
# puts results in the current folder!
# Adding --mutation_rates will give a mutation rates table for each of the alignes markers and a summary table for the concatenated MSA
# Removing the --print_clades only will actually run it 

# PARAMS
N=1 # --marker_in_n_samples
S=10 # --sample_with_n_markers
DB=/data/umcg-tifn/rgacesa/conda_biobakery4/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl 

# purge modules
module purge
# load conda
ml Miniconda3/4.8.3
# load conda env
source activate /data/umcg-tifn/rgacesa/conda_biobakery4
# run clade profiling
strainphlan -s *.pkl --database /data/umcg-tifn/rgacesa/conda_biobakery4/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --marker_in_n_samples ${N} --sample_with_n_markers ${S} --print_clades_only --phylophlan_mode accurate --output_dir . > strainphlan4_clades_${N}.txt

```
### Execution 

```
sbatch ./profileClades.sh 

```
Mon Mar 20 17:47:58 2023: Start StrainPhlAn 4.0.6 execution
Mon Mar 20 17:47:58 2023: Loading MetaPhlAn mpa_vJan21_CHOCOPhlAnSGB_202103 database...
Mon Mar 20 17:48:20 2023: Done.
Mon Mar 20 17:48:23 2023: Detecting clades...
Mon Mar 20 18:40:14 2023: Done.
Mon Mar 20 18:40:14 2023: Detected clades: 
Mon Mar 20 18:40:14 2023:       t__SGB10068: in 103 samples.
Mon Mar 20 18:40:14 2023:       t__SGB8007_group: in 86 samples.
Mon Mar 20 18:40:14 2023:       t__SGB6936: in 82 samples.
Mon Mar 20 18:40:14 2023:       t__SGB17248: in 77 samples.
Mon Mar 20 18:40:14 2023:       t__SGB6952: in 74 samples.
Mon Mar 20 18:40:14 2023:       t__SGB6939: in 73 samples.
Mon Mar 20 18:40:14 2023:       t__SGB9712_group: in 60 samples.
Mon Mar 20 18:40:14 2023:       t__SGB10120: in 58 samples.
Mon Mar 20 18:40:14 2023:       t__SGB10115: in 58 samples.
Mon Mar 20 18:40:14 2023:       t__SGB17247: in 56 samples.
Mon Mar 20 18:40:14 2023:       t__SGB7962: in 55 samples.


We next process this output file to select only the clade names

```
cat strainphlan4_clades_1.txt | grep t__ | cut -f 2 | cut -f 1 -d ':' > LLNEXT_sp_clades_names.txt

```

This will give us the names of each species found: 
t__SGB10068
t__SGB8007_group
t__SGB6936
t__SGB17248
t__SGB6952
t__SGB6939
t__SGB9712_group
t__SGB10120
t__SGB10115
t__SGB17247
t__S6_group
t__SGB2303


## Step 3: Build the multiple sequence alignment (sp4_runMarkerComparison.sh)


```
#!/bin/bash
#SBATCH --mem=40gb
#SBATCH --time=2-23:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# NOTES:
# > $1 is input folder
# > $2 is clade name

echo "Invoking runMarkerComparison.sh"
echo "CLs: ${1} ${2}"

# HELP / USE
if [[ $# -ne 2 ]]; 
    then 
    echo "ERROR: script requires two command line parameters:"
    echo " <input folder with pkl files> <clade name>"
    exit 2
fi

# PARAMS
# ===================
N=10 # --marker_in_n_samples
S=10 # --sample_with_n_markers 
MODE=accurate # {accurate,fast}
CONDA=/scratch/hb-tifn/condas/conda_biobakery4/
#CM= # clade markers

# purge modules
module purge

# load conda
ml Anaconda3/2022.05
# load conda env
source deactivate
source activate ${CONDA}

# prep results folder (where clade result goes)
mkdir ${2}
# run strainphlan for that clade
echo "strainphlan -s ${1}/*.pkl --output_dir ${1}/${2} --clade ${2} --marker_in_n_samples ${N} --sample_with_n_markers ${S} --nprocs 8 --phylophlan_mode ${MODE}"
strainphlan -s ${1}/*.pkl --database /scratch/hb-tifn/condas/conda_biobakery3/lib/python3.9/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --output_dir ${1}/${2} --clade ${2} --marker_in_n_samples ${N} --sample_with_n_markers ${S} --nprocs 8 --phylophlan_mode ${MODE} #--tmp ${OUT_TMP}

```

### Execution

```
for i in $(cat LLNEXT_VG_BM_GUT_sp_clades_names.txt); do sbatch sp4_runMarkerComparison.sh  /scratch/p280306/LLNEXT_BREASTMILK_VAGINAL/ALL_STRAINPHLAN_BM_VG_GUT $i; done
```
Where LLNEXT_VG_BM_GUT_sp_clades_names.txt contains a list of names of all (sub) species identified in the previous step or a list of names of specific species that you want to run 
This will perform MSA and create .tre files and .aln files for each of the (sub)species you feed it in 



## Step 4: Make distance matrix from MSA file (makeDistMatAll.sh)

`
