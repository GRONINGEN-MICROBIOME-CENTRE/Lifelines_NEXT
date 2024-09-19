
# Garmaeva et al. (Nature communications, 2024) Data Analysis Scripts

This repository contains scripts used for the analysis of sequencing data from the Garmaeva et al. (Nature communications, 2024) study.

## Table of Contents
1. [Initial QC for Sequencing Data](#initial-qc-for-sequencing-data)
2. [Assembly](#assembly)
3. [Virome Discovery](#virome-discovery)
4. [Contig Extension](#contig-extension)
5. [Plasmid Removal, Prophage Pruning, and Quality Assessment](#plasmid-removal-prophage-pruning-and-quality-assessment)
6. [Read Alignment and Mapping](#read-alignment-and-mapping)
7. [Auxiliary Scripts](#auxiliary-scripts)

---

## Initial QC for Sequencing Data
Scripts used for quality control (QC) of the sequencing data.

- **01.reads_QC.sh**  
  Script for QC of all samples.
  
- **01.reads_QC_altREC.sh**  
  QC script for the negative control (NC) sample due to `metaSPAdes` failing to perform read error correction for this sample.  
  *Note*: This script was only used for the NC, while the other samples used `01.reads_QC.sh`.

### Features:
- Adapter trimming
- Quality trimming
- Human read removal

*Since no MDA or DNA amplification was used prior to sonication, read deduplication was not performed in this study.*

---

## Assembly
Scripts for assembling the sequencing data.

- **02.sc_assembly.sh**  
  Assembly script for all samples.

- **02.sc_assembly_altREC.sh**  
  Assembly script for the NC sample, adapted to handle issues with `metaSPAdes`.

### Assembly Details:
- Performed using `metaSPAdes` in meta mode.

---

## Virome Discovery
Discovery of viral sequences from assembled contigs larger than 1 kb using several tools.

- **03.vd_deepvirfinder.sh**  
  Virome discovery using DeepVirFinder.

- **03.vd_virsorter.sh**  
  Virome discovery using VirSorter.  
  **Note**: Download `prodigal-gv` prior to running this script (link in the script).

- **03.vd_vibrant.sh**  
  Virome discovery using VIBRANT.

- **03.vd_genomad.sh**  
  Virome discovery using geNomad.

---

## Parsing Virome Discovery Results
- **04.vd_parsing.sh**  
  Script for parsing the results from the virome discovery tools.

---

## Contig Extension
Extension of contigs using COBRA to create longer contigs.

- **05.vd_cobra.sh**  
  Contig extension script for all samples.

- **05.vd_cobra_GNC.sh**  
  Contig extension script for the NC sample.

---

## Plasmid Removal, Prophage Pruning, and Quality Assessment
- **06.vd_prophage_pruning_fin.sh**  
  Script to remove plasmids (identified by geNomad), prune prophages, and perform quality assessment using Check-V.  
  **Note**: Automatically launches the R script `New_contigs_ID_and_metadata_fin.R`.

---

## Read Alignment and Mapping
### Read Alignment
- **07.pd_readalign.sh**  
  Script for aligning sample reads to the viral contigs.

- **07.pd_readalign_GNS.sh**  
  Script for aligning NC reads to the viral contigs.

### Read Mapping to Dereplicated Contigs
- **09.at_readmapping_NCP.sh**  
  Script for mapping reads to dereplicated viral contigs for samples.

- **09.at_readmapping_NCP_GNC.sh**  
  Script for mapping reads to dereplicated viral contigs for NC.

### Read Mapping to Decontaminated Contigs
- **09.at_readmapping_decontam_99_NCP.sh**  
  Script for mapping reads to dereplicated viral contigs (99% breadth of coverage) for samples.

- **09.at_readmapping_decontam_99_NCP_GNC.sh**  
  Script for mapping reads to dereplicated viral contigs (99% breadth of coverage) for NC.

---

## Auxiliary Scripts
- **runAllSamples_02.bash**  
  Wrapper script for running `01.reads_QC.sh` and `02.sc_assembly.sh` for all samples.

- **runAllSamples_02_altREC.bash**  
  Wrapper script for running `01.reads_QC_altREC.sh` and `02.sc_assembly_altREC.sh` for the NC sample.

- **New_contigs_ID_and_metadata_fin.R**  
  R script to create new unified IDs and contig metadata.  
  **Note**: This script is automatically launched when running `06.vd_prophage_pruning_fin.sh`.

- **filter_contigs.pl**  
  Script to filter out contigs shorter than 1 kb.  
  **Note**: This script is automatically launched within `02.sc_assembly.sh` and `02.sc_assembly_altREC.sh`.

- **table_of_origin.R**  
  R script to generate contig metadata after viral discovery.
