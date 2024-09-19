# Combined Analysis of Viral Contigs from Four Studies

This repository contains scripts and instructions for the combined analysis of viral contigs from four different studies. The general workflow consists of the following steps:

## Table of Contents
1. [Initial Dereplication of Viral Contigs](#1-initial-dereplication-of-viral-contigs)
2. [Building BOWTIE2 Index](#2-building-bowtie2-index)
3. [Parsing Read Mapping Results](#3-parsing-read-mapping-results)
4. [Initial RPKM Creation](#4-initial-rpkm-creation)
5. [Host Prediction](#5-host-prediction)
6. [Dereplication for Decontamination (99% ANI and 85% AF)](#6-dereplication-for-decontamination-99-ani-and-85-af)
7. [Dereplication of Decontaminated Contigs](#7-dereplication-of-decontaminated-contigs)
8. [Building BOWTIE2 Index for Dereplicated Decontaminated Contigs](#8-building-bowtie2-index-for-dereplicated-decontaminated-contigs)
9. [Parsing Read Mapping Results (Decontaminated Contigs)](#9-parsing-read-mapping-results-decontaminated-contigs)
10. [Decontaminated RPKM Creation](#10-decontaminated-rpkm-creation)
    - [RPKM Table with 95% Breadth of Coverage](#a-rpkm-table-with-95-breadth-of-coverage)
    - [RPKM Table with 99% Breadth of Coverage](#b-rpkm-table-with-99-breadth-of-coverage)
    - [Final Decontaminated RPKM Table](#c-final-decontaminated-rpkm-table)

## Step-by-Step Description

### 1. Initial Dereplication of Viral Contigs

**Script:** `07.pd_dereplication.sh`  
This step performs dereplication on viral contigs, obtained after the "Plasmid Removal, Prophage Pruning, and Quality Assessment" step, using 95% ANI and 85% AF. 

**Note:**  
- Fasta files with viral contigs from each study were manually merged before this step.
- Ensure to download `anicalc.py` and `aniclust.py` (links are provided in the script).

### 2. Building BOWTIE2 Index

**Script:** `08.pd_virusindex.sh`  
Builds a BOWTIE2 index for the dereplicated viral contigs.

### 3. Parsing Read Mapping Results

**Scripts:**  
- `10.pd_combine_coverage.sh`
- `combine_coverage.R` (launched by the `.sh` script)

This step combines the read mapping results from all samples across the studies to generate `coverage_table` and `count_table`.

**Note:**  
Run the read mapping (Step 9) for each study separately before executing this step.

### 4. Initial RPKM Creation

**Scripts:**  
- `11.01a.get_RPKM.sh`
- `01a.get_RPKM.R` (launched by the `.sh` script)

Creates the initial RPKM table based on 75% breadth of coverage.

### 5. Host Prediction

**Scripts:**  
- `12.a.contigs_split.sh` (splits contigs for performance reasons)
- `12.b.host_prediction.sh`

Predicts the host of viral contigs.

### 6. Dereplication for Decontamination (99% ANI and 85% AF)

**Scripts:**  
- `07.pd_dereplication_99.sh`
- `Extracting_contics_fom_NCs_VCs.R`

This step performs dereplication at the strain level (99% ANI and 85% AF).  
**Note:**  
Use `SeqKit` to manually remove contigs listed in the output from `Extracting_contics_fom_NCs_VCs.R`.

### 7. Dereplication of Decontaminated Contigs

**Script:** `07.pd_dereplication_decontam_99.sh`  
Dereplication of decontaminated viral contigs using 95% ANI and 85% AF.

### 8. Building BOWTIE2 Index for Dereplicated Decontaminated Contigs

**Script:** `08.pd_virusindex_decontam_99.sh`  
Builds a BOWTIE2 index for dereplicated, decontaminated viral contigs.

### 9. Parsing Read Mapping Results (Decontaminated Contigs)

**Scripts:**  
- `10.pd_combine_coverage_decontam_99.sh`
- `combine_coverage.R` (reused from previous step with different paths)

This step parses the read mapping results from all samples to the dereplicated, decontaminated viral contigs to create the `coverage_table` and `count_table`.

**Note:**  
Run the read mapping separately for each study before executing this step.

### 10. Decontaminated RPKM Creation

This process consists of three main steps:

#### a. RPKM Table with 95% Breadth of Coverage
**Scripts:**  
- `02.get_RPKM_decontam_99.R`
- `13.02.decontam_99_get_RPKM.sh`

Generates the RPKM table with 95% breadth of coverage.

#### b. RPKM Table with 99% Breadth of Coverage
**Scripts:**  
- `02.get_RPKM_decontam_99_99_coverage.R`
- `13.02.decontam_99_get_RPKM_99_coverage.sh`

Generates the RPKM table with 99% breadth of coverage.

#### c. Final Decontaminated RPKM Table
**Script:** `decontaminating_RPKMs.R`

- Sets values to 0 for contigs detected simultaneously in negative controls and samples (based on 99% breadth of coverage).
- Creates the final decontaminated RPKM table (based on 95% breadth).
- Generates an RPKM table decontaminated at the species level.

## Additional Notes

- **R Scripts**: The `.R` scripts handle specific parts of the analysis and are called within shell scripts where necessary.

---
