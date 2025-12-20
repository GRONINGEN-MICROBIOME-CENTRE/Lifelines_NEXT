# Pregnancy and Early-Life Gut Virome in the LLNEXT cohort: Origin, Persistence, Influencing Factors and Health Implications

Gut virome analysis from metagenomic data on the **Lifelines NEXT Cohort**.

<p align="center">
  <img src="images/Project_overview.png" alt="PentaChiliadal Virome Workflow" width="650">
</p>

---

## Table of Contents
- [Overview](#overview)
- [Methods](#methods)
  - [Viral Detection](#viral-detection)
  - [Contig Extension & Filtering](#contig-extension--filtering)
  - [Quality Control & Host Removal](#quality-control--host-removal)
  - [Dereplication & Abundance Estimation](#dereplication--abundance-estimation)
- [Data Availability](#data-availability)
- [Code Structure](#code-structure)
- [Citation](#citation)

---

## Overview
The **PentaChiliadal Virome (PCV)** project aims to characterize the gut virome landscape from metagenomic data within the *Lifelines NEXT* cohort.  
We provide a reproducible pipeline for the detection, filtering, dereplication, and abundance estimation of viral sequences from metagenomic assemblies.

---

## Methods

### Viral Detection
We used multiple complementary tools to identify viral sequences from assembled contigs:
- **VirSorter2**, **DeepVirFinder**, and **geNomad** for initial viral prediction.

### Contig Extension & Filtering
- Predicted viral contigs were extended using **COBRA**.
- Bacterial regions in prophages were trimmed using **geNomad** filters.

### Quality Control & Host Removal
- **CheckV** was used to:
  - Estimate genome completeness.
  - Remove host contamination from prophages.
  - Exclude sequences with ≤50% completeness.

### Dereplication & Abundance Estimation
- Viral sequences were dereplicated:
  - (A) Internally within the dataset.
  - (B) Against public viral genome databases.
- Viral abundance was estimated by mapping reads back to the dereplicated viral genome set.

---

## Data Availability
All resulting viral genome sets, abundance tables, and associated metadata will be made available through:
- [Zenodo DOI link (coming soon)](https://zenodo.org/)
- Lifelines (to be added)

---

## Code Structure
Scripts and workflows used in this project are organized as follows:
