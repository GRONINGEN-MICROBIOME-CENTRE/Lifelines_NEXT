# General read quality control and assembly pipeline
Gut virome analysis from metagenomic data on Lifelines NEXT Cohort.

## Content

We here describe the methods used to perform:

- Quality control of the raw metagenomic reads (FastQC)
- Removal of adapter sequences (BBDuk)
- Removal of contaminant reads (KneadData)
- Assembly of quality-filtered and trimmed reads (metaSPAdes)
- Assessment of assembly quality (QUAST)
