# PentaChiliadal-Virome
Gut virome analysis from metagenomic data on Lifelines NEXT Cohort.

## Content

We here describe the methods used to:

- Perform the viral detection from assembled metagenomic data (VirSorter2, DeepVirFinder, geNomad)
- Extend predicted viral ccontigs (COBRA)
- Do an additional filtering of viral sequences and trim bacterial regions from prophages (geNomad)
- Further remove host contamination from prophages and filter viral sequences with low completeness (≤50%) (CheckV)
- Dereplicate the predicted viral sequences A) by themselves and B) with viral genomes from public databases
- Perform the abundance estimation of the viral sequences
- Downstream analyses 
