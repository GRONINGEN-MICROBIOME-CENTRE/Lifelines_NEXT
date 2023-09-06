Motivation: Chiliadal virome project used [xGen ssDNA & low-input library preparation kit](https://eu.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-library-preparation/dna-library-preparation/ssdna-low-input-dna-library-prep-kit)  for sequencing virus DNA.
It enabled us to quantitatively annotate dsDNA, ssDNA, RNA virome, but also imposed an issue of read duplication and uneven coverage of assembled contigs.
Based on the paper by [Roux et al., PeerJ, 2019](http://dx.doi.org/10.7717/peerj.6902), a pre-assembly read quality control was chosen for these samples,
but since [Clumpify guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/) from [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/) implied the unpaired read input, a few tests were performed to find the optimal settings.

Here are the scripts for experiments of different settings for:
- read deduplication step (```02.reads_QC_dedup.EXP*.sh```)
- read assembly (input, ```03.assembly.EXP*.sh```)

```runAllSamples.bash``` is a wrapper around sbatching the jobs

**The conclusion from the experiments:**
- read deduplication should be run in paired-end read (does not matter if reads are supplemented in separate of interleaved format)
- supplementing unmatched reads slightly improves the assembly characteristics 
