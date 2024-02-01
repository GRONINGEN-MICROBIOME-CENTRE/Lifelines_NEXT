Directories contain attempts at code optimization and often contain their own README (this feature will probably be removed to avoid confusion).

The scripts outside directories are the actual scripts used for the upstream data analysis and are enumerated in the order of execution (the same numbers are used for scripts that can be executed in parallel). 

To execute R scripts in the upstream analysis:

```
# DIR: /scratch/p282752/ANALYSIS_CHILIADAL/scripts
# DIR exists: /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/table_of_origin
module load R
Rscript table_of_origin.R
```
