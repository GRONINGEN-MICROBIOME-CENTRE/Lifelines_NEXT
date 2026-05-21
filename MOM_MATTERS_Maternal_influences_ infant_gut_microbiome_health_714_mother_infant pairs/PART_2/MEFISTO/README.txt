This directory contains all code used run the MEFISTO factor analysis on LLNEXT data.

First, the <stability_analysis> directory contains the code used to select a stable number of MEFISTO factors, following the approach described in the manuscript. Individuals were repeatedly subsampled by selecting 80% of individuals without replacement across 15 subsampling runs. MEFISTO models were then compared across candidate factor numbers, and factor stability was assessed by correlating matched factor weights and factor scores across subsampling runs.

Second, the <three_factor_model> directory contains the code used to run MEFISTO on the complete dataset, including all individuals. Full-data models were run for selected numbers of factors, including 3, 6, 9, and 12 factors. Based on the stability analysis, the three-factor model was selected as the most stable model and was used for publication and for downstream association analyses between factors and phenotypes.

Third, the <penalized_regression> directory contains the code used to associate phenotypes with the MEFISTO factors using penalized regression. The final three-factor MEFISTO model was used to extract factor values, which were then modeled as outcomes in lasso-penalized mixed models with phenotype and covariate predictors.
