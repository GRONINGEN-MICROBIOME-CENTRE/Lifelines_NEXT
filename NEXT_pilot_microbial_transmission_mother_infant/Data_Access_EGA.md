# Instructions on sequencing data and phenotypes access from EGA

To access sample information, basic phenotypes, family structure, and quality-trimmed sequencing reads from EGA, you need to follow the instructions provided below:

## Prerequisites
1. Read the [LLNEXT DATA ACCESS AGREEMENT](https://groningenmicrobiome.org/?page_id=2598). Access to the LLNEXT Project data will be granted 
to all qualified researchers and will be governed by the provisions laid out in this data access agreement.

2. Fill out the [Application Form for Access to Lifelines NEXT Data](https://docs.google.com/forms/d/e/1FAIpQLScUaLZk6Smz66EAqgb0JmzyXLPF3V9mHdvWEuL98qT4yF1j5g/viewform) to indicate that you agree with LLNEXT DATA ACCESS AGREEMENT.

3. To request access to the data, [register an account at EGA](https://ega-archive.org/register/). Validate your account. The validation time for new accounts by EGA varies, but it usually takes around 18-24 hours.

## Requesting data via EGA
1. Login to EGA
   
2. Request access to datasets from [LifeLines-NEXT pilot study](https://ega-archive.org/studies/EGAS00001005969). To do so, press the 'Request Access' button on the page of the dataset of interest ([LLNEXT Pilot MGS sequencing](https://ega-archive.org/datasets/EGAD00001011293) or [LLNEXT Pilot VLP sequencing](https://ega-archive.org/datasets/EGAD00001011291), or both).

3. After receiving your filled **Application Form for Access to Lifelines NEXT Data** and your **Access Request**, the Data Access Committee will evaluate your request and grant access to the data.

### Metadata (phenotypes):
Upon receiving access to a dataset, you will be able to download the metadata containing phenotypes. Below is an example of parsing these files.

```
# necessary libraries:
library(jsonlite)
library(data.table)

# importing the metadata:
DF <- read.table("samples.tsv", header=T, sep='\t', quote = "")

# Extracting the JSON-like strings from the DataFrame
json_strings <- as.character(DF[, "extra_attributes"])

# Parsing each JSON-like string
parsed_json <- lapply(json_strings, function(x) jsonlite::fromJSON(x))

# Reformatting each data frame
result_list <- lapply(parsed_json, function(json_obj) {
  buffer <- as.data.frame(t(as.data.table(json_obj)))
  colnames(buffer) <- buffer[1, ]
  buffer <- buffer[3, , drop = FALSE]  # Keeping it as a data frame
  return(buffer)
})

# Extracted values and columns
combined_df <- do.call(rbind, result_list)

# Resulting data frame
result_df <- cbind(DF, combined_df)

# Remove row names
row.names(result_df) <- NULL

```

(Sana to finish the parsing instructions and following steps)

### Sequencing data:

To download the sequencing data, you need to use [EGA download client: pyEGA3](https://github.com/EGA-archive/ega-download-client). The pyEGA3 has an extensive README.md on installation and usage. Below are details on installation and usage at the HPC of the University of Groningen:

### Installation of pyEGA3:

```
git clone https://github.com/EGA-archive/ega-download-client.git

module load Python
# Python version: 3.11.3

cd ega-download-client/
# edit red_hat_dependency_install.sh in the following way:
# remove all "sudo"
sh red_hat_dependency_install.sh

# testing installation:
python -m pyega3.pyega3 --help

wget https://raw.githubusercontent.com/EGA-archive/ega-download-client/master/pyega3/config/default_credential_file.json

# edit default_credential_file.json in the following way:
# "username": "e-mail used for EGA account registration",
# "password": "password to the EGA account"

mv default_credential_file.json ega-download-client/credential_file.json

# to check which datasets are available:

python -m pyega3.pyega3 datasets

```

