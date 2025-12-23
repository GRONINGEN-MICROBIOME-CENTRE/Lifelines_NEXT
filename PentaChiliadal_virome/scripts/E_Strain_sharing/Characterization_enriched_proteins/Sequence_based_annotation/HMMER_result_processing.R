################################################################################
##### Viral protein annotation: PHROGS HMM  results processing
### Author(s):Asier Fernández-Pato
### Last updated: 2nd October, 2025
################################################################################

# Define the directory paths as arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if four arguments are provided
if (length(args) != 3) {
  stop("Usage: Metadata_annotation_DBs Hmmsearch_results_PHROGs output_dir")
}

#****************
# Load libraries
#****************
library(data.table)
library(dplyr)

# Load metadata and HMMER result)
metadata <- read.delim(args[1])
HMMER_result <- fread(args[2])
output_dir <- args[3]

# Combine all hmmsearch results in one table
colnames(HMMER_result) [colnames(HMMER_result) == "Description"] <- "Protein_Info"

# Combine HMMER result and metadata to get the annotations for each protein (order by E-value)
Protein_annotations <- left_join(HMMER_result, metadata,  by = c("Match" = "NAME"))
Protein_annotations <- Protein_annotations[order(Protein_annotations$E_value), ]
Protein_annotations_simple <- Protein_annotations[, c("Protein_ID","Match", "Protein_Info", "DESCRIPTION", "CATEGORY")]

# Get unninformative annotations and filter them out
Unique_annotations <- data.frame(sort(table(Protein_annotations_simple$DESCRIPTION), decreasing = T))
Unninformative_annotations <- grep("unknown|hypothetical", Unique_annotations$Var1, ignore.case = TRUE, value = TRUE) 
Protein_annotations_simple_informative <- Protein_annotations_simple[!Protein_annotations_simple$DESCRIPTION %in% (Unninformative_annotations),]

#************************************************************
# Generate the annotation table with one row per protein
#************************************************************
# A) Raw Annotation table: includes all the annotations from all the HMM matches (including NAs and unninformative)
# Therefore, each entry in DESCRIPTION or CATEGORY corresponds to the HMM match in the same order
Protein_functional_annotation <- Protein_annotations_simple %>%
  group_by(Protein_ID) %>%
  summarize(
    Match = paste(Match, collapse = ";"),
    Protein_Info = first(Protein_Info),
    DESCRIPTION = paste(DESCRIPTION, collapse = "; "),
    CATEGORY = paste(CATEGORY, collapse = "; "),
  )

# B) Clean Annotation table: includes the processed annotations (condensing repeated annotations and removing NAs)
# Here, unninformative annotations are kept (some proteins will only have these annotations)
Protein_functional_annotation_clean <- Protein_annotations_simple %>%
  group_by(Protein_ID) %>%
  summarize(
    Match = paste(Match, collapse = ";"),
    Protein_Info = first(Protein_Info),
    DESCRIPTION = paste(na.omit(unique(DESCRIPTION)), collapse = "; "),
    CATEGORY = paste(na.omit(unique(CATEGORY)), collapse = "; ") 
  ) %>%
  filter(DESCRIPTION != "")  

# C) Clean and Informative Annotation table: includes the processed annotations (removing unninformative annotations)
Protein_functional_annotation_clean_informative <- Protein_annotations_simple_informative %>%
  group_by(Protein_ID) %>%
  summarize(
    Match = paste(Match, collapse = ";"),
    Protein_Info = first(Protein_Info),
    DESCRIPTION = paste(na.omit(unique(DESCRIPTION)), collapse = "; "),
    CATEGORY = paste(na.omit(unique(CATEGORY)), collapse = "; ")
  ) %>%
  filter(DESCRIPTION != "")  

# Generate list of annotated proteins
annotated_proteins <- Protein_functional_annotation_clean_informative$Protein_ID

#************************************************************
# Save output files
#************************************************************
write.table(Protein_functional_annotation,
            file = file.path(output_dir, "Protein_functional_annotation.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Protein_functional_annotation_clean,
            file = file.path(output_dir, "Protein_functional_annotation_clean.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(Protein_functional_annotation_clean_informative,
            file = file.path(output_dir, "Protein_functional_annotation_clean_informative.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)   
write.table(annotated_proteins,
            file = file.path(output_dir, "Annotated_proteins.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)   
