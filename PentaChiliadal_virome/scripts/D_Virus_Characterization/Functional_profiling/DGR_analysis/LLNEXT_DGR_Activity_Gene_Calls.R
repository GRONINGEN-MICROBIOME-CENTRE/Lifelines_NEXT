#****************
# Load libraries
#****************
library(stringr)
library(seqinr)

#****************
# Input arguments
#****************
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: script.R <prodigal_gv_output.csv> <prodigal_proteins.faa>")
}

input_file <- args[1]   # Prodigal-gv .csv file
faa_file   <- args[2]   # Prodigal-gv protein fasta
output_file <- "external_gene_calls.txt"

#***********************************************
# 1. Parse prodigal CSV to get gene call table
#***********************************************
lines <- readLines(input_file)

gene_calls <- data.frame()
gene_id <- 1
current_contig <- NA

for (i in seq_along(lines)) {
  line <- lines[i]

  # Extract contig name
  if (grepl("^DEFINITION", line)) {
    match <- str_match(line, 'seqhdr=\\"([^\\"]+)\\"')
    current_contig <- match[2]

  # CDS line
  } else if (grepl("^\\s*CDS\\s", line)) {

    # Reverse strand case
    if (grepl("complement\\(", line)) {
      coords <- str_match(line, "complement\\(<?(\\d+)\\.\\.>?(\\d+)\\)")
      start <- as.numeric(coords[2]) - 1   # <-- Convert to 0-based
      stop  <- as.numeric(coords[3])
      direction <- "r"
    } else {
      coords <- str_match(line, "\\s*<?(\\d+)\\.\\.>?(\\d+)")
      start <- as.numeric(coords[2]) - 1   # <-- Convert to 0-based
      stop  <- as.numeric(coords[3])
      direction <- "f"
    }

    # Extract partial info from next line (/note=...)
    note_line <- lines[i + 1]
    partial_match <- str_match(note_line, "partial=([01]{2})")
    partial_value_raw <- if (!is.na(partial_match[2])) partial_match[2] else "00"
    partial_value <- if (partial_value_raw == "00") 0 else 1

    # Append to table
    gene_calls <- rbind(
      gene_calls,
      data.frame(
        gene_callers_id = gene_id,
        contig = current_contig,
        start = start,
        stop = stop,
        direction = direction,
        partial = partial_value,
        call_type = 1,
        source = "Prodigal",
        version = "2.11.0",
        stringsAsFactors = FALSE
      )
    )
    gene_id <- gene_id + 1
  }
}

#***********************************************
# 2. Parse amino acid fasta and attach sequences
#***********************************************
faa <- read.fasta(faa_file, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
if (length(faa) != nrow(gene_calls)) {
  stop("Number of sequences in FAA does not match number of CDS in CSV!")
}

# Remove trailing asterisk (stop codon)
gene_calls$aa_sequence <- sub("\\*$", "", unlist(faa))

#***********************************************
# 3. Write output
#***********************************************
write.table(
  gene_calls,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Finished writing", nrow(gene_calls), "gene calls to", output_file, "\n")

