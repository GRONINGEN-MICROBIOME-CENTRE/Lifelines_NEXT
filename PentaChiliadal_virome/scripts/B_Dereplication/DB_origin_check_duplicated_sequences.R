# Load R packages
library(stats)
library(Biostrings)
library(data.table)

# Get FASTA file of viral genomes from command line arguments
args <- commandArgs(trailingOnly = TRUE)
viral_sequences <- readDNAStringSet(args[1]) # FASTA file with the viral sequences 

# Get sequence IDs
viral_DB_origin <- data.frame(names(viral_sequences))

# Add DB information
# Important things to consider:
# Yutin DB has some overlap with RefSeq and Benler DBs (not in ID, but in sequence)
# 8 duplicated sequences (in ID) between Yutin-Guerin DBs:ERR844030_ms_1, ERR975045_s_1, Fferm_ms_11, HvCF_A6_ms_4, IAS_virus_KJ003983
# Inf125_s_2, Sib1_ms_5, SRR4295175_ms_5
colnames(viral_DB_origin)[1] <- "viral_seq"
viral_DB_origin$DB <- NA
viral_DB_origin[,"DB"] [grep("^[ABCDY].*NODE", viral_DB_origin$viral_seq)] <- "LLNEXT" # 
viral_DB_origin[,"DB"] [grep("MGV-GENOME", viral_DB_origin$viral_seq)] <- "MGV" #189,680
viral_DB_origin[,"DB"] [grep("uvig|ivig", viral_DB_origin$viral_seq)] <- "GPD" # 82,621
viral_DB_origin[,"DB"] [grep("^NL", viral_DB_origin$viral_seq)] <- "Gulyaeva" #637
viral_DB_origin[,"DB"] [grep("^[UQPOC].*\\.1\\s[^>]*shotgun sequence", viral_DB_origin$viral_seq)] <- "Benler" #1480
viral_DB_origin[,"DB"] [grep("OTU_", viral_DB_origin$viral_seq)] <- "Shah" #4627
viral_DB_origin[,"DB"] [grep("IMGVR_", viral_DB_origin$viral_seq)] <- "IMG_VR" #36,064
viral_DB_origin[,"DB"] [grep("^VIBRANT|^VirFinder|^VirSorter2", viral_DB_origin$viral_seq)] <- "ELGV" #21,295
viral_DB_origin[,"DB"] [grep("^FSK", viral_DB_origin$viral_seq)] <- "Neg_control" #2
viral_DB_origin[,"DB"] [grep("^NC.*\\.|^AC.*\\.", viral_DB_origin$viral_seq)] <- "RefSeq" #5,199
viral_DB_origin[,"DB"][(nrow(viral_DB_origin)-248):nrow(viral_DB_origin)] <- "Guerin" #249
viral_DB_origin[,"DB"][which(is.na(viral_DB_origin$DB))] <- "Yutin" #656

# Check for duplicated values
virus_duplicated <- viral_DB_origin[duplicated(viral_DB_origin$viral_seq), "viral_seq"]
deduplicated_viral_DB_origin <- viral_DB_origin[!duplicated(viral_DB_origin$viral_seq), ]
deduplicated_viral_DB_origin$DB [which(deduplicated_viral_DB_origin$viral_seq %in% virus_duplicated)] <- "Guerin"
cat("The duplicated sequences are:", virus_duplicated)

# Add new simplified viral IDs (without prophage information from CheckV)
deduplicated_viral_DB_origin$new_viral_ID <- ave(deduplicated_viral_DB_origin$DB, deduplicated_viral_DB_origin$DB, 
                                                 FUN = function(x) paste0(x, "_", seq_along(x)))

# Remove sequences with duplicate IDs from the FASTA file
viral_sequences_no_dup <- viral_sequences[unique(names(viral_sequences))]

# Replace the IDs in the FASTA file with the new simplified IDs
names(viral_sequences_no_dup) <- deduplicated_viral_DB_origin$new_viral_ID

# Save final table with the X sequences used as input for STEP5 dereplication and their DB of origin (including 5 NEG-CONTROLS)
write.table(deduplicated_viral_DB_origin,"Dereplication_input_sequences_nodup_DB_origin.txt", sep = "\t", 
            row.names = F, col.names = F, quote = FALSE)

# Save the updated FASTA file
writeXStringSet(viral_sequences_no_dup, "Dereplication_combined_sequences_nodup_renamed.fa")
