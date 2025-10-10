###     SCRIPT: MICROBIAL TAXA BM and VG SUMMARY STATISTICS
###     AUTHOR(S): TRISHLA SINHA
###     DESCRIPTION: SUMMARY STATISTICS MICROBIAL TAXA
###     PROJECT: LL-NEXT
###     LAST UPDATED: 20TH JUNE 2025

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/taxa/")

#libraries
library(tidyverse)

#Contents 
## 0. LOAD FILES
## 1. CROSS-SECTIONAL SUBSETS
## 2. LONGITUDINAL SUBSETS

#functions
#summary statistics 
generate_summary_statistics <- function(df) {
  
  #subset numeric variables
  numeric_vars <- df %>%
    select_if(is.numeric)
  
  #generate summary statistics for numeric variables
  numeric_summary <- numeric_vars %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      Total_N = sum(!is.na(value)),
      percentage_of_total_n = sum(!is.na(value)) / nrow(df) * 100,
      total_number_missings = sum(is.na(value)),
      percentage_of_total_number_missings = sum(is.na(value)) / nrow(df) * 100,
      Mean = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Min = min(value, na.rm = TRUE),
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      Max = max(value, na.rm = TRUE)
      
    ) %>%
    ungroup()
  
  #subset categorical variables
  factor_vars <- df %>%
    select_if(is.factor)
  
  #generate summary statistics for categorical variables
  factor_summary <- factor_vars %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      n_factor_levels = n_distinct(value, na.rm = TRUE),
      total_n = sum(!is.na(value)),
      percentage_total_n = sum(!is.na(value)) / nrow(df) * 100,
      total_number_missings = sum(is.na(value)),
      percentage_of_total_number_missings = sum(is.na(value)) / nrow(df) * 100)
  factor_summary_levels <- #this script was updated here for factor_summary_levels and how it was defined
    factor_vars %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    group_by(variable, value) %>%
    summarise(
      n = n(),
      perc = n() / nrow(factor_vars) * 100
    ) %>%
    drop_na()
  factor_summary_levels <- factor_summary_levels %>%
    mutate(perc=n/sum(n)*100) %>%
    mutate(
      level = paste0("level_", row_number())
    ) %>%
    filter(value!="NA") %>%
    pivot_wider(names_from = level,
                values_from = c(value, n, perc),
                names_glue = "{level}_{.value}")
  
  #determine the maximum number of levels in the factor_summary_levels (for colname reordering)
  level_cols <- grep("^level_\\d+_value$", names(factor_summary_levels), value = TRUE)
  n <- max(as.integer(gsub("^level_(\\d+)_value$", "\\1", level_cols)))
  
  #define order of colnames for categorical variables
  colnames_order <- c("variable",
                      unlist(lapply(1:n, function(i) paste0("level_", i, "_", c("value", "n", "perc")))))
  
  #reorder the columns of the dataframe
  factor_summary_levels <- factor_summary_levels[, colnames_order]
  
  #combine summary statistics data for factor summary
  factor_summary <- left_join(factor_summary,factor_summary_levels, by = "variable") %>%
    ungroup()
  
  #export the new table
  write.table (numeric_summary, file = "./meta_data_numeric_summary_stats.txt" , quote = F, sep = "\t", row.names = F)
}

##### =========================== 0. LOAD FILES =========================== #####

#metadata associated with microbial data
metadata <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/Metadata_BM_VG_Gut_20_05_2025.txt")
str(metadata, list.len=ncol(metadata))
summary(metadata)

#microbial taxa

taxa_all<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/taxa/ALL_METAPHLAN_BM_VG_GUT_metaphlan_12_05_2025.txt")
row.names(taxa_all)<-taxa_all$clade_name
taxa_all$clade_name=NULL
taxa_all<-as.data.frame(t(taxa_all))
taxa_all$NG_ID <-row.names(taxa_all)
taxa_all$NG_ID <- str_sub(taxa_all$NG_ID, 1, 12)
row.names(taxa_all)<-taxa_all$NG_ID
taxa_all$NG_ID=NULL
taxa=taxa_all
colnames(taxa) <- gsub("\\|", ".", colnames(taxa))


##### =========================== 1. Breast milk =========================== #####
####1.1 mother cross-sectional summary statistics####
metadata_BM <- metadata %>% filter (
  Type == "Milk")
#91  10

metadata_BM$NEXT_ID <- as.factor(as.character(metadata_BM$NEXT_ID))
str(metadata_BM) #91 unique mothers
 

taxa_BM <- taxa %>%
  filter(row.names(.) %in% metadata_BM$NG_ID)
taxa_BM <- taxa_BM %>%
  #select(where(~sum(.x) != 0)) %>% #retrospectively omitted this, to be able to compare the presence/absence of these microbes between mothers and infants (original script: remove columns (i.e. taxa) where values for all participants == 0; probably present in infants)
  mutate(NG_ID = as.factor(rownames(taxa_BM))) %>%
  relocate(NG_ID, .before=UNCLASSIFIED)
str(taxa_BM)


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/New_supplementary_files")


#generate summary statistics
generate_summary_statistics(taxa_BM) 

#generate prevalence
taxa_BM_prev <- taxa_BM %>%
  summarise(across(-NG_ID, ~ sum(. > 0) / n())) %>%
  pivot_longer(cols = everything(),
               names_to = "taxon_full_name",
               values_to = "Prevalence")

#load and format summary statistics 
taxa_BM_sumstats <- read.delim("meta_data_BM_numeric_summary_stats.txt", sep = "\t", header = TRUE, stringsAsFactors = T)
str(taxa_BM_sumstats, list.len=ncol(taxa_BM_sumstats))
#6938 obs. of  12 variables

taxa_BM_sumstats <- taxa_BM_sumstats %>% 
  #select(-c(total_n, percentage_of_total_n,
  #          total_number_missings, percentage_of_total_number_missings)) %>%
  rename(taxon_full_name=variable) %>%
  mutate(taxon_full_name = as.character(taxon_full_name))%>%
  mutate(taxon_short_name = ifelse(grepl("\\.", taxon_full_name),
                                   str_extract(taxon_full_name, "(?<=\\.)[^.]*$"),
                                   taxon_full_name)) %>%
  relocate(taxon_short_name, .after = taxon_full_name) %>%
  mutate(taxonomic_level = case_when(
    startsWith(taxon_short_name, "k__") ~ "kingdom",
    startsWith(taxon_short_name, "p__") ~ "phylum",
    startsWith(taxon_short_name, "c__") ~ "class",
    startsWith(taxon_short_name, "o__") ~ "order",
    startsWith(taxon_short_name, "f__") ~ "family",
    startsWith(taxon_short_name, "g__") ~ "genus",
    startsWith(taxon_short_name, "s__") ~ "species",
    startsWith(taxon_short_name, "t__") ~ "SGB",
    TRUE ~ "unclassified")) %>%
  relocate(taxonomic_level, .after=taxon_short_name) %>%
  mutate_if(is.character, as.factor)

summary(taxa_BM_sumstats)

#join with summary statistics with prevalence data
taxa_BM_sumstats_prev <- full_join(taxa_BM_sumstats, taxa_BM_prev) #Joining, by = "taxon_full_name"


taxa_BM_sumstats_prev$type <- "Breast milk"
taxa_BM_sumstats_prev  <- taxa_BM_sumstats_prev  %>% 
  relocate(type, .after=taxonomic_level)

write.table(taxa_BM_sumstats_prev, "supplementary_table_BM_taxa_summary_stats.txt", sep = "\t", row.names = F)




##### =========================== 2. Vaginal ========================== #####
####1.1 mother cross-sectional summary statistics####
metadata_VG <- metadata %>% filter (
  Type == "VG")
#91  10

metadata_VG$NEXT_ID <- as.factor(as.character(metadata_VG$NEXT_ID))
str(metadata_VG) #91 unique mothers


taxa_VG <- taxa %>%
  filter(row.names(.) %in% metadata_VG$NG_ID)
taxa_VG <- taxa_VG %>%
  #select(where(~sum(.x) != 0)) %>% #retrospectively omitted this, to be able to compare the presence/absence of these microbes between mothers and infants (original script: remove columns (i.e. taxa) where values for all participants == 0; probably present in infants)
  mutate(NG_ID = as.factor(rownames(taxa_VG))) %>%
  relocate(NG_ID, .before=UNCLASSIFIED)
str(taxa_VG)


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/New_supplementary_files")


#generate summary statistics
generate_summary_statistics(taxa_VG) 

#generate prevalence
taxa_VG_prev <- taxa_VG %>%
  summarise(across(-NG_ID, ~ sum(. > 0) / n())) %>%
  pivot_longer(cols = everything(),
               names_to = "taxon_full_name",
               values_to = "Prevalence")

#load and format summary statistics 
taxa_VG_sumstats <- read.delim("meta_data_VG_numeric_summary_stats.txt", sep = "\t", header = TRUE, stringsAsFactors = T)
str(taxa_VG_sumstats, list.len=ncol(taxa_VG_sumstats))
#6938 obs. of  12 variables

taxa_VG_sumstats <- taxa_VG_sumstats %>% 
  #select(-c(total_n, percentage_of_total_n,
  #          total_number_missings, percentage_of_total_number_missings)) %>%
  rename(taxon_full_name=variable) %>%
  mutate(taxon_full_name = as.character(taxon_full_name))%>%
  mutate(taxon_short_name = ifelse(grepl("\\.", taxon_full_name),
                                   str_extract(taxon_full_name, "(?<=\\.)[^.]*$"),
                                   taxon_full_name)) %>%
  relocate(taxon_short_name, .after = taxon_full_name) %>%
  mutate(taxonomic_level = case_when(
    startsWith(taxon_short_name, "k__") ~ "kingdom",
    startsWith(taxon_short_name, "p__") ~ "phylum",
    startsWith(taxon_short_name, "c__") ~ "class",
    startsWith(taxon_short_name, "o__") ~ "order",
    startsWith(taxon_short_name, "f__") ~ "family",
    startsWith(taxon_short_name, "g__") ~ "genus",
    startsWith(taxon_short_name, "s__") ~ "species",
    startsWith(taxon_short_name, "t__") ~ "SGB",
    TRUE ~ "unclassified")) %>%
  relocate(taxonomic_level, .after=taxon_short_name) %>%
  mutate_if(is.character, as.factor)

summary(taxa_VG_sumstats)

#join with summary statistics with prevalence data
taxa_VG_sumstats_prev <- full_join(taxa_VG_sumstats, taxa_VG_prev) #Joining, by = "taxon_full_name"


taxa_VG_sumstats_prev$type <- "Vaginal"
taxa_VG_sumstats_prev  <- taxa_VG_sumstats_prev  %>% 
  relocate(type, .after=taxonomic_level)

write.table(taxa_VG_sumstats_prev, "supplementary_table_VG_taxa_summary_stats.txt", sep = "\t", row.names = F)

