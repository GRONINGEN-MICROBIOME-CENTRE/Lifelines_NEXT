################### CAZYMES PERMANOVA BACTERIA ##########################
# Author: T.Sinha 
# Last updated: 1st September, 2025

library(tidyverse)
library(vegan)
library(foreach)



# Load functions

Call_adonis <- function(phenotypes, Phenotype_list, Distance, perm = 10000, cores = 8) {
  Distance <- as.matrix(Distance)
  adonis_results <- data.frame(matrix(ncol = 5, nrow = length(Phenotype_list)))
  colnames(adonis_results) <- c("Phenotype", "Df", "F", "R2", "p-value")
  
  Phenotype_list <- Phenotype_list[!is.na(Phenotype_list)]
  
  adon <- foreach(i = 1:length(Phenotype_list), .combine = rbind) %do% {
    Pheno <- Phenotype_list[i]
    print(Pheno)
    
    Values <- as.vector(as_vector(phenotypes[[Pheno]]))
    
    # Filter NAs
    r <- row.names(phenotypes)[!is.na(Values)]
    phenos2 <- phenotypes[r, , drop = FALSE]
    distmat_cleaned <- Distance[r, r]
    
    # Formula without covariates
    FR <- as.formula(paste0("distmat_cleaned ~ ", Pheno))
    
    # Run adonis2
    ad1 <- adonis2(FR, data = phenos2, permutations = perm, parallel = cores, na.action = na.fail)
    
    # Extract results (row 1 if only 1 term in model)
    adonis_results[i, ] <- c(Pheno, ad1$Df[1], ad1$F[1], ad1$R2[1], ad1$`Pr(>F)`[1])
  }
  
  return(adonis_results %>% drop_na())
}


##############################
# Infant 
##############################

# Function to prepare cazyme data 
Prepare_cazymes = function(File = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/NEXT_cayman_infant.txt", Prevalence_min = 0.4, Transform = "log" ){
  caz = read_tsv(File)
  caz %>% filter(! feature %in% c("total_reads", "filtered_reads")) -> caz
  #Get unique counts + ambigous counts normalized
  caz %>% select(ID, feature, combined_rpkm ) %>% spread(key = feature, value = combined_rpkm) -> caz
  #Change ID
  caz$ID = caz$ID 
  #Get summary
  summary_as_df <- function(x) {
    result <- summary(x)
    data.frame(stat = names(result), value = as.vector(result))
  }
  summary_df <- caz %>%
    select(-ID) %>%
    purrr::map(summary_as_df) %>%
    bind_rows(.id = "Cazyme") %>% spread(stat, value) %>% as_tibble()
  summary_df = summary_df %>% mutate( `NA's` = ifelse( is.na(`NA's`), 0,  `NA's` ), N = dim(caz)[1], Present= N - `NA's`)
  summary_df = summary_df %>% mutate(Present_perc = Present/N ) %>% mutate(Keep = ifelse(Present_perc>=Prevalence_min, T, F ) )
  #Make NA into 0
  caz %>% mutate_all(~ replace(., is.na(.), 0)) -> caz
  #Remove Cazymes with high NA prop
  caz %>% select(c("ID", summary_df$Cazyme[summary_df$Keep==T] )) -> caz_anal
  if (Transform == "log"){
    Log_transform(caz_anal)  ->caz_anal_tranf
  } else if (Transform == "clr"){
    CLR_transformation(caz_anal)  -> caz_anal_tranf
  } else { caz_anal_tranf = NULL}
  
  return(list("caz" = caz, "caz_filtered" = caz_anal, "caz_tranf"=caz_anal_tranf  ,"Summary_caz" = summary_df ) )
}


Log_transform = function(df){
  colnames(df)[1] = "ID"
  PS = Pseudocount(df)
  df %>% select(-ID) -> df2
  log10(df2 + PS) %>% as_tibble() %>% mutate(ID = df$ID, .before=1) %>% return()
}
Pseudocount = function(df_tara){
  df_tara %>% select(-ID) %>% as.matrix() -> PS
  min(PS[PS!=0])/2 %>% return()
}
CLR_transformation = function( df_tara ){
  df_2 = select(df_tara, -ID)
  df_2 = df_2 + Pseudocount(df_tara)
  matrix_data = as.matrix(df_2)
  
  clr_matrix <- apply(matrix_data, 1, function(row){
    exp(mean(log(row)) ) -> geometric_mean
    log(row / geometric_mean)
  } ) %>% t()
  clr_matrix %>% as_tibble() %>% mutate(ID =df_tara$ID  , .before=1) %>% return()
  
}
Filter_prevalence = function(df_tara, Prevalence = 0.4 ){
  df_tara %>% select(-ID) %>% apply(2, function(x){ mean(x!=0) }  ) -> Prevalences
  Prevalences[Prevalences>=Prevalence] %>% names() -> KEEP
  df_tara %>% select(c("ID",KEEP) ) %>% return()
}


# Load metadata 
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_infants<-metadata[metadata$Type=="infant", ]
names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))


# Loading & processing cazyme data 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
caz_tables = Prepare_cazymes(File = "NEXT_cayman_infant.txt", Prevalence_min = 0.3, Transform = "NULL" )
caz_filtered = as.data.frame (caz_tables[["caz_filtered"]])
row.names(caz_filtered)<-caz_filtered$ID
caz_filtered$ID<-NULL



infant_taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
infant_taxa<-infant_taxa[match(row.names(caz_filtered),row.names(infant_taxa)),]
row.names (metadata_infants)<-metadata_infants$NG_ID
infant_taxa_metadata<-merge(metadata_infants, infant_taxa, by="row.names")
row.names (infant_taxa_metadata)<-infant_taxa_metadata$Row.names
infant_taxa_metadata$Row.names=NULL

infant_taxa_metadata<-infant_taxa_metadata[row.names(infant_taxa_metadata)%in% rownames(caz_filtered),] 
caz_filtered<-caz_filtered[row.names(caz_filtered)%in% rownames(infant_taxa_metadata),] 

my_pseudocount_normal=min(caz_filtered[caz_filtered!=0])/2# 
distance=vegdist(caz_filtered, method = "aitchison", pseudocount=my_pseudocount_normal) 

# Define metadata columns to exclude from testing (i.e., not bacterial taxa)
not_to_test <- c("NG_ID", "NEXT_ID", "FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                 "next_id_infant", "next_id_mother", "Modified_NEXT_ID_without_preg_number",
                 "days_from_first_collection", "human_reads_FQ_1", "human_reads_FQ_2",
                 "clean_reads_FQ_1", "clean_reads_FQ_2", "reads_lost_QC", "dna_conc",
                 "isolation_method", "NG_ID_short", "exact_age_days_at_collection",
                 "exact_age_months_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                 "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                 "contaminant_2_Phyllobacterium_myrsinacearum",
                 "metaphlan4_unclassified_with_contaminants", "shannon", "BATCH_NUMBER",
                 "next_id_partner", "sibling_number", "timepoint", "infant_relations",
                 "sequence_control", "isolation_control", "exact_age_years_at_collection")

# Prepare empty results dataframe
ResultsAdonis <- data.frame()

# Loop over each unique timepoint
for (tp in unique(infant_taxa_metadata$Timepoint_categorical)) {
  message("\n🔄 Processing timepoint: ", tp)
  
  # Subset metadata and distance matrix for the current timepoint
  sub_meta <- infant_taxa_metadata %>%
    filter(Timepoint_categorical == tp)
  
  samples <- rownames(sub_meta)
  sub_dist <- as.dist(as.matrix(distance)[samples, samples])
  
  # Select bacterial taxa columns only
  taxon_cols <- sub_meta %>%
    select(where(is.numeric)) %>%
    select(-any_of(not_to_test)) %>%
    colnames()
  
  message("Found ", length(taxon_cols), " taxa to test at this timepoint.")
  
  # Loop over each taxon
  for (taxon in taxon_cols) {
    message("  ➤ Testing taxon: ", taxon)
    
    # Drop NA values for current taxon
    sub_data <- sub_meta %>%
      select(all_of(taxon)) %>%
      drop_na()
    
    valid_samples <- rownames(sub_data)
    
    # Skip if fewer than 4 valid samples
    if (length(valid_samples) < 4) {
      message("Skipping (fewer than 4 samples with data).")
      next
    }
    
    # Subset distance matrix
    sub_dist2 <- as.dist(as.matrix(sub_dist)[valid_samples, valid_samples])
    
    # Build formula and run adonis2
    formula_str <- as.formula(paste("sub_dist2 ~", taxon))
    
    ad_result <- adonis2(formula_str, data = sub_data, permutations = 10000)
    
    # Store results
    ResultsAdonis <- bind_rows(ResultsAdonis, tibble(
      Timepoint = tp,
      Taxon = taxon,
      R2 = ad_result$R2[1],
      p_value = ad_result$`Pr(>F)`[1]
    ))
  }
}

print(ResultsAdonis)

write.table(ResultsAdonis,
            file = "adonis_cazyme_taxa_by_timepoint_results.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


results_adonis_infants <- read.delim("adonis_cazyme_taxa_by_timepoint_results.txt")
significant <- results_adonis_infants[results_adonis_infants$p_value < .0001, ]

# Clean taxon names
significant$Taxon <- as.character(significant$Taxon)
significant$Taxon <- gsub("\\.t__", " ", significant$Taxon)
significant$Taxon <- gsub("_", " ", significant$Taxon)

# Clean taxon order names
taxon_order <- names(infant_taxa)
clean_taxon_order <- gsub("\\.t__", " ", taxon_order)
clean_taxon_order <- gsub("_", " ", clean_taxon_order)

# Set factor levels
significant$Taxon <- factor(significant$Taxon, levels = clean_taxon_order, ordered = TRUE)

# Timepoint and R2 filtering
significant$Timepoint <- factor(significant$Timepoint,
                                levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"),
                                ordered = TRUE)

significant <- significant %>%
  mutate(R2 = round(R2, 3)) %>%
  filter(R2 > 0.06)

# Segment data for lines
segments <- significant %>%
  mutate(Timepoint_num = as.numeric(Timepoint)) %>%
  group_by(Taxon) %>%
  summarise(Timepoint2 = min(Timepoint_num),
            xend = max(Timepoint_num)) %>%
  mutate(Timepoint = levels(significant$Timepoint)[Timepoint2],
         xend = levels(significant$Timepoint)[xend])
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/submission_Nature_2025/supp_figures/Parts_of_figures")
pdf("pcos_bacteria_cazyme_profile_infants.pdf", width = 8, height = 6)
ggplot() +
  geom_segment(data = segments, aes(x = Timepoint, y = Taxon,
                                    xend = xend, yend = Taxon),
               color = "black", linewidth = 0.2) +
  geom_point(data = significant,
             aes(x = Timepoint, y = Taxon,
                 size = R2, fill = Timepoint),
             shape = 21, color = "black") +
  xlab("Timepoint") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_manual(values = c("#B5DD88", "#41C0B4", "#4397BB",
                               "#ECA4C9", "#cb4563", "#A42097", "#390962"))

dev.off()


##############################
# Mother
##############################


# Function to prepare cazyme data 
Prepare_cazymes = function(File = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/NEXT_cayman_mothers.txt", Prevalence_min = 0.4, Transform = "log" ){
  caz = read_tsv(File)
  caz %>% filter(! feature %in% c("total_reads", "filtered_reads")) -> caz
  #Get unique counts + ambigous counts normalized
  caz %>% select(ID, feature, combined_rpkm ) %>% spread(key = feature, value = combined_rpkm) -> caz
  #Change ID
  caz$ID = caz$ID 
  #Get summary
  summary_as_df <- function(x) {
    result <- summary(x)
    data.frame(stat = names(result), value = as.vector(result))
  }
  summary_df <- caz %>%
    select(-ID) %>%
    purrr::map(summary_as_df) %>%
    bind_rows(.id = "Cazyme") %>% spread(stat, value) %>% as_tibble()
  summary_df = summary_df %>% mutate( `NA's` = ifelse( is.na(`NA's`), 0,  `NA's` ), N = dim(caz)[1], Present= N - `NA's`)
  summary_df = summary_df %>% mutate(Present_perc = Present/N ) %>% mutate(Keep = ifelse(Present_perc>=Prevalence_min, T, F ) )
  #Make NA into 0
  caz %>% mutate_all(~ replace(., is.na(.), 0)) -> caz
  #Remove Cazymes with high NA prop
  caz %>% select(c("ID", summary_df$Cazyme[summary_df$Keep==T] )) -> caz_anal
  if (Transform == "log"){
    Log_transform(caz_anal)  ->caz_anal_tranf
  } else if (Transform == "clr"){
    CLR_transformation(caz_anal)  -> caz_anal_tranf
  } else { caz_anal_tranf = NULL}
  
  return(list("caz" = caz, "caz_filtered" = caz_anal, "caz_tranf"=caz_anal_tranf  ,"Summary_caz" = summary_df ) )
}


Log_transform = function(df){
  colnames(df)[1] = "ID"
  PS = Pseudocount(df)
  df %>% select(-ID) -> df2
  log10(df2 + PS) %>% as_tibble() %>% mutate(ID = df$ID, .before=1) %>% return()
}
Pseudocount = function(df_tara){
  df_tara %>% select(-ID) %>% as.matrix() -> PS
  min(PS[PS!=0])/2 %>% return()
}
CLR_transformation = function( df_tara ){
  df_2 = select(df_tara, -ID)
  df_2 = df_2 + Pseudocount(df_tara)
  matrix_data = as.matrix(df_2)
  
  clr_matrix <- apply(matrix_data, 1, function(row){
    exp(mean(log(row)) ) -> geometric_mean
    log(row / geometric_mean)
  } ) %>% t()
  clr_matrix %>% as_tibble() %>% mutate(ID =df_tara$ID  , .before=1) %>% return()
  
}
Filter_prevalence = function(df_tara, Prevalence = 0.4 ){
  df_tara %>% select(-ID) %>% apply(2, function(x){ mean(x!=0) }  ) -> Prevalences
  Prevalences[Prevalences>=Prevalence] %>% names() -> KEEP
  df_tara %>% select(c("ID",KEEP) ) %>% return()
}


# Load metadata 
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
metadata_mothers<-metadata[metadata$Type=="mother", ]
names (metadata_mothers)[2]<-"next_id_mother"
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
metadata_mothers$Timepoint_categorical <- factor(metadata_mothers$Timepoint_categorical, 
                                                 levels = c("P12", "P28", "B", "M3"))
# Loading & processing cazyme data 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
caz_tables = Prepare_cazymes(File = "NEXT_cayman_mothers.txt", Prevalence_min = 0.3, Transform = "NULL" )
caz_filtered = as.data.frame (caz_tables[["caz_filtered"]])
row.names(caz_filtered)<-caz_filtered$ID
caz_filtered$ID<-NULL



mother_taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/NEXT_metaphlan_4_CLR_transformed_fil_SGB_mothers_03_08_2023.txt")
mother_taxa<-mother_taxa[match(row.names(caz_filtered),row.names(mother_taxa)),]
row.names (metadata_mothers)<-metadata_mothers$NG_ID
mother_taxa_metadata<-merge(metadata_mothers, mother_taxa, by="row.names")
row.names (mother_taxa_metadata)<-mother_taxa_metadata$Row.names
mother_taxa_metadata$Row.names=NULL

mother_taxa_metadata<-mother_taxa_metadata[row.names(mother_taxa_metadata)%in% rownames(caz_filtered),] 
caz_filtered<-caz_filtered[row.names(caz_filtered)%in% rownames(mother_taxa_metadata),] 

my_pseudocount_normal=min(caz_filtered[caz_filtered!=0])/2# 
distance=vegdist(caz_filtered, method = "aitchison", pseudocount=my_pseudocount_normal) 

# Define metadata columns to exclude from testing (i.e., not bacterial taxa)
not_to_test <- c("NG_ID", "NEXT_ID", "FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                 "next_id_mother", "next_id_mother", "Modified_NEXT_ID_without_preg_number",
                 "days_from_first_collection", "human_reads_FQ_1", "human_reads_FQ_2",
                 "clean_reads_FQ_1", "clean_reads_FQ_2", "reads_lost_QC", "dna_conc",
                 "isolation_method", "NG_ID_short", "exact_age_days_at_collection",
                 "exact_age_months_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                 "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                 "contaminant_2_Phyllobacterium_myrsinacearum",
                 "metaphlan4_unclassified_with_contaminants", "shannon", "BATCH_NUMBER",
                 "next_id_partner", "sibling_number", "timepoint", "mother_relations",
                 "sequence_control", "isolation_control", "exact_age_years_at_collection")

# Prepare empty results dataframe
ResultsAdonis <- data.frame()

# Loop over each unique timepoint
for (tp in unique(mother_taxa_metadata$Timepoint_categorical)) {
  message("\n🔄 Processing timepoint: ", tp)
  
  # Subset metadata and distance matrix for the current timepoint
  sub_meta <- mother_taxa_metadata %>%
    filter(Timepoint_categorical == tp)
  
  samples <- rownames(sub_meta)
  sub_dist <- as.dist(as.matrix(distance)[samples, samples])
  
  # Select bacterial taxa columns only
  taxon_cols <- sub_meta %>%
    select(where(is.numeric)) %>%
    select(-any_of(not_to_test)) %>%
    colnames()
  
  message("Found ", length(taxon_cols), " taxa to test at this timepoint.")
  
  # Loop over each taxon
  for (taxon in taxon_cols) {
    message("  ➤ Testing taxon: ", taxon)
    
    # Drop NA values for current taxon
    sub_data <- sub_meta %>%
      select(all_of(taxon)) %>%
      drop_na()
    
    valid_samples <- rownames(sub_data)
    
    # Skip if fewer than 4 valid samples
    if (length(valid_samples) < 4) {
      message("Skipping (fewer than 4 samples with data).")
      next
    }
    
    # Subset distance matrix
    sub_dist2 <- as.dist(as.matrix(sub_dist)[valid_samples, valid_samples])
    
    # Build formula and run adonis2
    formula_str <- as.formula(paste("sub_dist2 ~", taxon))
    
    ad_result <- adonis2(formula_str, data = sub_data, permutations = 10000)
    
    # Store results
    ResultsAdonis <- bind_rows(ResultsAdonis, tibble(
      Timepoint = tp,
      Taxon = taxon,
      R2 = ad_result$R2[1],
      p_value = ad_result$`Pr(>F)`[1]
    ))
  }
}

print(ResultsAdonis)

write.table(ResultsAdonis,
            file = "adonis_cazyme_taxa_by_timepoint_results_mothers.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

results_adonis_mothers <- read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/adonis_cazyme_taxa_by_timepoint_results_mothers.txt")
significant <- results_adonis_mothers[results_adonis_mothers$p_value < .0001, ]

# Clean Taxon names
significant$Taxon <- as.character(significant$Taxon)
significant$Taxon <- gsub("\\.t__", " ", significant$Taxon)
significant$Taxon <- gsub("_", " ", significant$Taxon)

# Clean taxon order
taxon_order <- names(mother_taxa)
clean_taxon_order <- gsub("\\.t__", " ", taxon_order)
clean_taxon_order <- gsub("_", " ", clean_taxon_order)

# Set factor levels
significant$Taxon <- factor(significant$Taxon, levels = clean_taxon_order, ordered = TRUE)

# Timepoint and R2
significant$Timepoint <- factor(significant$Timepoint,
                                levels = c("P12", "P28", "B", "M3"),
                                ordered = TRUE)

significant <- significant %>%
  mutate(R2 = round(R2, 3)) %>%
  filter(R2 > 0.03)

# Segment data
segments <- significant %>%
  mutate(Timepoint_num = as.numeric(Timepoint)) %>%
  group_by(Taxon) %>%
  summarise(Timepoint2 = min(Timepoint_num),
            xend = max(Timepoint_num)) %>%
  mutate(Timepoint = levels(significant$Timepoint)[Timepoint2],
         xend = levels(significant$Timepoint)[xend])

# Save plot to PDF
pdf("pcos_bacteria_cazyme_profile_mothers.pdf", width = 8, height = 6)

ggplot() +
  geom_segment(data = segments, aes(x = Timepoint, y = Taxon,
                                    xend = xend, yend = Taxon),
               color = "black", linewidth = 0.2) +
  geom_point(data = significant,
             aes(x = Timepoint, y = Taxon,
                 size = R2, fill = Timepoint),
             shape = 21, color = "black") +
  xlab("Timepoint") +
  ylab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_size_continuous(range = c(1, 4)) +
  scale_fill_manual(values = c("#f90404", "#f78310", "#fbd123", "#eca4c9"))

dev.off()
