######## MATERNAL HMO & INFANT MICROBIOME ##############
### AUTHOR:  TRISHLA SINHA
### ORIGINAL SCRIPT: 15th November, 2023
### LAST UPDATE: 1st of August 2024

library(tidyverse)


##### FUNCTIONS ##########

# To remove phenotypes with too many NA's and too little variance 
drop_phenotypes <- function(phenotype_table, non_NA_values, minimum_variance) {
  cleaned_phenotypes <- phenotype_table #create a variable with phenotypic data
  discarded_phenotypes <- NULL
  for (i in 1:ncol(cleaned_phenotypes)) { 
    # check if the number of non-NA values is below the threshold
    if (length(which(!is.na(cleaned_phenotypes[,i]))) < non_NA_values) { 
      print(paste(colnames(cleaned_phenotypes)[i], "has less than n non-NA values"))
      discarded_phenotypes <- c(discarded_phenotypes, colnames(cleaned_phenotypes[i]))
      # check if the variance is below the threshold
    } else if (as.vector(sort(table(cleaned_phenotypes[,i]), decreasing = T) [1]) > length(which(!is.na(cleaned_phenotypes[,i]))) * (1-minimum_variance*0.01)) {
      print(paste(colnames(cleaned_phenotypes)[i], "has less than n% variance"))
      discarded_phenotypes <- c(discarded_phenotypes, colnames(cleaned_phenotypes[i]))
      # the phenotypes that do not fullfill the previous statements should be included
    } else {
      print(paste(colnames(cleaned_phenotypes)[i], ": accepted"))
    }
  } 
  cleaned_phenotypes[, discarded_phenotypes] <- NULL #drop non-selected phenotypes
  return(cleaned_phenotypes)
}

# To inverse rank tansform numeric data 
run_invrank_dataFrame = function(data){
  invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))} 
  data_types = unlist(lapply(data,class))
  newdata = data
  for(i in 1:ncol(data)){
    if(data_types[i]=="numeric"|data_types[i]=="integer") {
      newdata[,i] = invrank(data[,i])
    }
  }
  newdata
}


gam_function <- function(taxa, phenotypes) {
  metadata_columns <- c("NG_ID", "NEXT_ID","Modified_NEXT_ID_without_preg_number", "days_from_first_collection","FAMILY", "raw_reads_FQ_1", "raw_reads_FQ_2",
                        "human_reads_FQ_1", "human_reads_FQ_2", "clean_reads_FQ_1", "clean_reads_FQ_2",
                        "reads_lost_QC", "dna_conc", "isolation_method", "NG_ID_short",
                        "exact_age_days_at_collection", "exact_age_months_at_collection",
                        "exact_age_years_at_collection", "Timepoint_categorical", "SAMPLE_ID",
                        "metaphlan4_unclassified", "contaminant_1_Sphingomonas_sp_FARSPH",
                        "contaminant_2_Phyllobacterium_myrsinacearum", "metaphlan4_unclassified_with_contaminants",
                        "shannon", "BATCH_NUMBER", "next_id_mother", "next_id_partner", "sibling_number", "timepoint")
  
  phenotypes2correlate <- phenotypes[, !(colnames(phenotypes) %in% metadata_columns)] 
  phenotypes2correlate[sapply(phenotypes2correlate, is.character)] <- lapply(phenotypes2correlate[sapply(phenotypes2correlate, is.character)], as.factor)
  
  covariates <- phenotypes[, c("clean_reads_FQ_1","dna_conc", "BATCH_NUMBER")] 
  
  phenotypes2correlate <- run_invrank_dataFrame(phenotypes2correlate) 
  covariates <- run_invrank_dataFrame(covariates)
  
  result = foreach(i = 1:ncol(taxa), .combine = rbind) %:% 
    foreach(j = 1:ncol(phenotypes2correlate), .combine = rbind) %do% {
      taxa_name = colnames(taxa)[i]
      phenotype_name = colnames(phenotypes2correlate)[j]
      print(paste("Processing taxa:", taxa_name, "- Phenotype:", phenotype_name))
      data.fit <- data.frame(bac = taxa[,i],
                             trait = phenotypes2correlate[,j],
                             Time = phenotypes$exact_age_months_at_collection,
                             RD = covariates$clean_reads_FQ_1,
                             Batch = covariates$BATCH_NUMBER,
                             DNAcon = covariates$dna_conc,
                             NEXT_ID = phenotypes$NEXT_ID)
      data.fit$NEXT_ID <- factor(data.fit$NEXT_ID)
      data.fit <- data.fit[complete.cases(data.fit),]
      
      mod_gam1 <- tryCatch({
        bam(bac ~ trait + RD + DNAcon + Batch + s(Time) + s(NEXT_ID, bs = 're'), data = data.fit, discrete = TRUE)
      }, warning = function(w) {
        return(NULL)  
      }, error = function(e) {
        return(NULL)  
      })
      
      if (is.null(mod_gam1)) {
        statistics <- NULL
        converged <- FALSE
        beta <- NA
        SE <- NA
        P <- NA
        N <- nrow(data.fit)
      } else {
        summary.results <- summary(mod_gam1)
        statistics <- summary.results$p.table[2:(which(rownames(summary.results$p.table) == 'RD') - 1),, drop = FALSE]
        converged <- TRUE
        beta <- statistics[,1]
        SE <- statistics[,2]
        P <- statistics[,4]
        N <- nrow(data.fit)
      }
      
      output <- data.frame(
        species = colnames(taxa)[i],
        variable = colnames(phenotypes2correlate)[j],
        effect.level = if(!is.null(statistics)) sub("^trait", "", rownames(statistics)) else NA,
        beta = beta,
        SE = SE,
        P = P,
        N = N,
        converged = converged
      )
      
      output
    }
  result
}


# Load metadata and phenotypes
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/hmo")
#hmo<-read.delim("hmo_RA_20_07_2024.txt")
#metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
#linkage<-read.delim("~/Desktop/LLNEXT/Analysis/linkage/20230414_Linkingfile_fam_longID_formatted_TS_SB.txt")
#names (hmo)[2]<-"next_id_mother"
#names (hmo)[1]<-"sample_id_mother"

# Modify HMO data to get infant SAMPLE_ID's 
#hmo_linkage<-left_join(hmo, linkage)
#hmo_linkage$SAMPLE_ID<-paste0(hmo_linkage$next_id_infant,"_", hmo_linkage$time_point )

# Merging files and selecting infant relevant data
#metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                  # as.factor)
#metadata_infants<-metadata[metadata$Type=="infant", ]
#names (metadata_infants)[2]<-"next_id_infant"
#metadata_infants$BATCH_NUMBER<-as.factor(metadata_infants$BATCH_NUMBER)

# For this specific data remove duplicates for categorical timepoints 
#metadata_infants <- metadata_infants %>%
  #distinct(SAMPLE_ID, .keep_all = TRUE) 

# Merging HMO with microbiome metadata 
#hmo_met<-left_join(metadata_infants, hmo_linkage)
#write.table(hmo_met, "HMOs_RA_matched_infant_microbiome.txt", sep = "\t", row.names = F)


# Remove phenotypes with too many NA's 
hmo_met <-read.delim("HMOs_RA_matched_infant_microbiome_NA_outliers_removed.txt")
names (hmo_met)[2]<-"NEXT_ID"
hmo_met$NEXT_ID=as.factor(hmo_met$NEXT_ID)
row.names(hmo_met)<-hmo_met$NG_ID

hmo_met$twin_pair=NULL
hmo_met$infant_relations=NULL
hmo_met$Type=NULL
hmo_met$Sequenced=NULL
hmo_met$sequence_control=NULL
hmo_met$isolation_control=NULL
hmo_met$metaphlan4_unclassified_high_contaminants_factor_75=NULL
hmo_met$NEXT_ID_mother_simple=NULL
hmo_met$pregnancy=NULL
hmo_met$sample_id_mother=NULL
hmo_met$time_point=NULL
hmo_met$mother_milk_HMO_UPLC_ID=NULL
hmo_met <- hmo_met[!(hmo_met$Timepoint_categorical %in% c("M6", "M9", "M12")), ]

hmo_met$mother_milk_HMO_milk_group <- factor(hmo_met$mother_milk_HMO_milk_group, levels = c("Le+Se+", "Le+Se-", "Le-Se+", "Le-Se-"))

# Remove mixed feeding and only do breastfeeding 
dynamic_phenotypes<-read.delim("~/Desktop/LLNEXT/Analysis/phenotypes/masterfile_longitudinal_2023_09_29.txt")
feeding <-dynamic_phenotypes[, c("SAMPLE_ID", "infant_ffq_feeding_mode_simple")]

hmo_met_feed<-left_join(hmo_met, feeding)
hmo_met_feed <- hmo_met_feed[(hmo_met_feed$infant_ffq_feeding_mode_simple %in% c("excl_BF")), ]
hmo_met_feed <- hmo_met_feed[!is.na(hmo_met_feed$mother_milk_HMO_milk_group), ]
row.names(hmo_met_feed)<-hmo_met_feed$NG_ID
num_unique_infants <- length(unique(hmo_met_feed$NEXT_ID))
# 277

#Taxa
taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")

# Analysis

# Select outcome 
alpha<-hmo_met_feed %>% select(shannon)

# Select merged phenotype & metadata file 
result_alpha <- gam_function(alpha, hmo_met_feed)
result_alpha$FDR<-p.adjust (result_alpha$P, method = "fdr")
result_alpha$trait_group <- sapply(strsplit(as.character(result_alpha$variable), "_"), function(x) paste(x[1:2], collapse = "_"))

table (hmo_met_feed$mother_milk_HMO_milk_group, hmo_met_feed$Timepoint_categorical)
table (hmo_met_feed$mother_milk_HMO_milk_group, hmo_met_feed$infant_ffq_feeding_mode_simple, hmo_met_feed$Timepoint_categorical)


ggplot(hmo_met_feed, aes(mother_milk_HMO_milk_group, y = shannon, fill = mother_milk_HMO_milk_group, color = mother_milk_HMO_milk_group)) +
  scale_fill_manual(values = wes_palette("BottleRocket2", n = 4)) +
  scale_color_manual(values = wes_palette("BottleRocket2", n = 4)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() +
  labs(x = "", y = "Shannon Diversity Index", fill = "Mother Milk HMO Milk Group", color = "Mother Milk HMO Milk Group") +
  theme(
    plot.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1)
  )

taxa_1 <- taxa[match(rownames(hmo_met_feed), rownames(taxa)),]

result_taxa <- gam_function(taxa_1, hmo_met_feed)
result_taxa$FDR<-p.adjust (result_taxa$P, method = "fdr")
result_taxa$trait_group <- sapply(strsplit(as.character(result_taxa$variable), "_"), function(x) paste(x[1:2], collapse = "_"))

pathways<-read.delim("~/Desktop/LLNEXT/Analysis/pathways/NEXT_humann_pathways_INFANTS_CLR_transformed_fil_30_percent_10_07_2024.txt")

pathways_1<- pathways[match(rownames(hmo_met_feed), rownames(pathways)),]
pathways_1<-na.omit(pathways_1)

result_pathways <- gam_function(pathways_1, hmo_met_feed)
result_pathways$FDR<-p.adjust (result_pathways$P, method = "fdr")
result_pathways$trait_group <- sapply(strsplit(as.character(result_pathways$variable), "_"), function(x) paste(x[1:2], collapse = "_"))


setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/hmo")

write.table(result_alpha, "hmo_microbiome/HMO_infant_microbiome_alpha_changed_levels.txt", sep = "\t", row.names = F)
write.table(result_taxa, "hmo_microbiome/HMO_infant_microbiome_taxa_changed_levels.txt", sep = "\t", row.names = F)
write.table(result_pathways, "hmo_microbiome/HMO_infant_microbiome_pathways.txt", sep = "\t", row.names = F)

all<-merge(hmo_met_feed, taxa, by="row.names")
#all2<-merge(hmo_met, alpha, by="row.names")

ggplot(all, aes(mother_milk_HMO_milk_group, y = Clostridium_perfringens.t__SGB6191, fill = mother_milk_HMO_milk_group, color = mother_milk_HMO_milk_group)) +
  scale_fill_manual(values = wes_palette("BottleRocket2", n = 4)) +
  scale_color_manual(values = wes_palette("BottleRocket2", n = 4)) +
  geom_boxplot(alpha = 0.4, outlier.colour = NA) +
  geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  ggtitle("") +
  theme_bw() +
  labs(x = "", y = "Clostridium_perfringens.t__SGB6191", fill = "Mother Milk HMO Milk Group", color = "Mother Milk HMO Milk Group") +
  theme(
    plot.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1)
  )


hmo_results <-read.delim("hmo_microbiome/HMO_infant_microbiome_taxa.txt")
hmo_results <- na.omit(hmo_results)
hmo_results_nom_sig<-hmo_results[hmo_results$P<0.05,]
hmo_results_nom_sig$hmo<-paste0(hmo_results_nom_sig$variable,hmo_results_nom_sig$effect.level)
hmo_results_nom_sig$t.value <- hmo_results_nom_sig$beta/ hmo_results_nom_sig$SE

wide_data <- dcast(hmo_results_nom_sig, species ~ hmo, value.var = "t.value")
row.names(wide_data)<-wide_data$species
wide_data$species=NULL
wide_data[is.na(wide_data)] <- 0

# Convert the data frame to a matrix

heatmap_data <- as.matrix(t(wide_data))

# Plot clustered heatmap

pdf ("supplementary_figure_1_B.pdf", 
     width = 12, height = 5)
hmo<-pheatmap(heatmap_data,
                 cluster_rows = F, 
                 clustering_distance_cols = "euclidean", 
                 clustering_method = "complete",
                 color = colorRampPalette(c("blue", "white", "red"))(50),
                 display_numbers = F,
                 main = "Heatmap of t values of HMOs with infant microbiome", 
                 angle_col = 90,
                 fontsize_col = 8)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/supplementary")
infant


all<-merge(hmo_met_feed, pathways, by="row.names")


hmo_results <- na.omit(result_pathways)
hmo_results_nom_sig<-hmo_results[hmo_results$P<0.05,]
hmo_results_nom_sig$hmo<-paste0(hmo_results_nom_sig$variable,hmo_results_nom_sig$effect.level)
hmo_results_nom_sig$t.value <- hmo_results_nom_sig$beta/ hmo_results_nom_sig$SE

wide_data <- dcast(hmo_results_nom_sig, species ~ hmo, value.var = "t.value")
row.names(wide_data)<-wide_data$species
wide_data$species=NULL
wide_data[is.na(wide_data)] <- 0

# Convert the data frame to a matrix

heatmap_data <- as.matrix(t(wide_data))

# Plot clustered heatmap

pdf ("supplementary_figure_1_B.pdf", 
     width = 12, height = 5)
hmo<-pheatmap(heatmap_data,
              cluster_rows = F, 
              clustering_distance_cols = "euclidean", 
              clustering_method = "complete",
              color = colorRampPalette(c("blue", "white", "red"))(50),
              display_numbers = F,
              main = "Heatmap of t values of HMOs with infant pathways", 
              angle_col = 90,
              fontsize_col = 4,
              fontsize_row = 4)
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/figures/supplementary")
infant

