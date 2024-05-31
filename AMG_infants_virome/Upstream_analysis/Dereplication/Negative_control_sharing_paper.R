## Code description
#######################################################################################################################################
## Script for the negative control sharing paper
## 
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
library(readr)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
RPKM <- read.delim('RPKM_filtered_v2.txt')

extended_tof <- read.delim('MERGED_Extended_TOF')
row.names(extended_tof) <- extended_tof$New_CID

meta_all_with_qc_curated <- as.data.frame(read_tsv("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\metadata\\metadata_with_qc_v3.tsv"))
#######################################################################################################################################

## Filtering out the GFMice, formula virome and positive control samples from the metadata & RPKM filtered
## Filtering the samples by the presence of viral genes
#######################################################################################################################################
meta_working <- meta_all_with_qc_curated[!(meta_all_with_qc_curated$Type %in%
                                             c("Pos_ctlr", "GFMice", "Formula_virome", "Induced_phages", NA)), 
                                         !(names(meta_all_with_qc_curated) %in% 
                                             c("infant_sex", "infant_gestational_age", "infant_birthweight", "mother_age_years", 
                                               "mother_bmi_category", "infant_place_delivery", "Cows_milk", "infant_zygosity", "Race"))]

RPKM_filtered <- RPKM[, names(RPKM) %in% meta_working$Sample_name]

virus_selection <- extended_tof$New_CID[extended_tof$viral_genes > 0 | (extended_tof$viral_genes == 0 & extended_tof$host_genes == 0)]

RPKM_filtered <- RPKM_filtered[row.names(RPKM_filtered) %in% virus_selection, ]
RPKM_filtered <- RPKM_filtered[rowSums(RPKM_filtered) > 0, colSums(RPKM_filtered) > 0]
#######################################################################################################################################

## Analysis for the vOTUs found in negative controls: Bray-curtis heat-map
#######################################################################################################################################
negative_controls <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" &
                                                meta_working$Sample_name %in% colnames(RPKM_filtered)]

RPKM_filtered_nc <- RPKM_filtered[, colnames(RPKM_filtered) %in% negative_controls]
RPKM_filtered_nc <- RPKM_filtered_nc[rowSums(RPKM_filtered_nc) > 0, colSums(RPKM_filtered_nc) > 0]

# preapring phenotypes
meta_nc_nmds <- meta_working[meta_working$Sample_name %in% colnames(RPKM_filtered_nc), c("Sample_name", "cohort", "clean_reads_comb", "total_viruses_discovered")]
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- merge(meta_nc_nmds, as.data.frame(t(as.data.frame(lapply(RPKM_filtered_nc, diversity, index = "shannon")))),
                      by="row.names")
row.names(meta_nc_nmds) <- meta_nc_nmds$Sample_name
meta_nc_nmds <- meta_nc_nmds[c("cohort", "clean_reads_comb", "total_viruses_discovered", "V1")]
names(meta_nc_nmds) <- c("cohort", "clean_reads_comb", "total_viruses_discovered", "shannon_div")

## Bray-curtis heat map analysis
nc_bray_dist_matrix <- as.matrix(vegdist(t(RPKM_filtered_nc), method="bray"))

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pheatmap(nc_bray_dist_matrix, 
         annotation = annotation, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         display_numbers = TRUE, 
         color = colorRampPalette(c("#D2042D", "#FAA0A0"))(50),
         main = "Negative controls Bray-Curtis Distance")
#######################################################################################################################################

## Analysis for the vOTUs found in negative controls: NMDS
#######################################################################################################################################
ord <- metaMDS(t(RPKM_filtered_nc), distance = "bray", k=2, parallel=8)

en = envfit(ord, meta_nc_nmds, permutations = 999, na.rm = TRUE)
en$factors #r2 0.0031; p-value = 0.988
en$vectors


centroids <- as.data.frame(scores(en, "factors"))
centroids$cohort <- c(gsub('cohort', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, meta_nc_nmds, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=cohort)) + 
  geom_point(size = 2, alpha=0.8) + 
  # xlim(-15, 5) +  # Several samples from liang and maqsood were severe outliers, thus they were excluded
  # ylim(-15, 10) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = cohort, fill=cohort), linetype = 2) +
  geom_point(data=centroids, aes(fill=cohort),shape=23, size=4, color='black', ) + 
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") 
#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts
#######################################################################################################################################
negative_controls_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "garmaeva" & 
                                                            meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "liang" &
                                                         meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "maqsood" &
                                                           meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "shah" &
                                                        meta_working$Sample_name %in% colnames(RPKM_filtered)]


negative_controls_wo_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort != "garmaeva" & 
                                                            meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_wo_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort != "liang" &
                                                         meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_wo_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort != "maqsood" &
                                                           meta_working$Sample_name %in% colnames(RPKM_filtered)]

negative_controls_wo_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort != "shah" &
                                                        meta_working$Sample_name %in% colnames(RPKM_filtered)]

RPKM_filtered_w_dummy <- RPKM_filtered

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy$dummy_wo_garmaeva <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_wo_garmaeva])
RPKM_filtered_w_dummy$dummy_wo_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_wo_liang])
RPKM_filtered_w_dummy$dummy_wo_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_wo_maqsood])
RPKM_filtered_w_dummy$dummy_wo_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_wo_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

garmaeva_samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "garmaeva"]
maqsood_samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "maqsood"]
liang_samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "liang"]
shah_samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "shah"]

temp_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% garmaeva_samples] + RPKM_filtered_w_dummy$dummy_wo_garmaeva
temp_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% maqsood_samples] + RPKM_filtered_w_dummy$dummy_wo_maqsood
temp_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% liang_samples] + RPKM_filtered_w_dummy$dummy_wo_liang
temp_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% shah_samples] + RPKM_filtered_w_dummy$dummy_wo_shah

temp_garmaeva_maqsood <- merge(temp_garmaeva, temp_maqsood, by="row.names", all = T)
temp_liang_shah <- merge(temp_liang, temp_shah, by="row.names", all = T)
temp_all <- merge(temp_garmaeva_maqsood, temp_liang_shah, by="Row.names", all = T)
row.names(temp_all) <- temp_all$Row.names
temp_all$Row.names <- NULL
temp_all[temp_all < 2] <- 0

table_for_plot_a <- as.data.frame(colSums(temp_all) / 2)

temp_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% garmaeva_samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% maqsood_samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% liang_samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% shah_samples] + RPKM_filtered_w_dummy$dummy_shah

temp_garmaeva_maqsood <- merge(temp_garmaeva, temp_maqsood, by="row.names", all = T)
temp_liang_shah <- merge(temp_liang, temp_shah, by="row.names", all = T)
temp_all <- merge(temp_garmaeva_maqsood, temp_liang_shah, by="Row.names", all = T)
row.names(temp_all) <- temp_all$Row.names
temp_all$Row.names <- NULL
temp_all[temp_all < 2] <- 0

table_for_plot_b <- as.data.frame(colSums(temp_all) / 2)

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_c <- as.data.frame(colSums(RPKM_count))

table_for_plot <- merge(table_for_plot_a, table_for_plot_b, by='row.names')
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL
table_for_plot <- merge(table_for_plot, table_for_plot_c, by='row.names', all.x = T)


colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_other_cohorts", "Shared_w_nc_same_cohort", "Viral_richness")

table_for_plot$Perc_shared_w_other_cohorts <- table_for_plot$Shared_w_nc_other_cohorts / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_other_cohorts <- table_for_plot$Shared_w_nc_same_cohort / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_other_cohorts", "Shared_w_nc_same_cohort", "Viral_richness",
                                                   "Type", "Timepoint", "cohort"))

ggplot(table_for_plot_melt, aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=2) +
  facet_grid(cohort ~ Type) +
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        strip.text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))




## To do
#######################################################################################################################################
# add the savings for the plot
# add walters to the last analysis
# now I only have sharedness with the negative controls from other cohorts, but I need one more for sharedness with the same cohort
# color the dots in the violin plot by timepoint
#######################################################################################################################################

