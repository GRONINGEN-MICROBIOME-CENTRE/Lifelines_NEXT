## Code description
#######################################################################################################################################
## Script for the negative control sharing paper
## 
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
.libPaths("/home1/p309176/R_libraries")
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
setwd("/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/plots")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
RPKM <- read.delim("../RPKM_counts_VLP_NCP.txt")
dim(RPKM)  # 306275   1301

extended_tof <- read.delim('../../VIR_DB/virus_contigs/MERGED_Extended_TOF_NCP')
row.names(extended_tof) <- extended_tof$New_CID
dim(extended_tof)  # 971583     41

meta_all_with_qc_curated <- as.data.frame(read_tsv('../../metadata_with_qc_NCP.tsv'))
dim(meta_all_with_qc_curated)  # 1376   28
#######################################################################################################################################

##############################################################################################################################################################################################################################################################################


# THIS BLOCK IS DEDICATED TO EXPLORING THE RESULTS, WHEN THE RPKM TABLE IS NOT FILTERED


##############################################################################################################################################################################################################################################################################

#######################################################################################################################################
meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM), ]
#######################################################################################################################################

## Analysis for the vOTUs found in negative controls: Bray-curtis heat-map
#######################################################################################################################################
negative_controls <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" &
                                                meta_working$Sample_name %in% colnames(RPKM)]

RPKM_filtered_nc <- RPKM[, colnames(RPKM) %in% negative_controls]
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
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix

# heatmap.2(x=nc_bray_dist_matrix, dendrogram=F)

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()



#######################################################################################################################################

## Analysis for the vOTUs composition
#######################################################################################################################################

extended_tof <- extended_tof %>%
  mutate(checkv_quality_and_plasmid = ifelse(plasmid == "No", checkv_quality, "plasmid"))


RPKM_filtered_nc_count <- RPKM_filtered_nc
RPKM_filtered_nc_count[RPKM_filtered_nc_count > 0] <- 1


RPKM_filtered_nc_count <- merge(RPKM_filtered_nc_count, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
RPKM_filtered_nc_count$Row.names <- NULL

summarized_df <- RPKM_filtered_nc_count %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

long_df <- summarized_df %>%
  pivot_longer(cols = -checkv_quality_and_plasmid, 
               names_to = "sample", 
               values_to = "count")

total_counts <- long_df %>%
  group_by(sample) %>%
  summarise(total_count = sum(count))


# Create the plot and save it as a PDF
pdf('summary_per_nc_absolute.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df, aes(x = sample, y = count, fill = checkv_quality_and_plasmid)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#f78c6b", "#06d6a0", "#ffd166", "#118ab2", "#073b4c")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Virus Type") +
  geom_text(data = total_counts, aes(x = sample, y = total_count, label = total_count), 
            angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()

long_df <- long_df + 1

pdf('summary_per_nc_log.pdf', width = 12/2.54, height = 8/2.54)
ggplot(long_df, aes(x = sample, y = count, fill = checkv_quality_and_plasmid)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#ef476f", "#f78c6b", "#06d6a0", "#ffd166", "#118ab2", "#073b4c")) +
  scale_y_log10() +  # Set y-axis to log scale
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),  # Rotate x-ticks and decrease font size
    axis.text.y = element_text(size = 4),  # Decrease y-axis font size
    axis.title = element_text(size = 5),  # Decrease axis titles font size
    legend.title = element_text(size = 5),  # Decrease legend title font size
    legend.text = element_text(size = 4),  # Decrease legend text font size
    legend.key.size = unit(0.3, "cm")  # Decrease the size of the color boxes in the legend
  ) +
  labs(x = "Sample", y = "Count", fill = "Virus Type") +
  geom_text(data = total_counts, aes(x = sample, y = total_count, label = total_count), 
            angle = 45, size = 1, inherit.aes = FALSE)  # Position text above the bars on a single line
dev.off()

#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts
#######################################################################################################################################
negative_controls_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "garmaeva" & 
                                                            meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "liang" &
                                                         meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "maqsood" &
                                                           meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "shah" &
                                                        meta_working$Sample_name %in% colnames(RPKM)]

RPKM_filtered_w_dummy <- RPKM

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl"]

temp_shared_w_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_shared_w_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_shared_w_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shared_w_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_shah

temp_shared_w_garmaeva[temp_shared_w_garmaeva < 2] <- 0
temp_shared_w_maqsood[temp_shared_w_maqsood < 2] <- 0
temp_shared_w_liang[temp_shared_w_liang < 2] <- 0
temp_shared_w_shah[temp_shared_w_shah < 2] <- 0

table_for_plot_g <- as.data.frame(colSums(temp_shared_w_garmaeva) / 2)
table_for_plot_m <- as.data.frame(colSums(temp_shared_w_maqsood) / 2)
table_for_plot_l <- as.data.frame(colSums(temp_shared_w_liang) / 2)
table_for_plot_s <- as.data.frame(colSums(temp_shared_w_shah) / 2)

table_for_plot_g_m <- merge(table_for_plot_g, table_for_plot_m, by="row.names", all = T)
table_for_plot_l_s <- merge(table_for_plot_l, table_for_plot_s, by="row.names", all = T)
table_for_plot <- merge(table_for_plot_g_m, table_for_plot_l_s, by="Row.names", all = T)
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_add <- as.data.frame(colSums(RPKM_count))
table_for_plot <- merge(table_for_plot, table_for_plot_add, by='row.names', all.x = T)

colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                              "Shared_w_nc_shah", "Viral_richness")

table_for_plot$Perc_shared_w_garmaeva <- table_for_plot$Shared_w_nc_garmaeva / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_maqsood <- table_for_plot$Shared_w_nc_maqsood / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_liang <- table_for_plot$Shared_w_nc_liang / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_shah <- table_for_plot$Shared_w_nc_shah / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                                                   "Shared_w_nc_shah", "Viral_richness", "Type", "Timepoint", "cohort"))

table_for_plot_melt$variable <- factor(table_for_plot_melt$variable, levels = c("Perc_shared_w_garmaeva", "Perc_shared_w_liang", 
                                                                                "Perc_shared_w_maqsood", "Perc_shared_w_shah"))

## Plotting only sharedness for infants
pdf('shared_inf.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ ., scales = "free") +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

## Plotting only sharendess for mothers
pdf('shared_moms.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Mother", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ .) +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

variable.labs <- c("Garmaeva NCs", "Maqsood NCs", "Liang NCs", "Shah NCs")
names(variable.labs) <- c("Perc_shared_w_garmaeva", "Perc_shared_w_maqsood", "Perc_shared_w_liang", "Perc_shared_w_shah")

# Convert named vector to labeller
variable_labeller <- as_labeller(variable.labs)

## Plotting sharedness depending on the timepoints
pdf('shared_timepoints.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant" & table_for_plot_melt$cohort %in% c("garmaeva", "liang"), ], aes(x = Timepoint, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.13, size = 0.5) +
  ylim(0, 100) +
  facet_grid(rows = vars(variable), cols = vars(cohort), scales = "free", labeller = labeller(variable = variable_labeller)) +
  scale_x_discrete() + 
  labs(y = "%", x = "Timepoint") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()
#######################################################################################################################################




##############################################################################################################################################################################################################################################################################


# THIS BLOCK IS DEDICATED TO EXPLORING THE RESULTS, WHEN THE RPKM TABLE IS FILTERED FOR PLASMIDS


##############################################################################################################################################################################################################################################################################


## File preparation
#######################################################################################################################################
virus_selection <- extended_tof$New_CID[extended_tof$plasmid == "No"]
length(virus_selection)  # 948998 (22,585 identified redundant "viruses" are actually plasmids)

RPKM_filtered <- RPKM[row.names(RPKM) %in% virus_selection, ]
RPKM_filtered <- RPKM_filtered[rowSums(RPKM_filtered) > 0, colSums(RPKM_filtered) > 0]
dim(RPKM_filtered)  # 298689   1300

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM_filtered), ]
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
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix

# heatmap.2(x=nc_bray_dist_matrix, dendrogram=F)

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis_no_plasmids.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()


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

RPKM_filtered_w_dummy <- RPKM_filtered

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl"]

temp_shared_w_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_shared_w_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_shared_w_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shared_w_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_shah

temp_shared_w_garmaeva[temp_shared_w_garmaeva < 2] <- 0
temp_shared_w_maqsood[temp_shared_w_maqsood < 2] <- 0
temp_shared_w_liang[temp_shared_w_liang < 2] <- 0
temp_shared_w_shah[temp_shared_w_shah < 2] <- 0

table_for_plot_g <- as.data.frame(colSums(temp_shared_w_garmaeva) / 2)
table_for_plot_m <- as.data.frame(colSums(temp_shared_w_maqsood) / 2)
table_for_plot_l <- as.data.frame(colSums(temp_shared_w_liang) / 2)
table_for_plot_s <- as.data.frame(colSums(temp_shared_w_shah) / 2)

table_for_plot_g_m <- merge(table_for_plot_g, table_for_plot_m, by="row.names", all = T)
table_for_plot_l_s <- merge(table_for_plot_l, table_for_plot_s, by="row.names", all = T)
table_for_plot <- merge(table_for_plot_g_m, table_for_plot_l_s, by="Row.names", all = T)
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL

RPKM_count <- RPKM_filtered
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_add <- as.data.frame(colSums(RPKM_count))
table_for_plot <- merge(table_for_plot, table_for_plot_add, by='row.names', all.x = T)

colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                              "Shared_w_nc_shah", "Viral_richness")

table_for_plot$Perc_shared_w_garmaeva <- table_for_plot$Shared_w_nc_garmaeva / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_maqsood <- table_for_plot$Shared_w_nc_maqsood / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_liang <- table_for_plot$Shared_w_nc_liang / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_shah <- table_for_plot$Shared_w_nc_shah / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                                                   "Shared_w_nc_shah", "Viral_richness", "Type", "Timepoint", "cohort"))

table_for_plot_melt$variable <- factor(table_for_plot_melt$variable, levels = c("Perc_shared_w_garmaeva", "Perc_shared_w_liang", 
                                                                                "Perc_shared_w_maqsood", "Perc_shared_w_shah"))

## Plotting only sharedness for infants
pdf('shared_inf_no_plasmids.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ ., scales = "free") +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

## Plotting only sharendess for mothers
pdf('shared_moms_no_plasmids.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Mother", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ .) +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

variable.labs <- c("Garmaeva NCs", "Maqsood NCs", "Liang NCs", "Shah NCs")
names(variable.labs) <- c("Perc_shared_w_garmaeva", "Perc_shared_w_maqsood", "Perc_shared_w_liang", "Perc_shared_w_shah")

# Convert named vector to labeller
variable_labeller <- as_labeller(variable.labs)

## Plotting sharedness depending on the timepoints
pdf('shared_timepoints_no_plasmids.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant" & table_for_plot_melt$cohort %in% c("garmaeva", "liang"), ], aes(x = Timepoint, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.13, size = 0.5) +
  ylim(0, 100) +
  facet_grid(rows = vars(variable), cols = vars(cohort), scales = "free", labeller = labeller(variable = variable_labeller)) +
  scale_x_discrete() + 
  labs(y = "%", x = "Timepoint") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()





##############################################################################################################################################################################################################################################################################


# THIS BLOCK IS DEDICATED TO EXPLORING THE RESULTS, WHEN THE RPKM TABLE IS FILTERED FOR UNDETERMINED LOW QUALITY VIRUSES AND PLASMIDS


##############################################################################################################################################################################################################################################################################


#######################################################################################################################################
virus_selection <- extended_tof$New_CID[extended_tof$plasmid == "No" & extended_tof$checkv_quality %in% c("High-quality", "Medium-quality", "Complete")]
length(virus_selection)  # 38988 (932,595 identified redundant "viruses" are either viruses of low- undetermined quality or plasmids)

RPKM_filtered <- RPKM[row.names(RPKM) %in% virus_selection, ]
RPKM_filtered <- RPKM_filtered[rowSums(RPKM_filtered) > 0, colSums(RPKM_filtered) > 0]
dim(RPKM_filtered)  # 14519  1185

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM_filtered), ]
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
nc_bray_dist_matrix_rev <- 1 - nc_bray_dist_matrix

# heatmap.2(x=nc_bray_dist_matrix, dendrogram=F)

annotation <- meta_nc_nmds %>%
  mutate(clean_reads_comb = NULL,
         total_viruses_discovered = NULL,
         shannon_div = NULL)

pdf('dist_mat_bray_curtis_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=10/2.54)
pheatmap(nc_bray_dist_matrix_rev, 
         annotation = annotation, 
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = F, 
         color = colorRampPalette(c("white", "darkblue"))(100),
         fontsize = 5,
         main = "Negative controls Bray-Curtis Similarity")
dev.off()

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

RPKM_filtered_w_dummy <- RPKM_filtered

RPKM_filtered_w_dummy$dummy_all <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls])

RPKM_filtered_w_dummy$dummy_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1

samples <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl"]

temp_shared_w_garmaeva <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_garmaeva
temp_shared_w_maqsood <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_maqsood
temp_shared_w_liang <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_liang
temp_shared_w_shah <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% samples] + RPKM_filtered_w_dummy$dummy_shah

temp_shared_w_garmaeva[temp_shared_w_garmaeva < 2] <- 0
temp_shared_w_maqsood[temp_shared_w_maqsood < 2] <- 0
temp_shared_w_liang[temp_shared_w_liang < 2] <- 0
temp_shared_w_shah[temp_shared_w_shah < 2] <- 0

table_for_plot_g <- as.data.frame(colSums(temp_shared_w_garmaeva) / 2)
table_for_plot_m <- as.data.frame(colSums(temp_shared_w_maqsood) / 2)
table_for_plot_l <- as.data.frame(colSums(temp_shared_w_liang) / 2)
table_for_plot_s <- as.data.frame(colSums(temp_shared_w_shah) / 2)

table_for_plot_g_m <- merge(table_for_plot_g, table_for_plot_m, by="row.names", all = T)
table_for_plot_l_s <- merge(table_for_plot_l, table_for_plot_s, by="row.names", all = T)
table_for_plot <- merge(table_for_plot_g_m, table_for_plot_l_s, by="Row.names", all = T)
row.names(table_for_plot) <- table_for_plot$Row.names
table_for_plot$Row.names <- NULL

RPKM_count <- RPKM_filtered
RPKM_count[RPKM_count > 0] <- 1

table_for_plot_add <- as.data.frame(colSums(RPKM_count))
table_for_plot <- merge(table_for_plot, table_for_plot_add, by='row.names', all.x = T)

colnames(table_for_plot) <- c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                              "Shared_w_nc_shah", "Viral_richness")

table_for_plot$Perc_shared_w_garmaeva <- table_for_plot$Shared_w_nc_garmaeva / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_maqsood <- table_for_plot$Shared_w_nc_maqsood / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_liang <- table_for_plot$Shared_w_nc_liang / table_for_plot$Viral_richness *100
table_for_plot$Perc_shared_w_shah <- table_for_plot$Shared_w_nc_shah / table_for_plot$Viral_richness *100

table_for_plot <- merge(table_for_plot, meta_working[c("Sample_name", "Type", "Timepoint", "cohort")], by="Sample_name", all.x = T)

table_for_plot_melt <- melt(table_for_plot, id = c("Sample_name", "Shared_w_nc_garmaeva", "Shared_w_nc_maqsood", "Shared_w_nc_liang", 
                                                   "Shared_w_nc_shah", "Viral_richness", "Type", "Timepoint", "cohort"))

table_for_plot_melt$variable <- factor(table_for_plot_melt$variable, levels = c("Perc_shared_w_garmaeva", "Perc_shared_w_liang", 
                                                                                "Perc_shared_w_maqsood", "Perc_shared_w_shah"))

## Plotting only sharedness for infants
pdf('shared_inf_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ ., scales = "free") +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

## Plotting only sharendess for mothers
pdf('shared_moms_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Mother", ], aes(x=variable, y=value)) +
  geom_violin(trim=FALSE) +
  geom_jitter(width = 0.13, size=0.5) +
  ylim(0, 100) +
  facet_grid(cohort ~ .) +
  scale_x_discrete(labels=c("Perc_shared_w_garmaeva" = "% shared with nc \n from Garmaeva et al",
                            "Perc_shared_w_maqsood" = "% shared with nc \n from Maqsood et al", 
                            "Perc_shared_w_liang" = "% shared with nc \n from Liang et al", 
                            "Perc_shared_w_shah" = "% shared with nc \n from Shah et al")) + 
  labs(y = "%") +
  theme_bw() +
  theme(plot.title = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()

variable.labs <- c("Garmaeva NCs", "Maqsood NCs", "Liang NCs", "Shah NCs")
names(variable.labs) <- c("Perc_shared_w_garmaeva", "Perc_shared_w_maqsood", "Perc_shared_w_liang", "Perc_shared_w_shah")

# Convert named vector to labeller
variable_labeller <- as_labeller(variable.labs)

## Plotting sharedness depending on the timepoints
pdf('shared_timepoints_no_plasmids_no_lowundetermined_quality.pdf', width=12/2.54, height=8/2.54)
ggplot(table_for_plot_melt[table_for_plot_melt$Type == "Infant" & table_for_plot_melt$cohort %in% c("garmaeva", "liang"), ], aes(x = Timepoint, y = value)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.13, size = 0.5) +
  ylim(0, 100) +
  facet_grid(rows = vars(variable), cols = vars(cohort), scales = "free", labeller = labeller(variable = variable_labeller)) +
  scale_x_discrete() + 
  labs(y = "%", x = "Timepoint") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))
dev.off()




















## Analysis for the vOTUs sharedness with different samples
#######################################################################################################################################

test <- RPKM_filtered_w_dummy[RPKM_filtered_w_dummy$dummy_all > 0, ]
test1 <- test[c("dummy_garmaeva", "dummy_liang", "dummy_maqsood", "dummy_shah")]
test1$shared_neg_ctrl <- rowSums(test1)
test1$shared_all_samples <- rowSums(test[colnames(test) %in% 
                                           meta_working$Sample_name[meta_working$Type != "Neg_ctrl"] ])
test1$shared_infants <- rowSums(test[colnames(test) %in% 
                                       meta_working$Sample_name[meta_working$Type == "Infant"] ])
test1$shared_mothers <- rowSums(test[colnames(test) %in% 
                                       meta_working$Sample_name[meta_working$Type == "Mother"] ])


test1$shared_infants_l <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "liang"] ])

test1$shared_infants_s <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "shah"] ])

test1$shared_infants_m <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "maqsood"] ])

test1$shared_infants_g <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Infant" & meta_working$cohort == "garmaeva"] ])

test1$shared_mother_m <- rowSums(test[colnames(test) %in% 
                                        meta_working$Sample_name[meta_working$Type == "Mother" & meta_working$cohort == "maqsood"] ])

test1$shared_mothers_g <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type == "Mother" & meta_working$cohort == "garmaeva"] ])


test1$shared_samples_l <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "liang"] ])

test1$shared_samples_s <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "shah"] ])

test1$shared_samples_m <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "maqsood"] ])

test1$shared_samples_g <- rowSums(test[colnames(test) %in% 
                                         meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "garmaeva"] ])



## To do
#######################################################################################################################################
# add the savings for the plot
# add walters to the last analysis
# now I only have sharedness with the negative controls from other cohorts, but I need one more for sharedness with the same cohort
# color the dots in the violin plot by timepoint
#######################################################################################################################################


# Redo all the analysis for the following










