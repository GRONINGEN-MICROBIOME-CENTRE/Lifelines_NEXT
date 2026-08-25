# ─────────────────────────────────────────────
# Import libraries
# ─────────────────────────────────────────────

library(reshape2)
library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggbeeswarm)
library(foreach)

# ─────────────────────────────────────────────
# Import plant FoodSeq data
# ─────────────────────────────────────────────

asv_binary <- read.csv("foodseq_plant_asv_binary.csv", row.names = 1)
metadata_selected <- read.csv("foodseq_metadata_clean.csv")


# ─────────────────────────────────────────────
# Plant PCoA
# ─────────────────────────────────────────────

asv_bin2 = asv_binary

asv_bin2 = asv_binary[rowSums(asv_binary)>0,]
dist1 = vegdist(asv_bin2,method = "jaccard",binary= T)
data.pcoa = capscale(dist1 ~ 0)

lvls <- c("Infant (M9)", "Infant (M12)", "Mother")

idx  <- match(rownames(asv_bin2), metadata_selected$NG_ID)  # positions in metadata
keep <- !is.na(idx)  

metadata_selected <- metadata_selected[idx[keep], , drop = FALSE]

metadata_selected$Timepoint_categorical[metadata_selected$Timepoint_categorical == "M3"] <- "Mother"
metadata_selected$Timepoint_categorical[metadata_selected$Timepoint_categorical == "M9"] <- "Infant (M9)"
metadata_selected$Timepoint_categorical[metadata_selected$Timepoint_categorical == "M12"] <- "Infant (M12)"

site_scores <- scores(data.pcoa, display = "sites") %>% as.data.frame()
site_scores$SampleID <- rownames(site_scores)

eig  <- data.pcoa$CA$eig
varp <- eig / sum(eig)
pc1_lab <- paste0("PCoA1 (", round(100 * varp[1], 1), "%)")
pc2_lab <- paste0("PCoA2 (", round(100 * varp[2], 1), "%)")

# Join scores to metadata
plot_df <- site_scores %>%
  cbind(metadata_selected %>% select(Timepoint_categorical))

fill_pal <- c("Infant (M9)"="#A6D8FF","Infant (M12)"="#F2C77C","Mother"="#66C2A3")
line_pal <- c("Infant (M9)"="#2E7DB2","Infant (M12)"="#B77700","Mother"="#006D5B")

p_pca <- ggplot(plot_df, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    type        = "t",
    geom        = "polygon",
    alpha       = 0.18,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = Timepoint_categorical,
        color = Timepoint_categorical),
    size   = 2,
    alpha  = 1,
    shape  = 21,
    stroke = 0.7
  ) +
  coord_equal() +
  scale_fill_manual(values = fill_pal,  breaks = lvls) +
  scale_color_manual(values = line_pal, breaks = lvls) +
  labs(
    x     = pc1_lab,
    y     = pc2_lab,
    color = "Sample group",
    fill  = "Sample group",
    title = "Plant"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.4),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length  = unit(2.2, "mm"),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    legend.position    = "right",
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.margin        = margin(t = 16, r = 14, b = 14, l = 10)
  )

p_out <- cowplot::ggdraw(p_pca) +
  cowplot::draw_plot_label(
    label = "A", x = 0, y = 1, hjust = -0.5, vjust = 1.5,
    fontface = "bold", size = 20
  )

p_out

ggsave("plant_pcoa.tiff", plot = p_out,
       dpi = 600, width = 6, height = 4.5, units = "in",
       compression = "lzw", bg = "white")   # <-- force white background


# ─────────────────────────────────────────────
# PERMANOVA + betadisper
# ─────────────────────────────────────────────

# Mother and infant M9

mothers_babiesM9 <- metadata_selected$Timepoint_categorical %in% c("Mother", "Infant (M9)")

asv_m9 <- asv_bin2[mothers_babiesM9, ]
meta_m9 <- metadata_selected[mothers_babiesM9, ]

distances_m9 <- vegdist(asv_m9, method = "jaccard", binary = T)

adonis_m9 <- adonis2(
  distances_m9 ~ Timepoint_categorical,
  data = meta_m9,
  permutations = 999)
adonis_m9

betadisper_m9 <- betadisper(distances_m9, meta_m9$Timepoint_categorical)
anova(betadisper_m9)
# permutest(betadisper_m9, permutations = 999)


# Mother and infant M12
mothers_babiesM12 <- metadata_selected$Timepoint_categorical %in% c("Mother", "Infant (M12)")

asv_m12 <- asv_bin2[mothers_babiesM12, ]
meta_m12 <- metadata_selected[mothers_babiesM12, ]

distances_m12 <- vegdist(asv_m12, method = "jaccard", binary = T)

adonis_m12 <- adonis2(
  distances_m12 ~ Timepoint_categorical,
  data = meta_m12,
  permutations = 999)
adonis_m12

betadisper_m12 <- betadisper(distances_m12, meta_m12$Timepoint_categorical)
anova(betadisper_m12)
# permutest(betadisper_m12, permutations = 999)


# ─────────────────────────────────────────────
# Family distances (plant)
# ─────────────────────────────────────────────

asv_bin2 <- asv_binary[rowSums(asv_binary) > 0, ]

d1 = vegdist(asv_bin2, method = "jaccard", binary = TRUE)
d1 <- as.matrix(d1)

subset.matrix = d1[
  colnames(d1) %in% with(metadata_selected, NG_ID[Type == "mother"]),
  rownames(d1) %in% with(metadata_selected, NG_ID[Type == "infant"])]
subset.matrix.long = melt(subset.matrix)
subset.matrix.long$baby_FID = metadata_selected[match(subset.matrix.long$Var1, metadata_selected$NG_ID),"FAMILY"]
subset.matrix.long$mom_FID = metadata_selected[match(subset.matrix.long$Var2, metadata_selected$NG_ID),"FAMILY"]

values_sameFAM <- subset.matrix.long$value[subset.matrix.long$baby_FID == subset.matrix.long$mom_FID]

values_differentFAM <- subset.matrix.long$value[subset.matrix.long$baby_FID != subset.matrix.long$mom_FID]

# Unadjusted pvalue
pval = wilcox.test(values_sameFAM,values_differentFAM)$p.value

# Permutated pvalue
set.seed(123)
perm.P <- foreach(i = 1:1000, .combine = c) %do% {
  subset.matrix.perm = subset.matrix[,sample(1:ncol(subset.matrix))]
  subset.matrix.perm.long = melt(subset.matrix.perm)
  
  values_sameFAM_perm <- subset.matrix.perm.long[subset.matrix.long$baby_FID == subset.matrix.long$mom_FID, "value"]
  values_differentFAM_perm <- subset.matrix.perm.long[subset.matrix.long$baby_FID != subset.matrix.long$mom_FID, "value"]
  
  wilcox.test(values_sameFAM_perm, values_differentFAM_perm)$p.value
  }

test.real.p = sum(perm.P<=pval) / 1000 
test.real.p


#### Plot #### 

distances_df <- rbind(
  data.frame(Family_status = "sameFAM",      Distance = values_sameFAM),
  data.frame(Family_status = "differentFAM", Distance = values_differentFAM)
)
distances_df$Family_status <- factor(distances_df$Family_status,
                                     levels = c("sameFAM", "differentFAM"))

fill_pal <- c("sameFAM" = "#A6D8FF", "differentFAM" = "#F4A6A6")
line_pal <- c("sameFAM" = "#2E7DB2", "differentFAM" = "#C52233")

p_label_df <- data.frame(
  group1     = "sameFAM",
  group2     = "differentFAM",
  y.position = max(distances_df$Distance, na.rm = TRUE) * 1.08,
  label      = paste0("Perm. p = ", signif(test.real.p, 3))
)

p <- ggplot(distances_df, aes(x = Family_status, y = Distance)) +
  geom_boxplot(
    aes(fill = Family_status, color = Family_status),
    width = 0.55, outlier.shape = NA, size = 0.7, fatten = 1.15, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Family_status, color = Family_status),
    groupOnX = TRUE, width = 0.05,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = line_pal) +
  stat_pvalue_manual(p_label_df, label = "label", tip.length = 0.015, size = 3.6) +
  labs(x = "Family", y = "Distance value", title = "Plant") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  scale_x_discrete(labels = c(
    "sameFAM"      = "Within family",
    "differentFAM" = "Between families"
  )) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.4),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length  = unit(2.2, "mm"),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.margin        = margin(t = 16, r = 14, b = 14, l = 10)
  )

p_out <- ggdraw(p) +
  draw_plot_label(label = "A", x = 0, y = 1,
                  hjust = -0.5, vjust = 1.5,
                  fontface = "bold", size = 20)

p_out

ggsave(
  filename = "plant_permutations.tiff",
  plot = p_out,
  dpi = 600, width = 4, height = 4.5, units = "in",
  compression = "lzw",
  bg = "white"
)


## IQR same family
iqr_same_family <- c(median=median(values_sameFAM), 
                     Q1=quantile(values_sameFAM, 0.25), 
                     Q3=quantile(values_sameFAM, 0.75))
iqr_same_family

## IQR different family
iqr_diff_family <- c(median=median(values_differentFAM), 
                     Q1=quantile(values_differentFAM, 0.25), 
                     Q3=quantile(values_differentFAM, 0.75))
iqr_diff_family


# ─────────────────────────────────────────────
# ─────────────────────────────────────────────
# Animal FoodSeq data
# ─────────────────────────────────────────────
# ─────────────────────────────────────────────

rm(list=ls())

# Import data
asv_binary <- read.csv("foodseq_animal_asv_binary.csv", row.names = 1)
metadata_selected <- read.csv("foodseq_metadata_clean.csv")


# ─────────────────────────────────────────────
# Animal PCoA
# ─────────────────────────────────────────────

asv_bin2 = asv_binary[rowSums(asv_binary)>0,]
dist1 = vegdist(asv_bin2,method = "jaccard",binary= T)
data.pcoa = capscale(dist1 ~ 0)

lvls <- c("Infant (M9)", "Infant (M12)", "Mother")

idx  <- match(rownames(asv_bin2), metadata_selected$NG_ID)  # positions in metadata
keep <- !is.na(idx)  

metadata_selected <- metadata_selected[idx[keep], , drop = FALSE]

metadata_selected$Timepoint_categorical[metadata_selected$Timepoint_categorical == "M3"] <- "Mother"
metadata_selected$Timepoint_categorical[metadata_selected$Timepoint_categorical == "M9"] <- "Infant (M9)"
metadata_selected$Timepoint_categorical[metadata_selected$Timepoint_categorical == "M12"] <- "Infant (M12)"

site_scores <- scores(data.pcoa, display = "sites") %>% as.data.frame()
site_scores$SampleID <- rownames(site_scores)

eig  <- data.pcoa$CA$eig
varp <- eig / sum(eig)
pc1_lab <- paste0("PCoA1 (", round(100 * varp[1], 1), "%)")
pc2_lab <- paste0("PCoA2 (", round(100 * varp[2], 1), "%)")

# Join scores to metadata
plot_df <- site_scores %>%
  cbind(metadata_selected %>% select(Timepoint_categorical))

fill_pal <- c("Infant (M9)"="#A6D8FF","Infant (M12)"="#F2C77C","Mother"="#66C2A3")
line_pal <- c("Infant (M9)"="#2E7DB2","Infant (M12)"="#B77700","Mother"="#006D5B")

p_pca <- ggplot(plot_df, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    type        = "t",
    geom        = "polygon",
    alpha       = 0.18,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = Timepoint_categorical,
        color = Timepoint_categorical),
    size   = 2,
    alpha  = 1,
    shape  = 21,
    stroke = 0.7
  ) +
  coord_equal() +
  scale_fill_manual(values = fill_pal,  breaks = lvls) +
  scale_color_manual(values = line_pal, breaks = lvls) +
  labs(
    x     = pc1_lab,
    y     = pc2_lab,
    color = "Sample group",
    fill  = "Sample group",
    title = "Animal"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.4),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length  = unit(2.2, "mm"),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    legend.position    = "right",
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.margin        = margin(t = 16, r = 14, b = 14, l = 10)
  )

p_out <- cowplot::ggdraw(p_pca) +
  cowplot::draw_plot_label(
    label = "B", x = 0, y = 1, hjust = -0.5, vjust = 1.5,
    fontface = "bold", size = 20
  )

p_out

ggsave("animal_pcoa.tiff", plot = p_out,
       dpi = 600, width = 6, height = 4.5, units = "in",
       compression = "lzw", bg = "white")


# ─────────────────────────────────────────────
# PERMANOVA + betadisper
# ─────────────────────────────────────────────

# Mother and infant M9

mothers_babiesM9 <- metadata_selected$Timepoint_categorical %in% c("Mother", "Infant (M9)")

asv_m9 <- asv_bin2[mothers_babiesM9, ]
meta_m9 <- metadata_selected[mothers_babiesM9, ]

distances_m9 <- vegdist(asv_m9, method = "jaccard", binary = T)

adonis_m9 <- adonis2(
  distances_m9 ~ Timepoint_categorical,
  data = meta_m9,
  permutations = 999)
adonis_m9

betadisper_m9 <- betadisper(distances_m9, meta_m9$Timepoint_categorical)
anova(betadisper_m9)
# permutest(betadisper_m9, permutations = 999)


# Mother and infant M12
mothers_babiesM12 <- metadata_selected$Timepoint_categorical %in% c("Mother", "Infant (M12)")

asv_m12 <- asv_bin2[mothers_babiesM12, ]
meta_m12 <- metadata_selected[mothers_babiesM12, ]

distances_m12 <- vegdist(asv_m12, method = "jaccard", binary = T)

adonis_m12 <- adonis2(
  distances_m12 ~ Timepoint_categorical,
  data = meta_m12,
  permutations = 999)
adonis_m12

betadisper_m12 <- betadisper(distances_m12, meta_m12$Timepoint_categorical)
anova(betadisper_m12)
# permutest(betadisper_m12, permutations = 999)


# ─────────────────────────────────────────────
# Family distances (animal)
# ─────────────────────────────────────────────

asv_bin2 <- asv_binary[rowSums(asv_binary) > 0, ]

d1 = vegdist(asv_bin2, method = "jaccard", binary = TRUE)
d1 <- as.matrix(d1)

subset.matrix = d1[
  colnames(d1) %in% with(metadata_selected, NG_ID[Type == "mother"]),
  rownames(d1) %in% with(metadata_selected, NG_ID[Type == "infant"])]
subset.matrix.long = melt(subset.matrix)
subset.matrix.long$baby_FID = metadata_selected[match(subset.matrix.long$Var1, metadata_selected$NG_ID),"FAMILY"]
subset.matrix.long$mom_FID = metadata_selected[match(subset.matrix.long$Var2, metadata_selected$NG_ID),"FAMILY"]

values_sameFAM <- subset.matrix.long$value[subset.matrix.long$baby_FID == subset.matrix.long$mom_FID]

values_differentFAM <- subset.matrix.long$value[subset.matrix.long$baby_FID != subset.matrix.long$mom_FID]

# Unadjusted pvalue
pval = wilcox.test(values_sameFAM,values_differentFAM)$p.value

# Permutated pvalue
set.seed(123)
perm.P <- foreach(i = 1:1000, .combine = c) %do% {
  subset.matrix.perm = subset.matrix[,sample(1:ncol(subset.matrix))]
  subset.matrix.perm.long = melt(subset.matrix.perm)
  
  values_sameFAM_perm <- subset.matrix.perm.long[subset.matrix.long$baby_FID == subset.matrix.long$mom_FID, "value"]
  values_differentFAM_perm <- subset.matrix.perm.long[subset.matrix.long$baby_FID != subset.matrix.long$mom_FID, "value"]
  
  wilcox.test(values_sameFAM_perm, values_differentFAM_perm)$p.value
}

test.real.p = sum(perm.P<=pval) / 1000 
test.real.p



#### Plot #### 

distances_df <- rbind(
  data.frame(Family_status = "sameFAM",      Distance = values_sameFAM),
  data.frame(Family_status = "differentFAM", Distance = values_differentFAM)
)
distances_df$Family_status <- factor(distances_df$Family_status,
                                     levels = c("sameFAM", "differentFAM"))

fill_pal <- c("sameFAM" = "#A6D8FF", "differentFAM" = "#F4A6A6")
line_pal <- c("sameFAM" = "#2E7DB2", "differentFAM" = "#C52233")

p_label_df <- data.frame(
  group1     = "sameFAM",
  group2     = "differentFAM",
  y.position = max(distances_df$Distance, na.rm = TRUE) * 1.08,
  label      = paste0("Perm. p = ", signif(test.real.p, 3))
)

p <- ggplot(distances_df, aes(x = Family_status, y = Distance)) +
  geom_boxplot(
    aes(fill = Family_status, color = Family_status),
    width = 0.55, outlier.shape = NA, size = 0.7, fatten = 1.15, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Family_status, color = Family_status),
    groupOnX = TRUE, width = 0.1,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = line_pal) +
  stat_pvalue_manual(p_label_df, label = "label", tip.length = 0.015, size = 3.6) +
  labs(x = "Family", y = "Distance value", title = "Animal") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  scale_x_discrete(labels = c(
    "sameFAM"      = "Within family",
    "differentFAM" = "Between families"
  )) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.4),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length  = unit(2.2, "mm"),
    axis.title.x       = element_text(margin = margin(t = 8)),
    axis.title.y       = element_text(margin = margin(r = 8)),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    plot.margin        = margin(t = 16, r = 14, b = 14, l = 10)
  )


p_out <- ggdraw(p) +
  draw_plot_label(label = "B", x = 0, y = 1,
                  hjust = -0.5, vjust = 1.5,
                  fontface = "bold", size = 20)

p_out

ggsave(
  filename = "animal_permutations.tiff",
  plot = p_out,
  dpi = 600, width = 4, height = 4.5, units = "in",
  compression = "lzw",
  bg = "white"
)


## IQR same family
iqr_same_family <- c(median=median(values_sameFAM), 
                     Q1=quantile(values_sameFAM, 0.25), 
                     Q3=quantile(values_sameFAM, 0.75))
iqr_same_family

## IQR different family
iqr_diff_family <- c(median=median(values_differentFAM), 
                     Q1=quantile(values_differentFAM, 0.25), 
                     Q3=quantile(values_differentFAM, 0.75))
iqr_diff_family
