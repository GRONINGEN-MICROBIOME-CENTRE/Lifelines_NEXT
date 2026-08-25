# ─────────────────────────────────────────────
# Import libraries
# ─────────────────────────────────────────────

library(vegan)
library(rstatix)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)
library(coin)

# ─────────────────────────────────────────────
# Import PLANT data
# ─────────────────────────────────────────────

asv_table <- read.csv("foodseq_plant_asv_raw.csv", row.names = 1)
metadata_selected <- read.csv("foodseq_metadata_clean.csv")

# ─────────────────────────────────────────────
# Shannon diversity for plant data
# ─────────────────────────────────────────────

# Relative abundance (proportion)
asv_table <- decostand(asv_table, MARGIN=1, method="total")

metadata_selected$Timepoint_categorical <- factor(metadata_selected$Timepoint_categorical, levels = c("M9", "M12", "M3"),
                                                  labels = c("Infant (M9)", "Infant (M12)", "Mother"))

shannon_asv <- as.data.frame(diversity(asv_table, index = "shannon", MARGIN = 1))
colnames(shannon_asv)[1] <- "shannon"
shannon_asv$NG_ID <- rownames(shannon_asv)

shannon_merged <- left_join(shannon_asv, metadata_selected, by = "NG_ID")

lvls <- c("Infant (M9)","Infant (M12)","Mother")
shannon_merged$Timepoint_categorical <- factor(shannon_merged$Timepoint_categorical, levels = lvls)


# ─────────────────────────────────────────────
# Wilcox & BH correct for Shannon diversity
# ─────────────────────────────────────────────

compare_groups <- list(c("Infant (M9)", "Infant (M12)"), c("Infant (M9)", "Mother"), c("Infant (M12)", "Mother"))

shannon_tests <- pairwise_wilcox_test(
  shannon_merged,
  shannon ~ Timepoint_categorical,
  comparisons = compare_groups,
  p.adjust.method = "BH"
  ) %>%
  mutate(y.position = c(2.5, 2.75, 3),
         p.adj.label = paste0("adj. p = ", signif(p.adj)))


# Wilcox effect size (r)
shannon_effects <- wilcox_effsize(shannon_merged,
                                  shannon ~ Timepoint_categorical,
                                  comparisons = compare_groups)

shannon_results <- left_join(shannon_tests, shannon_effects, by = c("group1", "group2"))

shannon_results


# ─────────────────────────────────────────────
# Plot Shannon diversity for plants
# ─────────────────────────────────────────────

fill_pal <- c("Infant (M9)"="#A6D8FF","Infant (M12)"="#F2C77C","Mother"="#66C2A3")
line_pal <- c("Infant (M9)"="#2E7DB2","Infant (M12)"="#B77700","Mother"="#006D5B")

p <- ggplot(shannon_merged, aes(Timepoint_categorical, shannon)) +
  geom_boxplot(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    width = 0.55, outlier.shape = NA, size = 0.7, median.linewidth = 0.8, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    groupOnX = TRUE, width = 0.1, varwidth = FALSE,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  stat_pvalue_manual(
    shannon_tests,
    label = "p.adj.label",
    size = 3.6,
    tip.length = 0.01
  ) +
  scale_fill_manual(values=fill_pal, breaks=lvls) +
  scale_color_manual(values=line_pal, breaks=lvls) +
  labs(x="Sample group", y="Shannon Diversity Index", title="Plant") +
  scale_y_continuous(expand=expansion(mult=c(.02,.14))) +
  coord_cartesian(ylim=c(0,3),clip="off") +
  theme_minimal(14) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="#EFEAE3", linewidth=.4),
    axis.line=element_line(color="black", linewidth=.5),
    axis.ticks=element_line(color="black", linewidth=.5),
    legend.position="none",
    plot.title=element_text(hjust=.5, face="bold"),
    plot.margin=margin(t=16, r=14, b=14, l=10))

p_out <- ggdraw(p) + draw_plot_label("A", x=0, y=1, hjust=-.5, vjust=1.5, fontface="bold", size=20)

p_out

ggsave("plant_shannon.tiff", plot = p_out, dpi = 600, width = 4, height = 4.5, units = "in", compression = "lzw", bg = "white")



# ─────────────────────────────────────────────
# Simpson diversity for plant data
# ─────────────────────────────────────────────

simpson_asv <- as.data.frame(diversity(asv_table, index = "simpson"))
colnames(simpson_asv)[1] <- "simpson"
simpson_asv$NG_ID <- rownames(simpson_asv)

simpson_merged <- left_join(simpson_asv, metadata_selected, by = "NG_ID")

simpson_merged$Timepoint_categorical <- factor(simpson_merged$Timepoint_categorical, levels = lvls)


# ─────────────────────────────────────────────
# Wilcox & BH correct for Simpson diversity
# ─────────────────────────────────────────────

simpson_tests <- pairwise_wilcox_test(
  simpson_merged,
  simpson ~ Timepoint_categorical,
  comparisons = compare_groups,
  p.adjust.method = "BH"
) %>%
  mutate(y.position = c(0.89, 0.99, 1.09),
         p.adj.label = paste0("adj. p = ", signif(p.adj)))


# Wilcox effect size (r)
simpson_effects <- wilcox_effsize(simpson_merged,
                                  simpson ~ Timepoint_categorical,
                                  comparisons = compare_groups)

simpson_results <- left_join(simpson_tests, simpson_effects, by = c("group1", "group2"))

simpson_results


# ─────────────────────────────────────────────
# Plot Simpson diversity for plants
# ─────────────────────────────────────────────

p <- ggplot(simpson_merged, aes(Timepoint_categorical, simpson)) +
  geom_boxplot(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    width = 0.55, outlier.shape = NA, size = 0.7, median.linewidth = 0.8, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    groupOnX = TRUE, width = 0.1, varwidth = FALSE,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  scale_fill_manual(values=fill_pal, breaks=lvls) +
  scale_color_manual(values=line_pal, breaks=lvls) +
  stat_pvalue_manual(
    simpson_tests,
    label = "p.adj.label",
    size = 3.6,
    tip.length = 0.01
  ) +
  labs(x="Sample group", y="Simpson's Diversity Index", title="Plant") +
  scale_y_continuous(expand=expansion(mult=c(.02,.14))) +
  coord_cartesian(clip="off") +
  theme_minimal(14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="#EFEAE3", linewidth=.4),
        axis.line=element_line(color="black", linewidth=.5),
        axis.ticks=element_line(color="black", linewidth=.5),
        legend.position="none",
        plot.title=element_text(hjust=.5, face="bold"),
        plot.margin=margin(t=16, r=14, b=14, l=10))

p_out <- ggdraw(p) + draw_plot_label("B", x=0, y=1, hjust=-.5, vjust=1.5,
                                     fontface="bold", size=20)
p_out

ggsave("plant_simpson.tiff", plot = p_out, dpi = 600, width = 4, height = 4.5, units = "in", compression = "lzw", bg = "white")



# ─────────────────────────────────────────────
# ─────────────────────────────────────────────
# Import ANIMAL data
# ─────────────────────────────────────────────
# ─────────────────────────────────────────────

rm(list=ls())

asv_table <- read.csv("foodseq_animal_asv_raw.csv", row.names = 1)
metadata_selected <- read.csv("foodseq_metadata_clean.csv")

# ─────────────────────────────────────────────
# Shannon diversity for animal data
# ─────────────────────────────────────────────

asv_table <- asv_table[rowSums(asv_table) > 0, ]

# Relative abundance (proportion)
asv_table<-decostand(asv_table, MARGIN=1, method="total")


metadata_selected$Timepoint_categorical <- factor(metadata_selected$Timepoint_categorical, levels = c("M9", "M12", "M3"),
                                                  labels = c("Infant (M9)", "Infant (M12)", "Mother"))

shannon_asv <- as.data.frame(diversity(asv_table, index = "shannon", MARGIN = 1))
colnames(shannon_asv)[1] <- "shannon"
shannon_asv$NG_ID <- rownames(shannon_asv)

shannon_merged <- left_join(shannon_asv, metadata_selected, by = "NG_ID")

# To check the numbers of samples in each group
shannon_merged %>% count(Timepoint_categorical)

lvls <- c("Infant (M9)","Infant (M12)","Mother")
shannon_merged$Timepoint_categorical <- factor(shannon_merged$Timepoint_categorical, levels = lvls)


# ─────────────────────────────────────────────
# Wilcox & BH correct for Shannon diversity
# ─────────────────────────────────────────────

compare_groups <- list(c("Infant (M9)", "Infant (M12)"), c("Infant (M9)", "Mother"), c("Infant (M12)", "Mother"))

shannon_tests <- pairwise_wilcox_test(
  shannon_merged,
  shannon ~ Timepoint_categorical,
  comparisons = compare_groups,
  p.adjust.method = "BH"
) %>%
  mutate(y.position = c(1.22, 1.35, 1.48),
         p.adj.label = paste0("adj. p = ", signif(p.adj)))


# Wilcox effect size (r)
shannon_effects <- wilcox_effsize(shannon_merged,
                                  shannon ~ Timepoint_categorical,
                                  comparisons = compare_groups)

shannon_results <- left_join(shannon_tests, shannon_effects, by = c("group1", "group2"))

shannon_results

# ─────────────────────────────────────────────
# Plot Shannon diversity for animal data
# ─────────────────────────────────────────────

fill_pal <- c("Infant (M9)"="#A6D8FF","Infant (M12)"="#F2C77C","Mother"="#66C2A3")
line_pal <- c("Infant (M9)"="#2E7DB2","Infant (M12)"="#B77700","Mother"="#006D5B")

p <- ggplot(shannon_merged, aes(Timepoint_categorical, shannon)) +
  geom_boxplot(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    width = 0.55, outlier.shape = NA, size = 0.7, median.linewidth = 0.8, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    groupOnX = TRUE, width = 0.1, varwidth = FALSE,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  stat_pvalue_manual(
    shannon_tests,
    label = "p.adj.label",
    size = 3.6,
    tip.length = 0.01
  ) +
  scale_fill_manual(values=fill_pal, breaks=lvls) +
  scale_color_manual(values=line_pal, breaks=lvls) +
  labs(x="Sample group", y="Shannon Diversity Index", title="Animal") +
  scale_y_continuous(expand=expansion(mult=c(.02,.14))) +
  coord_cartesian(clip="off") +
  theme_minimal(14) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(color="#EFEAE3", linewidth=.4),
    axis.line=element_line(color="black", linewidth=.5),
    axis.ticks=element_line(color="black", linewidth=.5),
    legend.position="none",
    plot.title=element_text(hjust=.5, face="bold"),
    plot.margin=margin(t=16, r=14, b=14, l=10))

p_out <- ggdraw(p) + draw_plot_label("C", x=0, y=1, hjust=-.5, vjust=1.5, fontface="bold", size=20)

p_out

ggsave("animal_shannon.tiff", plot = p_out, dpi = 600, width = 4, height = 4.5, units = "in", compression = "lzw", bg = "white")


# ─────────────────────────────────────────────
# Simpson diversity for animal data
# ─────────────────────────────────────────────

simpson_asv <- as.data.frame(diversity(asv_table, index = "simpson"))
colnames(simpson_asv)[1] <- "simpson"
simpson_asv$NG_ID <- rownames(simpson_asv)

simpson_merged <- left_join(simpson_asv, metadata_selected, by = "NG_ID")

simpson_merged$Timepoint_categorical <- factor(simpson_merged$Timepoint_categorical, levels = lvls)


# ─────────────────────────────────────────────
# Wilcox & BH correct for Simpson diversity
# ─────────────────────────────────────────────

simpson_tests <- pairwise_wilcox_test(
  simpson_merged,
  simpson ~ Timepoint_categorical,
  comparisons = compare_groups,
  p.adjust.method = "BH"
) %>%
  mutate(y.position = c(0.73, 0.82, 0.91),
         p.adj.label = paste0("adj. p = ", signif(p.adj)))


# Wilcox effect size (r)
simpson_effects <- wilcox_effsize(simpson_merged,
                                  simpson ~ Timepoint_categorical,
                                  comparisons = compare_groups)

simpson_results <- left_join(simpson_tests, simpson_effects, by = c("group1", "group2"))

simpson_results


# ─────────────────────────────────────────────
# Plot Simpson diversity for animal data
# ─────────────────────────────────────────────

p <- ggplot(simpson_merged, aes(Timepoint_categorical, simpson)) +
  geom_boxplot(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    width = 0.55, outlier.shape = NA, size = 0.7, median.linewidth = 0.8, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    groupOnX = TRUE, width = 0.1, varwidth = FALSE,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  scale_fill_manual(values=fill_pal, breaks=lvls) +
  scale_color_manual(values=line_pal, breaks=lvls) +
  stat_pvalue_manual(
    simpson_tests,
    label = "p.adj.label",
    size = 3.6,
    tip.length = 0.01
  ) +
  labs(x="Sample group", y="Simpson's Diversity Index", title="Animal") +
  scale_y_continuous(expand=expansion(mult=c(.02,.14))) +
  coord_cartesian(clip="off") +
  theme_minimal(14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color="#EFEAE3", linewidth=.4),
        axis.line=element_line(color="black", linewidth=.5),
        axis.ticks=element_line(color="black", linewidth=.5),
        legend.position="none",
        plot.title=element_text(hjust=.5, face="bold"),
        plot.margin=margin(t=16, r=14, b=14, l=10))

p_out <- ggdraw(p) + draw_plot_label("D", x=0, y=1, hjust=-.5, vjust=1.5,
                                     fontface="bold", size=20)
p_out

ggsave("animal_simpson.tiff", plot = p_out, dpi = 600, width = 4, height = 4.5, units = "in", compression = "lzw", bg = "white")

