
# ─────────────────────────────────────────────
# Import libraries
# ─────────────────────────────────────────────

library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)
library(grid)
library(reshape2)
library(rstatix)
library(coin)

# ─────────────────────────────────────────────
# Import PLANT data
# ─────────────────────────────────────────────

foodseq_wide <- read.csv("foodseq_plant_asv_raw.csv", row.names = 1)
metadata_selected <- read.csv("foodseq_metadata_clean.csv")

# ─────────────────────────────────────────────
# Wilcoxon rank-sum tests with BH
# ─────────────────────────────────────────────

metadata_selected$Timepoint_categorical <- factor(metadata_selected$Timepoint_categorical, levels = c("M9", "M12", "M3"),
                                                  labels = c("Infant (M9)", "Infant (M12)", "Mother"))

foodseq_wide$raw_counts <- rowSums(foodseq_wide)
foodseq_wide$NG_ID <- rownames(foodseq_wide)

foodseq_merged <- left_join(foodseq_wide, metadata_selected, by = "NG_ID")


lvls <- c("Infant (M9)", "Infant (M12)", "Mother")
foodseq_merged$Timepoint_categorical <- factor(foodseq_merged$Timepoint_categorical, levels = lvls)

compare_groups <- list(c("Infant (M9)", "Infant (M12)"), c("Infant (M9)", "Mother"), c("Infant (M12)", "Mother"))

# Wilcoxon rank-sum tests with BH correction
raw_count_wilcox <- pairwise_wilcox_test(
  foodseq_merged,
  raw_counts ~ Timepoint_categorical,
  comparisons = compare_groups,
  p.adjust.method = "BH"
) %>%
  mutate(
    p.adj.label = paste0("adj. p = ", signif(p.adj)),
    y.position = c(90000, 100000, 110000)
  )


# ─────────────────────────────────────────────
# Plot raw sum counts for plant data
# ─────────────────────────────────────────────

fill_pal <- c("Infant (M9)"="#A6D8FF","Infant (M12)"="#F2C77C","Mother"="#66C2A3")
line_pal <- c("Infant (M9)"="#2E7DB2","Infant (M12)"="#B77700","Mother"="#006D5B")

p <- ggplot(
  foodseq_merged,
  aes(x = Timepoint_categorical, y = raw_counts)
) +
  geom_boxplot(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    width = 0.55, outlier.shape = NA, size = 0.7, fatten = 1.15, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    groupOnX = TRUE, width = 0.1, varwidth = FALSE,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  stat_pvalue_manual(
    raw_count_wilcox,
    label = "p.adj.label",
    size = 3.6,
    tip.length = 0.01
  ) +
  scale_fill_manual(values = fill_pal, breaks = lvls) +
  scale_color_manual(values = line_pal,  breaks = lvls) +
  labs(
    x = "Sample group",
    y = "Raw Read Sum",
    title = "Plant"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.4),
    axis.line  = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(2.2, "mm"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 16, r = 14, b = 14, l = 10)
  )

p_out <- ggdraw(p) +
  draw_plot_label(label = "A", x = 0, y = 1, hjust = -0.5, vjust = 1.5,
                  fontface = "bold", size = 20)

p_out

ggsave(
  filename = "plant_foodseq_rawcounts.tiff",
  plot = p_out,
  dpi = 600, width = 4, height = 4.5, units = "in",
  compression = "lzw",
  bg = "white"
)


# ─────────────────────────────────────────────
# ─────────────────────────────────────────────
# Import ANIMAL data
# ─────────────────────────────────────────────
# ─────────────────────────────────────────────

rm(list=ls())

foodseq_wide <- read.csv("foodseq_animal_asv_raw.csv", row.names = 1)
metadata_selected <- read.csv("foodseq_metadata_clean.csv")

# ─────────────────────────────────────────────
# Wilcoxon rank-sum tests with BH
# ─────────────────────────────────────────────

foodseq_wide <- foodseq_wide[rowSums(foodseq_wide) > 0, ]

metadata_selected$Timepoint_categorical <- factor(metadata_selected$Timepoint_categorical, levels = c("M9", "M12", "M3"),
                                                  labels = c("Infant (M9)", "Infant (M12)", "Mother"))

foodseq_wide$raw_counts <- rowSums(foodseq_wide)
foodseq_wide$NG_ID <- rownames(foodseq_wide)

foodseq_merged <- left_join(foodseq_wide, metadata_selected, by = "NG_ID")


lvls <- c("Infant (M9)", "Infant (M12)", "Mother")
foodseq_merged$Timepoint_categorical <- factor(foodseq_merged$Timepoint_categorical, levels = lvls)

compare_groups <- list(c("Infant (M9)", "Infant (M12)"), c("Infant (M9)", "Mother"), c("Infant (M12)", "Mother"))


# Wilcoxon rank-sum tests with BH correction
raw_count_wilcox <- pairwise_wilcox_test(
  foodseq_merged,
  raw_counts ~ Timepoint_categorical,
  comparisons = compare_groups,
  p.adjust.method = "BH"
) %>%
  mutate(
    p.adj.label = paste0("adj. p = ", signif(p.adj)),
    y.position = c(143000, 160000, 177000)
  )


# ─────────────────────────────────────────────
# Plot raw sum counts for animal data
# ─────────────────────────────────────────────

fill_pal <- c("Infant (M9)"="#A6D8FF","Infant (M12)"="#F2C77C","Mother"="#66C2A3")
line_pal <- c("Infant (M9)"="#2E7DB2","Infant (M12)"="#B77700","Mother"="#006D5B")

p <- ggplot(
  foodseq_merged,
  aes(x = Timepoint_categorical, y = raw_counts)
) +
  geom_boxplot(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    width = 0.55, outlier.shape = NA, size = 0.7, fatten = 1.15, alpha = 0.9
  ) +
  geom_quasirandom(
    aes(fill = Timepoint_categorical, color = Timepoint_categorical),
    groupOnX = TRUE, width = 0.1, varwidth = FALSE,
    shape = 21, size = 1, stroke = 0.7, alpha = 0.7
  ) +
  stat_pvalue_manual(
    raw_count_wilcox,
    label = "p.adj.label",
    size = 3.6,
    tip.length = 0.01
  ) +
  scale_fill_manual(values = fill_pal, breaks = lvls) +
  scale_color_manual(values = line_pal,  breaks = lvls) +
  labs(
    x = "Sample group",
    y = "Raw Read Sum",
    title = "Animal"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EFEAE3", linewidth = 0.4),
    axis.line  = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(2.2, "mm"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 16, r = 14, b = 14, l = 10)
  )

p_out <- ggdraw(p) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = -0.5, vjust = 1.5,
                  fontface = "bold", size = 20)

p_out

ggsave(
  filename = "animal_foodseq_rawcounts.tiff",
  plot = p_out,
  dpi = 600, width = 4, height = 4.5, units = "in",
  compression = "lzw",
  bg = "white"
)
