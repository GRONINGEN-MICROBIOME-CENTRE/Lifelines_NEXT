## Code description
#######################################################################################################################################
## Script for Patching the Figure 2
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
library(ggtext)
library(patchwork)
library(ggsignif)
library(ggtree)
library(lmerTest)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
df_for_figure2c <- as.data.frame(read_tsv('df_for_figure2c.tsv'))
df_for_figure2d <- as.data.frame(read_tsv('df_for_figure2d.tsv'))
df_for_figure2e <- as.data.frame(read_tsv('df_for_figure2e.tsv'))
df_for_figure2e$Timepoint <- factor(df_for_figure2e$Timepoint, levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12"))
df_for_figure2e <- df_for_figure2e[df_for_figure2e$Dataset == "RPKM_count", ]

df_for_figure2d$Type <- factor(df_for_figure2d$Type, levels = c("same", "different"))
df_for_figure2d$Timepoint_numeric

#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\plots\\neg_ctrl_sharing")
#######################################################################################################################################

## Patching the Figure 2
#######################################################################################################################################
figure_2A <- readRDS(file="burkholderia.rds")
figure_2A$labels$tag <- "a"

figure_2B <- readRDS(file="phix.rds")
figure_2B$labels$tag <- "b"

figure_2C <- readRDS(file="micro_bacteroides.rds")
figure_2C$labels$tag <- "c"

figure_2D <- ggplot(df_for_figure2c, aes(x = type_cohort, y = Distance)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +  # Apply log scale to the y-axis
  #  coord_cartesian(ylim = c(new_min_y, 1)) +  # Set limits in the original scale
  facet_grid(cohort_nc ~ category, scales = "free", labeller = labeller(
    cohort_nc = c(
      "garmaeva" = "NCs Garmaeva *et al.* <br>",
      "liang" = "NCs Liang *et al.* <br>",
      "maqsood" = "NCs Maqsood *et al.* <br>",
      "shah" = "NCs Shah *et al.* <br>"
    ),
    category = c(
      "NCs" = "NCs and NCs",
      "Samples" = "NCs and Samples"
    )
  )) +
  labs(y = "log10(1 - (Bray-Curtis dissimilarity))", x = "Type of cohort \n (samples from the same or different cohort compared)", tag="d") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )

figure_2E <- ggplot(df_for_figure2d, aes(x=different_cohort_NC_presence, y=same_cohort_NC_presence)) +
  geom_point(size = 0.5, color="#2E236C", alpha=0.75) +  # Adjusted point size and added transparency
  geom_smooth(method=lm, color="#2E236C", fill="#C8ACD6", se=TRUE, linewidth=0.5) +  # Use linewidth instead of size
  labs(x = "% shared vOTUs with NCs \n from different studies", y = "% shared vOTUs with NCs \n from same study") +
  labs(tag="e") +
  facet_wrap(~ cohort, nrow = 2, ncol = 2, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )

figure_2F <- ggplot(df_for_figure2e, aes(x = Timepoint, y = value)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>"
    )
  )) +
  labs(y = "log10(% shared vOTUs)", tag="f") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )



# Combine the plots using patchwork
combined_plot2 <- (figure_2A + figure_2B + figure_2C + plot_layout(nrow=1, guides = "collect")) / (figure_2D + figure_2E + figure_2F) +
  plot_layout(heights = c(2.2, 1))# Save the combined plot as a PDF

ggsave("combined_figure2.png", combined_plot2, width = 45/2.54, height = 30/2.54)
#######################################################################################################################################