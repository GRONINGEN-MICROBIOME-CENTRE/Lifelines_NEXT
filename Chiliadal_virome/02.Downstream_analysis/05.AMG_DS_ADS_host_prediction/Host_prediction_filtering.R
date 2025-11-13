## Load libraries
library(tidyverse)

## Load data
HP <- read.delim("../data/Host_prediction_to_genus_m90.csv", sep = ",", header=TRUE)

# Filtering assignments for general analysis
HP$Host.genus.short <- sub(".*;g__", "", HP$Host.genus)

HP_filtered <- HP %>%
  filter(Host.genus.short != "" & Host.genus != "Unknown") %>%
  group_by(Virus) %>%
  filter(Confidence.score == max(Confidence.score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Virus, Confidence.score, Host.genus.short, .keep_all = TRUE)

# Filtering the ones that are with the multiple hits as per Confidence.score 
HP_filtered <- HP_filtered %>%
  mutate(row_id = row_number()) %>%  # Add row ID to keep track of rows
  separate_rows(List.of.methods, sep = " ") %>%  # Split into method;score
  separate(List.of.methods, into = c("method", "score"), sep = ";") %>%  # Split method and score
  pivot_wider(names_from = method, values_from = score) %>%  # Reshape to wide format
  select(-row_id) 

# For unified filtering
set.seed(13)

HP_filtered <- HP_filtered %>%
  mutate(across(c(CRISPR, `iPHoP-RF`, RaFAH, blast), 
                ~ replace_na(as.numeric(.x), 0))) %>%
  group_by(Virus) %>%
  slice_max(order_by = blast, with_ties = TRUE) %>%
  slice_max(order_by = CRISPR, with_ties = TRUE) %>%
  slice_max(order_by = `iPHoP-RF`, with_ties = TRUE) %>%
  slice_max(order_by = RaFAH, with_ties = TRUE) %>%
  slice_sample(n = 1) %>%
  ungroup()

write.table(HP_filtered, "../results/Host_prediction_to_genus_filtered_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)  # final table output

# Basic plot to see what's in there
for_plot <- as.data.frame(table(HP_filtered$Host.genus))
for_plot$genus <- sub(".*;g__", "", for_plot$Var1)
for_plot$class <- sub(".*;c__", "", for_plot$Var1)
for_plot$class <- sub(";o__.*$", "", for_plot$class)

top <- for_plot %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 75)

ggplot(top, aes(x = reorder(genus, -Freq), y = Freq, fill = class)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Genus",
    y = "Frequency",
    title = "Top 75 genera colored by class"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_line()
  ) +
  scale_fill_brewer(palette = "Set3")  # Choose any palette you like