library(readr)
library(dplyr)
library(openxlsx)
library(stringr)
library(reshape2)
library(ggplot2)


meta_all_with_qc_curated <- as.data.frame(read_tsv("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\metadata\\metadata_with_qc_v1.tsv"))

meta_all_test_subset <- meta_all_with_qc_curated[!is.na(meta_all_with_qc_curated$Type) & !is.na(meta_all_with_qc_curated$cohort), c("Sample_name", "Type", "cohort", "Timepoint", "to_all_contigs",
                                                                                                                                    "to_1kb_contigs", "to_all_vir", "to_ext_vir", "to_ext_prun_vir")]

meta_all_test_subset_melt <- melt(meta_all_test_subset, id = c("Sample_name", "Type", "cohort", "Timepoint")) 

ggplot(meta_all_test_subset_melt, aes(x=Type, y=value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(x = "Sample type", y = "% total reads aligned") +
  facet_grid(variable ~ cohort) +
  scale_fill_manual(values = c("#ffd42f", "#138468", "#7849b8", "#f2609e", "#009dd6", "#ec111a")) +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 90),
        axis.text.y = element_text(size = 16))


ggplot(meta_all_test_subset_melt, aes(x=Type, y=value, fill = Type)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(x = "Sample type", y = "% total reads aligned") +
  facet_grid(Cohort ~ variable) +
  scale_fill_manual(values = c("#ffd42f", "#138468", "#7849b8", "#f2609e", "#009dd6", "#ec111a"))+
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 90),
        axis.text.y = element_text(size = 16))
