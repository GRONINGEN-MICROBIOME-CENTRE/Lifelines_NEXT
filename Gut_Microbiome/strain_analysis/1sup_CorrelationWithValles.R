library(tidyverse)

Threshold = "Results/Thresholds.tsv"

if (!file.exists(Threshold)){
	file_list <- list.files("Distance_tables", pattern = "*threshold.tsv", full.names = TRUE)
	# Initialize an empty data frame to store the merged data
	merged_data <- data.frame()
	# Loop through each file
	for (file in file_list) {
	  # Read the file using read_tsv
	  data <- read_tsv(file)
	  # Merge the data with the existing data frame
	  merged_data <- rbind(merged_data, data)
	}
	write_tsv(merged_data, Threshold)
} else { read_tsv(Threshold) -> merged_data }


print("Comparison with Valles-Colomer threhsolds on Jan21")
merged_data %>% drop_na() -> merged_data
cor.test(merged_data$final_threshold, merged_data$Threshold_valles) %>% print()
merged_data %>% ggplot(aes(x=final_threshold, y=Threshold_valles)) + geom_point(size=2) + theme_bw() + 
geom_smooth(method = "lm", se = FALSE, color = "blue") + labs(x = "Sinha_Jan21_thresholds", y = "VallesColomer_Jan21_thresholds") -> Plot
ggsave("Results/ComparisonThresholds.pdf", Plot)



