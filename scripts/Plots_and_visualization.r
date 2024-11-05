##Script for diverse plots and visualization of the data

# Load libraries
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ggrepel)
library(dplyr)

# Load data: DTE results
names_sheets <- c("KS_iPSCs", "HGA_iPSCs", "KS_NSCs", "HGA_NSCs", "KS_Neurons", "HGA_Neurons")

dte_results <- lapply(names_sheets, function(sheet) {
  read.xlsx("./results/DTE_results_2run.xlsx", sheet = sheet)
})

# Asign names to each data frame
ks_ipscs_res <- dte_results[[1]]
hga_ipscs_res <- dte_results[[2]]
ks_nscs_res <- dte_results[[3]]
hga_nscs_res <- dte_results[[4]]
ks_neurons_res <- dte_results[[5]]
hga_neurons_res <- dte_results[[6]]

# Plotting diverse statistics

#Pie chart yo show the proportion of DETs across chromosomes
create_dte_pie_chart <- function(data, file_name = "dte_pie_chart.png") {

  # Filter out transcripts categorized as "ns"
  filtered_data <- data %>%
    filter(DE != "NS")
  
  # Categorize transcripts by chromosome type
  dte_categories <- filtered_data %>%
    mutate(category = case_when(
      chromosome %in% as.character(1:22) ~ "Autosomes",
      chromosome == "X" ~ "Chr.X",
      chromosome == "Y" ~ "Chr.Y",
      TRUE ~ "Other"
    )) %>%
    filter(category != "Other") %>%
    count(category) %>%
    mutate(proportion = n / sum(n),
           label = paste0(category, " (", scales::percent(proportion, accuracy = 0.1), ")"))
  
  # Create pie chart with labels outside using geom_label_repel
  plot <- ggplot(dte_categories, aes(x = "", y = proportion, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    labs(title = "Proportion of DETs by Chromosome Category",
         x = NULL, y = NULL, fill = "Category") +
    theme_void() +
    geom_label_repel(aes(y = proportion / 2 + c(0, cumsum(proportion)[-length(proportion)]),
                         label = label),
                     nudge_x = 1,
                     show.legend = FALSE,
                     size = 4,
                     box.padding = 0.3,
                     point.padding = 0.3,
                     segment.size = 0.2)
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 6, height = 6, dpi = 300)
}

# ks ipscs
create_dte_pie_chart(ks_ipscs_res, "./Isoforms-Splicing_Bulk/plots/ks_ipscs_det_per_chr.png")
# hga ipscs
create_dte_pie_chart(hga_ipscs_res, "./Isoforms-Splicing_Bulk/plots/hga_ipscs_det_per_chr.png")
# ks nscs
create_dte_pie_chart(ks_nscs_res, "./Isoforms-Splicing_Bulk/plots/ks_nscs_det_per_chr.png")
# hga nscs
create_dte_pie_chart(hga_nscs_res, "./Isoforms-Splicing_Bulk/plots/hga_nscs_det_per_chr.png")
# ks neurons
create_dte_pie_chart(ks_neurons_res, "./Isoforms-Splicing_Bulk/plots/ks_neurons_det_per_chr.png")
# hga neurons
create_dte_pie_chart(hga_neurons_res, "./Isoforms-Splicing_Bulk/plots/hga_neurons_det_per_chr.png")

