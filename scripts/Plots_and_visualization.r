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

#PIE CHART TO SHOW THE PROPORTION OF DETs BY CHROMOSOME CATEGORY
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


#BAR PLOT TO SHOW THE NUMBER OF DETs PER CHROMOSOME
create_dte_bar_plot <- function(data, file_name = "dte_bar_plot.png") {

  # Filter out transcripts categorized as "ns" and count DETs per chromosome
  chromosome_counts <- data %>%
    filter(DE != "NS") %>%
    count(chromosome, name = "count") %>%
    arrange(desc(count))  # Order by count for a more informative plot
  
  # Create bar plot
  plot <- ggplot(chromosome_counts, aes(x = reorder(chromosome, -count), y = count, fill = chromosome)) +
    geom_bar(stat = "identity") +
    labs(title = "Number of DETs per Chromosome",
         x = "Chromosome",
         y = "Number of DETs") +
    theme_minimal() +
    theme(legend.position = "none")  # Hide legend as each bar is a chromosome
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 8, height = 5, dpi = 300)
}

# Usage example
create_dte_bar_plot(ks_ipscs_res, "./Isoforms-Splicing_Bulk/plots/ks_ipscs_bar_plot_dte.png")
create_dte_bar_plot(hga_ipscs_res, "./Isoforms-Splicing_Bulk/plots/hga_ipscs_bar_plot_dte.png")
create_dte_bar_plot(ks_nscs_res, "./Isoforms-Splicing_Bulk/plots/ks_nscs_bar_plot_dte.png")
create_dte_bar_plot(hga_nscs_res, "./Isoforms-Splicing_Bulk/plots/hga_nscs_bar_plot_dte.png")
create_dte_bar_plot(ks_neurons_res, "./Isoforms-Splicing_Bulk/plots/ks_neurons_bar_plot_dte.png")
create_dte_bar_plot(hga_neurons_res, "./Isoforms-Splicing_Bulk/plots/hga_neurons_bar_plot_dte.png")

#BAR PLOT OF TRANSCRIPT BIOTYPES

plot_biotype_distribution <- function(data, file_name = "biotype_distribution.png") {
  library(dplyr)
  library(ggplot2)
  
  # Filter out non-significant transcripts
  filtered_data <- data %>%
    filter(DE != "NS")
  
  # Count DETs by biotype and group low counts as "Others"
  biotype_counts <- filtered_data %>%
    count(tx_biotype) %>%
    mutate(tx_biotype = ifelse(n < 11, "Others", tx_biotype)) %>%
    group_by(tx_biotype) %>%
    summarize(count = sum(n)) %>%
    arrange(desc(count))
  
  # Bar plot for transcript biotypes
  plot <- ggplot(biotype_counts, aes(x = reorder(tx_biotype, count), y = count, fill = tx_biotype)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Distribution of DETs by Transcript Biotype",
         x = "Transcript Biotype",
         y = "Number of DETs") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 8, height = 6, dpi = 300)
}

# Usage example
plot_biotype_distribution(ks_ipscs_res, "./Isoforms-Splicing_Bulk/plots/biotype_distr_ks_ipscs.png")
plot_biotype_distribution(hga_ipscs_res, "./Isoforms-Splicing_Bulk/plots/biotype_distr_hga_ipscs.png")
plot_biotype_distribution(ks_nscs_res, "./Isoforms-Splicing_Bulk/plots/biotype_distr_ks_nscs.png")
plot_biotype_distribution(hga_nscs_res, "./Isoforms-Splicing_Bulk/plots/biotype_distr_hga_nscs.png")
plot_biotype_distribution(ks_neurons_res, "./Isoforms-Splicing_Bulk/plots/biotype_distr_ks_neurons.png")
plot_biotype_distribution(hga_neurons_res, "./Isoforms-Splicing_Bulk/plots/biotype_distr_hga_neurons.png")

#DENSITY PLOT OF THE LOG2FC DISTRIBUTION
plot_log2FC_distribution <- function(data, file_name = "log2FC_distribution.png") {

  # Filter out non-significant transcripts
  filtered_data <- data %>%
    filter(DE != "NS")
  
  # Density plot for log2FoldChange
  plot <- ggplot(filtered_data, aes(x = log2FoldChange, fill = DE)) +
    geom_density(alpha = 0.6) +
    labs(title = "Distribution of Log2 Fold Changes",
         x = "Log2 Fold Change",
         y = "Density",
         fill = "Regulation") +
    theme_minimal()
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 8, height = 5, dpi = 300)
}

# Usage example
plot_log2FC_distribution(ks_ipscs_res, "log2FC_distribution.png")

#BOXPLOT OF THE LOG2FC BY CHROMOSOME

plot_log2FC_by_chromosome <- function(data, file_name = "log2FC_by_chromosome.png") {

  # Filter out non-significant transcripts
  filtered_data <- data %>%
    filter(DE != "NS")
  
  # Define the chromosome order
  chromosome_levels <- c(as.character(1:22), "X", "Y")
  
  # Set chromosome as a factor with the specified levels
  filtered_data$chromosome <- factor(filtered_data$chromosome, levels = chromosome_levels)
  
  # Boxplot for log2FoldChange by chromosome
  plot <- ggplot(filtered_data, aes(x = chromosome, y = log2FoldChange)) +
    geom_boxplot(outlier.alpha = 0.2) +
    labs(title = "Log2 Fold Change by Chromosome",
         x = "Chromosome",
         y = "Log2FC") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 10, height = 6, dpi = 300)
}

# Usage example
plot_log2FC_by_chromosome(ks_ipscs_res, "./Isoforms-Splicing_Bulk/plots/log2FC_density_ks_ipscs.png")
plot_log2FC_by_chromosome(hga_ipscs_res, "./Isoforms-Splicing_Bulk/plots/log2FC_density_hga_ipscs.png")
plot_log2FC_by_chromosome(ks_nscs_res, "./Isoforms-Splicing_Bulk/plots/log2FC_density_ks_nscs.png")
plot_log2FC_by_chromosome(hga_nscs_res, "./Isoforms-Splicing_Bulk/plots/log2FC_density_hga_nscs.png")
plot_log2FC_by_chromosome(ks_neurons_res, "./Isoforms-Splicing_Bulk/plots/log2FC_density_ks_neurons.png")
plot_log2FC_by_chromosome(hga_neurons_res, "./Isoforms-Splicing_Bulk/plots/log2FC_density_hga_neurons.png")