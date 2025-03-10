##Script for diverse plots and visualization of the data

# Load libraries
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ggrepel)
library(dplyr)
library(forcats)
library(VennDiagram)
library(grid)
library(readr)

# Load data: DTE results
names_sheets <- c("KS_iPSCs", "HGA_iPSCs", "KS_NSCs", "HGA_NSCs", "KS_Neurons", "HGA_Neurons")

dte_results <- lapply(names_sheets, function(sheet) {
  read.xlsx("./results/DTE_results_annotated.xlsx", sheet = sheet)
})

# Asign names to each data frame
ks_ipscs_res <- dte_results[[1]]
hga_ipscs_res <- dte_results[[2]]
ks_nscs_res <- dte_results[[3]]
hga_nscs_res <- dte_results[[4]]
ks_neurons_res <- dte_results[[5]]
hga_neurons_res <- dte_results[[6]]

#Load data: DGE results
names_sheets <- c("iPSCs_47", "iPSCs_HGA", "NSCs_47", "NSCs_HGA", "Neurons_47", "Neurons_HGA")

dge_results <- lapply(names_sheets, function(sheet) {
  read.xlsx("../results_DGE/DGE_results_annotated.xlsx", sheet = sheet)
})

# Asign names to each data frame
ipscs_47_deg <- dge_results[[1]]
ipscs_hga_deg <- dge_results[[2]]
nscs_47_deg <- dge_results[[3]]
nscs_hga_deg <- dge_results[[4]]
neurons_47_deg <- dge_results[[5]]
neurons_hga_deg <- dge_results[[6]]

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

#PIE CHART TO SHOW THE PROPORTION OF XCI STATUS

create_xci_pie_chart <- function(data, file_name = "xci_pie_chart.png") {
  
  # Filter out transcripts categorized as "ns" and on chromosome X
  filtered_data <- data %>%
    filter(DE != "NS" & chromosome == "X")
  
  # Categorize transcripts by XCI status
  dte_categories <- filtered_data %>%
    mutate(category = case_when(
      XCI_status == "escapes" ~ "escapes",
      XCI_status == "inactive" ~ "inactive",
      XCI_status == "variable" ~ "variable",
      XCI_status == "PAR" ~ "PAR",
      is.na(XCI_status) ~ "NA",
      TRUE ~ "Other"
    )) %>%
    filter(category != "Other") %>%
    count(category, name = "value") %>%
    mutate(proportion = value / sum(value) * 100) # Calculate percentages
  
  # Calculate positions for the labels
  dte_categories <- dte_categories %>%
    mutate(csum = rev(cumsum(rev(proportion))), # Cumulative sum of proportions
           pos = proportion / 2 + lead(csum, 1), # Midpoint for each slice
           pos = if_else(is.na(pos), proportion / 2, pos)) # Handle NA for the last position
  
  # Create the pie chart
  plot <- ggplot(dte_categories, aes(x = "", y = proportion, fill = fct_inorder(category))) +
    geom_col(width = 1, color = 1) + # Pie chart slices
    geom_text(aes(label = paste0(round(proportion, 1), "%")),
              position = position_stack(vjust = 0.5)) + # Labels inside slices
    coord_polar(theta = "y") + # Convert to polar coordinates
    guides(fill = guide_legend(title = "Category")) + # Legend for categories
    scale_y_continuous(breaks = dte_categories$pos, labels = dte_categories$category) + # Align labels
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 12), 
          legend.position = "none", # Remove legend
          panel.background = element_rect(fill = "white")) # White background
  
  # Save the plot
  ggsave(filename = file_name, plot = plot, width = 6, height = 6, dpi = 300)
}

# ks ipscs
create_xci_pie_chart(res_KS_iPSCs_annotated, "./plots/ks_ipscs_xci_status.png")
# hga ipscs
create_xci_pie_chart(res_HGA_iPSCs_annotated, "./plots/hga_ipscs_xci_status.png")
# ks nscs
create_xci_pie_chart(res_KS_NSCs_annotated, "./plots/ks_nscs_xci_status.png")
# hga nscs
create_xci_pie_chart(res_HGA_NSCs_annotated, "./plots/hga_nscs_xci_status.png")
# ks neurons
create_xci_pie_chart(res_KS_Neurons_annotated, "./plots/ks_neurons_xci_status.png")
# hga neurons
create_xci_pie_chart(res_HGA_Neurons_annotated, "./plots/hga_neurons_xci_status.png")

#VENN DIAGRAM TO SHOW THE OVERLAP OF gDTE and DEGs

# Define the function
compare_gene_lists <- function(dte_data, dge_data, output_prefix) {
  
  # Extract unique gene IDs
  dte_genes <- unique(dte_data$gene_id)
  dge_genes <- unique(dge_data$ensembl_gene_id)

  # Find the intersection
  common_genes <- intersect(dte_genes, dge_genes)

  # Create a Venn diagram
  venn.plot <- draw.pairwise.venn(
    area1 = length(dte_genes),
    area2 = length(dge_genes),
    cross.area = length(common_genes),
    category = c("gDET", "DEG"),
    fill = c("blue", "yellow"),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.5,
    fontface = "bold"
  )

  # Save Venn diagram
  venn_file <- paste0(output_prefix, "_venn.png")
  png(venn_file)
  grid.draw(venn.plot)
  dev.off()

  # Save intersecting genes to a file
  intersect_file <- paste0(output_prefix, "_common_genes.tsv")
  write_tsv(data.frame(common_genes), intersect_file)

  # Print summary
  cat("Comparison for:", output_prefix, "\n")
  cat("Number of genes in DTE analysis:", length(dte_genes), "\n")
  cat("Number of genes in DGE analysis:", length(dge_genes), "\n")
  cat("Number of intersecting genes:", length(common_genes), "\n")
  cat("Venn diagram saved as:", venn_file, "\n")
  cat("Intersecting genes saved as:", intersect_file, "\n")
}

# Using it:
#ks ipscs
compare_gene_lists(ks_ipscs_res, ipscs_47_deg, "DTE_vs_DGE_ks_ipscs")
#hga ipscs
compare_gene_lists(hga_ipscs_res, ipscs_hga_deg, "DTE_vs_DGE_hga_ipscs")
#ks nscs
compare_gene_lists(ks_nscs_res, nscs_47_deg, "DTE_vs_DGE_ks_nscs")
#hga nscs
compare_gene_lists(hga_nscs_res, nscs_hga_deg, "DTE_vs_DGE_hga_nscs")
#ks neurons
compare_gene_lists(ks_neurons_res, neurons_47_deg, "DTE_vs_DGE_ks_neurons")
#hga neurons
compare_gene_lists(hga_neurons_res, neurons_hga_deg, "DTE_vs_DGE_hga_neurons")

