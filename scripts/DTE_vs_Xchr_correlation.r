# Correlation analysis between transcript expression levels and number of X chromosomes

# Load required libraries
library(ggplot2)
library(reshape2)
library(openxlsx)
library(dplyr)

# Load data
# Load the combat-seq corrected data - object iso_corrected
load("../Large_Files_No_repo/combat_seq_counts_isoforms_HTSFiltered.RData")
# Load the sample information, taking the sample_info object
load("../Large_Files_No_repo/isoforms_data.RData")

# save the iso_corrected object as a an excel file for easy access later 
#write.xlsx(as.data.frame(iso_corrected), "./results/normalized_corrected_counts.xlsx")

# Load the data and filter as needed
metadata <- sample_info
cell_type_of_interest <- "Neurons" # change as neeeded iPSCs, NSC, Neurons
metadata_filtered <- metadata %>% filter(Cell_Type == cell_type_of_interest)
# filter expression data to only include the cell type of interest
iso_corrected_filtered <- iso_corrected[, colnames(iso_corrected) %in% rownames(metadata_filtered)]
#check order
all(colnames(iso_corrected_filtered) == rownames(metadata_filtered))
# check dimensions
dim(iso_corrected_filtered)
dim(metadata_filtered)


# Extract the X chromosome counts
metadata_filtered$N_Xchr <- stringr::str_count(metadata_filtered$Karyotype, "X") # Count the number of X chromosomes

# Calculate the correlation between the number of X chromosomes and the expression of each gene
cor_results <- apply(iso_corrected_filtered, 1, function(gene_expr) {
  cor_test <- cor.test(gene_expr, metadata_filtered$N_Xchr, method = "spearman")  # Use "pearson" if linear
  c(correlation = cor_test$estimate, p_value = cor_test$p.value)
})

# Convert results to a data frame
cor_results_df <- as.data.frame(t(cor_results))

# Adjust p-values for multiple testing (FDR correction)
cor_results_df$adj_p_value <- p.adjust(cor_results_df$p_value, method = "fdr")

# Add gene names
cor_results_df$transcript <- rownames(iso_corrected_filtered)

# Step 3: Filter significant genes (FDR < 0.05)
significant_genes <- cor_results_df %>% filter(adj_p_value < 0.05)

# Sort by correlation coefficient (highest positive and negative correlations)
significant_genes <- significant_genes %>% arrange(desc(correlation.rho))

# Save results to a file
write.xlsx(significant_genes, file= "../results_DTE/Sperman_Corr_Xchr_vs_DETs_Neurons.xlsx")

# Step 4: Visualize an example gene
if (nrow(significant_genes) > 0) {
  gene_to_plot <- significant_genes$transcript[1]  # Choose the first significant transcript
  
  plot_df <- data.frame(
    X_chrom_count = metadata_filtered$N_Xchr,
    Expression = as.numeric(iso_corrected_filtered[gene_to_plot, ])
  )

  ggplot(plot_df, aes(x = X_chrom_count, y = Expression)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = paste("Expression vs Number of X in Neurons", gene_to_plot),
         x = "Number of X Chromosomes",
         y = "Normalized Expression") +
    theme_minimal()
}
# Save the plot
ggsave("./plots/Expression_vs_Number_of_Xchr_iPSCs.png", width = 6, height = 4, dpi = 300)