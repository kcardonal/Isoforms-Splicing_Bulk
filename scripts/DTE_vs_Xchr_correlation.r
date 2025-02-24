# Correlation analysis between transcript expression levels and number of X chromosomes

# Load required libraries
library(ggplot2)
library(reshape2)
library(openxlsx)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)


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

#__________________________PLOT ACP1 TRANSCRIPTS_____________________________________

# Step 4: Visualize all transcripts of ACP1
if (nrow(significant_genes) > 0) {
  # Filter for the transcripts of ACP1 manually using their transcript IDs
  gene_transcripts <- significant_genes %>% 
    filter(transcript %in% c("ENST00000272067", "ENST00000453390", "ENST00000272065"))  # Replace with ACP1 transcripts
  
  if (nrow(gene_transcripts) > 0) {
    # Prepare a data frame for all transcripts
    plot_df <- do.call(rbind, lapply(gene_transcripts$transcript, function(transcript) {
      data.frame(
        transcript_id = transcript,
        X_chrom_count = metadata_filtered$N_Xchr,  # Replace N_Xchr with your column for X chromosome counts
        Expression = as.numeric(iso_corrected_filtered[transcript, ])
      )
    }))
    
    # Create the multi-panel plot
    p <- ggplot(plot_df, aes(x = X_chrom_count, y = Expression)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      facet_wrap(~ transcript_id, scales = "free_y") +
      labs(title = "Expression vs Number of X Chromosomes for ACP1 Transcripts in Neurons",
           x = "Number of X Chromosomes",
           y = "Normalized Expression") +
      theme_minimal()
    
    # Save the plot
    ggsave(filename = "ACP1_MultiPanel_Plot.png", plot = p, width = 12, height = 8, dpi = 300)
  } else {
    print("No transcripts found for ACP1 in significant_genes.")
  }
} else {
  print("No significant genes found.")
}


#__________________________________________________________________________________________
##__________ Annotate the transcripts that correlates with the x dossage, with the gene information___________

# Load the correlation results
corr_ipscs <- read.xlsx("../results_DTE/Sperman_Corr_XchrDossage_vs_TE.xlsx" , sheet = "iPSCs")
corr_nscs <- read.xlsx("../results_DTE/Sperman_Corr_XchrDossage_vs_TE.xlsx" , sheet = "NSCs")
corr_neurons <- read.xlsx("../results_DTE/Sperman_Corr_XchrDossage_vs_TE.xlsx" , sheet = "Neurons")

#Make transcript_id column the rownames
rownames(corr_ipscs) <- corr_ipscs$transcript_id
rownames(corr_nscs) <- corr_nscs$transcript_id
rownames(corr_neurons) <- corr_neurons$transcript_id

#Remove the transcript_id column
corr_ipscs <- corr_ipscs[, -4]
corr_nscs <- corr_nscs[, -4]
corr_neurons <- corr_neurons[, -4]


# Get the transcript IDs from the results
transcript_ids_ipscs <- rownames(corr_ipscs)
transcript_ids_nscs <- rownames(corr_nscs)
transcript_ids_neurons <- rownames(corr_neurons)

# Connect to the Ensembl database
edb <- EnsDb.Hsapiens.v86

#for iPSCs

# Fetch data using ensembldb
annot_ipscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_ipscs), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
# Convert to data frame the annotation
annot_ipscs_ensembldb_df <- as.data.frame(annot_ipscs_ensembldb)
# Convert to data frame the results
corr_ipscs_df <- as.data.frame(corr_ipscs)
# Merge annotations with results
corr_ipscs_annotated <- merge(corr_ipscs_df, annot_ipscs_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
# Rename the row.names column to transcript_id for clarity
colnames(corr_ipscs_annotated)[1] <- "transcript_id"
# Rename the seqnames column to chromosome for clarity
colnames(corr_ipscs_annotated)[5] <- "chromosome"

# For NSCs
# Fetch data using ensembldb 
annot_nscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_nscs), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
# Convert to data frame the annotation
annot_nscs_ensembldb_df <- as.data.frame(annot_nscs_ensembldb)
# Convert to data frame the results
corr_nscs_df <- as.data.frame(corr_nscs)
# Merge annotations with results
corr_nscs_annotated <- merge(corr_nscs_df, annot_nscs_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
# Rename the row.names column to transcript_id for clarity
colnames(corr_nscs_annotated)[1] <- "transcript_id"
# Rename the seqnames column to chromosome for clarity
colnames(corr_nscs_annotated)[5] <- "chromosome"

# For Neurons
# Fetch data using ensembldb
annot_neurons_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_neurons), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
# Convert to data frame the annotation
annot_neurons_ensembldb_df <- as.data.frame(annot_neurons_ensembldb)
# Convert to data frame the results
corr_neurons_df <- as.data.frame(corr_neurons)
# Merge annotations with results
corr_neurons_annotated <- merge(corr_neurons_df, annot_neurons_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
# Rename the row.names column to transcript_id for clarity
colnames(corr_neurons_annotated)[1] <- "transcript_id"
# Rename the seqnames column to chromosome for clarity
colnames(corr_neurons_annotated)[5] <- "chromosome"

# Complete missing annotations with BiomaRt
#load biomart data. Modify each time as needed
biomart_data <- read.xlsx("../Large_Files_No_repo/MissingAnnot_corr.xlsx", sheet = "Neurons")
dim(biomart_data)

#for ipscs
# identify Rows with missing annotations
missing_annot <- corr_neurons_annotated[is.na(corr_neurons_annotated$chromosome), ]
dim(missing_annot)
# Merge the missing annotations with the biomart data
annotated_missing <- merge(missing_annot, biomart_data, by = "transcript_id", all.x = TRUE)
# Remove redundant columns - x
annotated_missing <- annotated_missing[, -c(5:13)]
# Remove the .y suffix from the columns but keep the original names
colnames(annotated_missing) <- gsub(".y", "", colnames(annotated_missing))
# change name of columns damaged by the .y subs
colnames(annotated_missing)[11] <- "gene_biotype"
colnames(annotated_missing)[12] <- "tx_biotype"
# Order columns as in the corr_ipscs_annotated data frame
annotated_missing <- annotated_missing[, colnames(corr_ipscs_annotated)]
# Save the annotated missing data
write.xlsx(annotated_missing, file = "../Large_Files_No_repo/MissingAnnot_neurons.xlsx")


#_________________Functional Annnotation of Genes with dossaage effect_________________

# Function to perform ORA analysis on the correlated tanscripts

# Define the function
perform_ora_analysis <- function(input_file, sheet_name, output_filename) {
  
  # Load correlation results
  corr_data <- read.xlsx(input_file, sheet = sheet_name)

  # Separate positive and negative correlation groups, modify correlation.rho to the needed threshold: abs(0.5), abs(0.6) and abs(0.7)
  positive_corr <- corr_data %>% filter(correlation.rho > 0.7)
  negative_corr <- corr_data %>% filter(correlation.rho < 0)

  # Extract unique gene IDs
  positive_genes <- unique(positive_corr$gene_id)
  negative_genes <- unique(negative_corr$gene_id)

  # Convert Ensembl Gene IDs to Entrez IDs
  positive_entrez <- na.omit(mapIds(org.Hs.eg.db, keys = positive_genes, 
                                    column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  negative_entrez <- na.omit(mapIds(org.Hs.eg.db, keys = negative_genes, 
                                    column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))

  # Define function for GO analysis (for BP, MF, CC)
  run_go_enrichment <- function(gene_list, ontology) {
    enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ontology, 
             pAdjustMethod = "fdr", qvalueCutoff = 0.05)
  }

  # Perform GO Enrichment Analysis
  go_positive_BP <- run_go_enrichment(positive_entrez, "BP")
  go_positive_MF <- run_go_enrichment(positive_entrez, "MF")
  go_positive_CC <- run_go_enrichment(positive_entrez, "CC")

  go_negative_BP <- run_go_enrichment(negative_entrez, "BP")
  go_negative_MF <- run_go_enrichment(negative_entrez, "MF")
  go_negative_CC <- run_go_enrichment(negative_entrez, "CC")

  # Perform KEGG Pathway Analysis
  kegg_positive <- enrichKEGG(gene = positive_entrez, organism = "hsa", keyType = "ncbi-geneid",
                              pAdjustMethod = "fdr", qvalueCutoff = 0.05)

  kegg_negative <- enrichKEGG(gene = negative_entrez, organism = "hsa", keyType = "ncbi-geneid",
                              pAdjustMethod = "fdr", qvalueCutoff = 0.05)

  # Save results into an Excel file with multiple sheets
  wb <- createWorkbook()

  # Add each result as a separate sheet
  addWorksheet(wb, "GO_BP_Positive")
  addWorksheet(wb, "GO_MF_Positive")
  addWorksheet(wb, "GO_CC_Positive")
  addWorksheet(wb, "GO_BP_Negative")
  addWorksheet(wb, "GO_MF_Negative")
  addWorksheet(wb, "GO_CC_Negative")
  addWorksheet(wb, "KEGG_Positive")
  addWorksheet(wb, "KEGG_Negative")

  # Write data to sheets
  writeData(wb, "GO_BP_Positive", as.data.frame(go_positive_BP))
  writeData(wb, "GO_MF_Positive", as.data.frame(go_positive_MF))
  writeData(wb, "GO_CC_Positive", as.data.frame(go_positive_CC))
  writeData(wb, "GO_BP_Negative", as.data.frame(go_negative_BP))
  writeData(wb, "GO_MF_Negative", as.data.frame(go_negative_MF))
  writeData(wb, "GO_CC_Negative", as.data.frame(go_negative_CC))
  writeData(wb, "KEGG_Positive", as.data.frame(kegg_positive))
  writeData(wb, "KEGG_Negative", as.data.frame(kegg_negative))

  # Save the Excel file
  saveWorkbook(wb, file = paste0(output_filename, ".xlsx"), overwrite = TRUE)

  # Visualization (optional)
  # Define output directory
  output_dir <- "./plots/ORA_Correlations_Plots_0.7"

  # Save dot plots for GO enrichment
  ggsave(filename = file.path(output_dir, "GO_BP_Positive.png"), 
       plot = dotplot(go_positive_BP, showCategory = 5, title = "GO BP - Positive Correlation"), 
       width = 8, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, "GO_BP_Negative.png"), 
       plot = dotplot(go_negative_BP, showCategory = 5, title = "GO BP - Negative Correlation"), 
       width = 8, height = 6, dpi = 300)
  
  ggsave(filename = file.path(output_dir, "GO_MF_Positive.png"), 
       plot = dotplot(go_positive_MF, showCategory = 5, title = "GO MF - Positive Correlation"), 
       width = 8, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, "GO_MF_Negative.png"), 
       plot = dotplot(go_negative_MF, showCategory = 5, title = "GO MF - Negative Correlation"), 
       width = 8, height = 6, dpi = 300)
  
  ggsave(filename = file.path(output_dir, "GO_CC_Positive.png"), 
       plot = dotplot(go_positive_CC, showCategory = 5, title = "GO CC - Positive Correlation"), 
       width = 8, height = 6, dpi = 300)
  
  ggsave(filename = file.path(output_dir, "GO_CC_Negative.png"), 
       plot = dotplot(go_negative_CC, showCategory = 5, title = "GO CC - Negative Correlation"), 
       width = 8, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, "KEGG_Positive.png"), 
       plot = dotplot(kegg_positive, showCategory = 5, title = "KEGG - Positive Correlation"), 
       width = 8, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, "KEGG_Negative.png"), 
       plot = dotplot(kegg_negative, showCategory = 5, title = "KEGG - Negative Correlation"), 
       width = 8, height = 6, dpi = 300)

}

# Run the function for ipscs
perform_ora_analysis("./results/Sperman_Corr_XchrDossage_vs_TE_annotated.xlsx", "iPSCs", "iPSCs_Corr_ORA_0.7")
# Run the function for nscs
perform_ora_analysis("./results/Sperman_Corr_XchrDossage_vs_TE_annotated.xlsx", "NSCs", "NSCs_Corr_ORA_0.7")
# Run the function for neurons
perform_ora_analysis("./results/Sperman_Corr_XchrDossage_vs_TE_annotated.xlsx", "Neurons", "Neurons_Corr_ORA_0.7")


#______________________Plot correlation of specific genes of interest________________________________
# Load the correlation results 
corr_data <- read.xlsx("./results/Sperman_Corr_XchrDossage_vs_TE_annotated.xlsx", sheet = "Neurons")

# Specify the gene of interest
gene_of_interest <- "ACP1"  # Replace with your desired gene_name

# Filter the correlation data for the gene of interest
gene_transcripts <- corr_data %>% filter(gene_name == gene_of_interest)

# Check if there are transcripts for the specified gene
if (nrow(gene_transcripts) > 0) {
  
  # Create a data frame for all transcripts belonging to the gene
  plot_data <- gene_transcripts %>%
    rowwise() %>%
    mutate(
      X_chrom_count = list(corr_data$X_chrom_count),  # Ensure this column exists
      Expression = list(corr_data[[transcript_id]])   # Add expression for each transcript
    ) %>%
    unnest(cols = c(X_chrom_count, Expression))
  
  # Generate a multi-panel plot using facet_wrap
  p <- ggplot(plot_data, aes(x = X_chrom_count, y = Expression)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    facet_wrap(~ transcript_id, scales = "free_y") +
    labs(title = paste("Correlation for Transcripts of", gene_of_interest),
         x = "Number of X Chromosomes",
         y = "Expression Level") +
    theme_minimal()
  
  # Save the multi-panel plot
  ggsave(filename = paste0("Correlation_", gene_of_interest, "_MultiPanel.png"), 
         plot = p, width = 12, height = 8, dpi = 300)
  
} else {
  print("No transcripts found for the specified gene.")
}