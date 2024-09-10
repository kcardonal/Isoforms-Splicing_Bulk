# Script to perform functional annotation of differentially expressed transcripts

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(ggplot2)

# Load the differentially expressed transcripts optimized
# Define the sheet names
sheet_names <- c("Common_up", "Common_down", "Common_up_to_down", "Common_down_to_up")  # Add all the sheet names you want to load
# Initialize an empty list to store the data frames
data_frames <- list()
# Column number to check for NA values
column_number <- 13  # Change this to the appropriate column number (gene_id)

# Loop through the sheet names and read each sheet into a data frame
for (sheet in sheet_names) {
  # Read the sheet into a data frame
  df <- read.xlsx("../results/upset_intersections_iPSCs.xlsx", sheet = sheet)
  # Remove rows with NA values in the specified column
  df <- df[!is.na(df[[column_number]]), ]
  # Store the data frame in the list
  data_frames[[sheet]] <- df
}

# Access the data frames using the sheet names
#IPSCs
det_ipscs_up <- data_frames[["Common_up"]]
det_ipscs_down <- data_frames[["Common_down"]]
det_ipscs_up_to_down <- data_frames[["Common_up_to_down"]]
det_ipscs_down_to_up <- data_frames[["Common_down_to_up"]]

#________________________________DEFINING BACKGROUND FOR ORA ANALYSIS______________________________________

#define the background gene set
#using lowly expressed genes as background
load("../Large_Files_No_repo/combat_seq_counts_isoforms.RData") # iso_corrected
# get the lowly expressed genes by calculating the mean of the counts of all samples
av_expr <- rowMeans(iso_corrected)
# set a threshold at the 10th percentile of the expression distribution
threshold <- quantile(av_expr, 0.2)
# Identify transcripts with low expression
lowly_expressed <- rownames(iso_corrected)[av_expr < threshold]
# convert the transcript ids to entrez ids
lowly_expressed_entrez <- mapIds(org.Hs.eg.db, keys = lowly_expressed, column = "ENTREZID", keytype = "ENSEMBLTRANS",multiVals = "first")
# Remove NA values
lowly_expressed_entrez <- lowly_expressed_entrez[!is.na(lowly_expressed_entrez)] # 40k transcripts mapped to 6577 genes

#________________________________FUNCTION TO PERFORM ENRICHMENT ANALYSIS______________________________________

perform_enrichment_analysis <- function(df, column_name, universe) {
  # Map the gene ids of the differentially expressed transcripts to entrez ids
  entrez_ids <- mapIds(org.Hs.eg.db, keys = df[[column_name]], column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  
  # Remove NA values
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Convert the entrez ids to use it in the enrichment analysis
  entrez_ids_df <- data.frame(ensembl_gene_id = names(entrez_ids), entrez_id = entrez_ids, stringsAsFactors = FALSE)
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(entrez_ids_df$entrez_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH", universe = universe, readable = TRUE)
  
  # Perform KEGG enrichment analysis
  kegg_enrichment <- enrichKEGG(entrez_ids_df$entrez_id, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  # Return the results as data frames
  return(list(go_results = as.data.frame(go_enrichment), kegg_results = as.data.frame(kegg_enrichment)))
}


results_up_ipscs <- perform_enrichment_analysis(det_ipscs_up, "gene_id_KS", lowly_expressed_entrez)
results_down_ipscs <- perform_enrichment_analysis(det_ipscs_down, "gene_id_KS", lowly_expressed_entrez)
results_up_to_down_ipscs <- perform_enrichment_analysis(det_ipscs_up_to_down, "gene_id_KS", lowly_expressed_entrez)
#results_down_to_up_ipscs <- perform_enrichment_analysis(det_ipscs_down_to_up, "gene_id_KS", lowly_expressed_entrez) #to few genes

# go_results <- results$go_results
# kegg_results <- results$kegg_results












#________________________________GO AND KEGG ENRICHMENT ANALYSIS______________________________________

# Map the gene ids of the differentially expressed transcripts to entrez ids
entrez_ids <- mapIds(org.Hs.eg.db, keys = det_ipscs_up$gene_id_KS, column = "ENTREZID", keytype = "ENSEMBL",multiVals = "first")
# Remove NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
# Convert the entrez ids to use it in the enrichment analysis
entrez_ids_df <- data.frame(ensembl_gene_id = names(entrez_ids), entrez_id = entrez_ids, stringsAsFactors = FALSE)

# Perform GO enrichment analysis
go_enrichment <- enrichGO(entrez_ids_df$entrez_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH", universe = lowly_expressed_entrez, readable = TRUE)
# visualise the results
dotplot(go_enrichment, showCategory = 10)
# Perform KEGG enrichment analysis
kegg_enrichment <- enrichKEGG(entrez_ids_df$entrez_id, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
# visualise the results
dotplot(kegg_enrichment, showCategory = 10)