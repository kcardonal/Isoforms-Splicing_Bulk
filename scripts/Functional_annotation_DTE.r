# Script to perform functional annotation of differentially expressed transcripts

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(ggplot2)
library(openxlsx)


# Load the differentially expressed transcripts optimized
# Define the sheet names, change as needed
sheet_names <- c("Common_up", "Common_down", "Up_to_down", "down_to_up")  # Add all the sheet names you want to load
#sheet_names <- c("Common_up", "Common_down") 
# Initialize an empty list to store the data frames
data_frames <- list()
# Column number to check for NA values
column_number <- 25  # Change this to the appropriate column number (gene_id)
# File with the intersections,change as needed
#file <- "./results/upset_intersections_iPSCs_run2.xlsx"
#file <- "./results/upset_intersections_NSCs_run2.xlsx"
#file <- "./results/upset_intersections_Neurons_run2.xlsx"
#file <- "./results/common_intersect_KS_all.xlsx"
file <- "./results/common_intersect_HGA_all.xlsx"

# Loop through the sheet names and read each sheet into a data frame
for (sheet in sheet_names) {
  # Read the sheet into a data frame
  df <- read.xlsx(file, sheet = sheet)
  # Remove rows with NA values in the specified column
  df <- df[!is.na(df[[column_number]]), ]
  # Store the data frame in the list
  data_frames[[sheet]] <- df
}

# Access the data frames using the sheet names
#IPSCs
det_ipscs_up <- data_frames[["Common_up"]]
det_ipscs_down <- data_frames[["Common_down"]]
det_ipscs_up_to_down <- data_frames[["Up_to_down"]]
det_ipscs_down_to_up <- data_frames[["down_to_up"]]

#NSCs
det_nscs_up <- data_frames[["Common_up"]]
det_nscs_down <- data_frames[["Common_down"]]
det_nscs_up_to_down <- data_frames[["Up_to_down"]]
det_nscs_down_to_up <- data_frames[["down_to_up"]]

#Neurons
det_neurons_up <- data_frames[["Common_up"]]
det_neurons_down <- data_frames[["Common_down"]]
det_neurons_up_to_down <- data_frames[["Up_to_down"]]
det_neurons_down_to_up <- data_frames[["down_to_up"]]

#KS
det_ks_up <- data_frames[["Common_up"]]
det_ks_down <- data_frames[["Common_down"]]

#HGA
det_hga_up <- data_frames[["Common_up"]]
det_hga_down <- data_frames[["Common_down"]]

#________________________________DEFINING BACKGROUND FOR ORA ANALYSIS______________________________________

#define the background gene set
#using lowly expressed genes as background: It is usually better to start with the ComBat-Seq corrected counts 
#since they account for technical variation such as batch effects, which can particularly affect lowly expressed genes.
#The corrected counts should provide a more reliable basis for determining the expression distribution.

load("../Large_Files_No_repo/combat_seq_counts_isoforms_HTSFiltered.RData") # iso_corrected
# get the lowly expressed genes by calculating the mean of the counts of all samples
av_expr <- rowMeans(iso_corrected)
# set a threshold at the 10th percentile of the expression distribution
threshold <- quantile(av_expr, 0.25)
# Identify transcripts with low expression
lowly_expressed <- rownames(iso_corrected)[av_expr < threshold]
#Check how many transcripts are lowly expressed
length(lowly_expressed) # 4670 transcripts
# convert the transcript ids to entrez ids
lowly_expressed_entrez <- mapIds(org.Hs.eg.db, keys = lowly_expressed, column = "ENTREZID", keytype = "ENSEMBLTRANS",multiVals = "first")
# Remove NA values
lowly_expressed_entrez <- lowly_expressed_entrez[!is.na(lowly_expressed_entrez)] # 4.6k transcripts mapped to 388 genes

#________________________________FUNCTION TO PERFORM ENRICHMENT ANALYSIS______________________________________

perform_enrichment_analysis <- function(df, column_name, universe) {
  # Map the gene ids of the differentially expressed transcripts to entrez ids
  entrez_ids <- mapIds(org.Hs.eg.db, keys = df[[column_name]], column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  
  # Remove NA values
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Perform GO enrichment analysis for all ontologies
  go_enrichment <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH", universe = universe, readable = TRUE)
  
  # Perform KEGG enrichment analysis
  kegg_enrichment <- enrichKEGG(entrez_ids, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
  
  # Convert the results to data frames
  go_results_df <- as.data.frame(go_enrichment)
  kegg_results_df <- as.data.frame(kegg_enrichment)
  
  # Separate the GO results by ontology
  go_bp <- go_results_df[go_results_df$ONTOLOGY == "BP", ]
  go_mf <- go_results_df[go_results_df$ONTOLOGY == "MF", ]
  go_cc <- go_results_df[go_results_df$ONTOLOGY == "CC", ]
  
  # Return the results as a list of data frames and original enrichment objects
  return(list(go_bp = go_bp, go_mf = go_mf, go_cc = go_cc, kegg_results = kegg_results_df, go_all =go_results_df,go_enrichment = go_enrichment, kegg_enrichment = kegg_enrichment))
}

# Enrichment results for iPCSs
results_up_ipscs <- perform_enrichment_analysis(det_ipscs_up, "gene_id_KS", lowly_expressed_entrez)
results_down_ipscs <- perform_enrichment_analysis(det_ipscs_down, "gene_id_KS", lowly_expressed_entrez)
results_up_to_down_ipscs <- perform_enrichment_analysis(det_ipscs_up_to_down, "gene_id_KS", lowly_expressed_entrez)
results_down_to_up_ipscs <- perform_enrichment_analysis(det_ipscs_down_to_up, "gene_id_KS", lowly_expressed_entrez)

# Enrichment results for NSCs
results_up_nscs <- perform_enrichment_analysis(det_nscs_up, "gene_id_KS", lowly_expressed_entrez)
results_down_nscs <- perform_enrichment_analysis(det_nscs_down, "gene_id_KS", lowly_expressed_entrez)
results_up_to_down_nscs <- perform_enrichment_analysis(det_nscs_up_to_down, "gene_id_KS", lowly_expressed_entrez)
results_down_to_up_nscs <- perform_enrichment_analysis(det_nscs_down_to_up, "gene_id_KS", lowly_expressed_entrez)

#Enrichment results for Neurons
results_up_neurons <- perform_enrichment_analysis(det_neurons_up, "gene_id_KS", lowly_expressed_entrez)
results_down_neurons <- perform_enrichment_analysis(det_neurons_down, "gene_id_KS", lowly_expressed_entrez)
results_up_to_down_neurons <- perform_enrichment_analysis(det_neurons_up_to_down, "gene_id_KS", lowly_expressed_entrez)
results_down_to_up_neurons <- perform_enrichment_analysis(det_neurons_down_to_up, "gene_id_KS", lowly_expressed_entrez)

#Enrichment results for KS
results_up_ks <- perform_enrichment_analysis(det_ks_up, "gene_id", lowly_expressed_entrez)
results_down_ks <- perform_enrichment_analysis(det_ks_down, "gene_id", lowly_expressed_entrez)

#Enrichment results for HGA
results_up_hga <- perform_enrichment_analysis(det_hga_up, "gene_id", lowly_expressed_entrez)
results_down_hga <- perform_enrichment_analysis(det_hga_down, "gene_id", lowly_expressed_entrez)


# Save the results of those that produced significant enrichment in iPCSs
write.xlsx(results_up_ipscs$kegg_results, file = "../results/KEGG_up_iPSCs_v2.xlsx")
write.xlsx(results_up_to_down_ipscs$kegg_results, file = "../results/KEGG_up_to_down_iPSCs_v2.xlsx")
# Plot the results from the KEGG enrichment analysis in the significant in iPCSs
dotplot(results_up_ipscs$kegg_enrichment, showCategory = 10)
# Save the plot
ggsave("./plots/KEGG_up_iPSCs_v2.png")
dotplot(results_up_to_down_ipscs$kegg_enrichment, showCategory = 10)
ggsave("./plots/KEGG_up_to_down_iPSCs_v2.png")


# Save the results of those that produced significant enrichment in NSCs
write.xlsx(results_up_nscs$kegg_results, file = "../results/KEGG_common_up_NSCs_v2.xlsx")
write.xlsx(results_up_to_down_nscs$kegg_results, file = "../results/KEGG_common_up_to_down_NSCs_v2.xlsx")
# Plot the results from the KEGG enrichment analysis in the significant in NSCs
dotplot(results_up_nscs$kegg_enrichment, showCategory = 10)
dotplot(results_up_to_down_nscs$kegg_enrichment, showCategory = 10)
# Save the plot
ggsave("./plots/KEGG_common_up_NSCs_v2.png")
ggsave("./plots/KEGG_common_up_to_down_NSCs_v2.png")

# Neurons
write.xlsx(results_up_neurons$kegg_results, file = "../results/KEGG_common_up_Neurons.xlsx")
write.xlsx(results_up_to_down_neurons$kegg_results, file = "../results/KEGG_common_up_to_down_Neurons.xlsx")
#plot the results from the KEGG enrichment analysis in the significant in Neurons
dotplot(results_up_neurons$kegg_enrichment, showCategory = 10)
ggsave("./plots/KEGG_common_up_Neurons_v2.png")
dotplot(results_up_to_down_neurons$kegg_enrichment, showCategory = 10)
ggsave("./plots/KEGG_common_up_to_down_Neurons_v2.png")


# Save the results of those that produced significant enrichment in KS
write.xlsx(results_up_ks$kegg_results, file = "../results/KEGG_common_up_KS.xlsx")
#plot the results from the KEGG enrichment analysis in the significant in KS
dotplot(results_up_ks$kegg_enrichment, showCategory = 10)
#save the plot
ggsave("./plots/KEGG_common_up_KS.png")

# Save the results of those that produced significant enrichment in HGA
write.xlsx(results_up_hga$kegg_results, file = "../results/KEGG_common_up_HGA.xlsx")