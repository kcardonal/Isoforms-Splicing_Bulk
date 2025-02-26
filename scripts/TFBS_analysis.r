# Preparing data for TFBS analysis


# Load the required libraries
library(openxlsx)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# ________________ Setup of the query list for the TFBS analysis________________

# Load the correlation results to extract the gene list, change sheet name as needed
corr_res <- read.xlsx("./results/Sperman_Corr_XchrDossage_vs_TE_annotated.xlsx", sheet = "Neurons")

# Extract the unique gene list from the positive correlation results
gene_list_pos <- unique(corr_res$gene_id[corr_res$correlation.rho > 0])
length(gene_list_pos)
# ipscs > full:159, 0.5:7 0.6:3 0.7:2
# nscs > full:2427 0.5:285 0.6:71 0.7:17
# neurons > full: 141 0.5:141 0.6:101 0.7:29

# Extract the unique gene list from the negative correlation results
gene_list_neg <- unique(corr_res$gene_id[corr_res$correlation.rho < 0])
length(gene_list_neg)
# ipscs > full: 76 -0.5:2 -0.6:0 0.7:0
# nscs > full:1944 -0.5:56 -0.6:6 0.7:0
# neurons > full:176 0 -0.5:176 -0.6:116 0.7:17

#Save gene lists to a text file
write.table(gene_list_pos, file = "../results_DTE/forTFBS/gene_list_pos_neurons.txt", row.names = FALSE, col.names = FALSE)
write.table(gene_list_neg, file = "../results_DTE/forTFBS/gene_list_neg_neurons.txt", row.names = FALSE, col.names = FALSE)


#________________ Setup of the background list for the TFBS analysis________________

# Load the DTE analysis results 
#dds_combat_seq <- readRDS("../Large_Files_No_repo/DET_analysis_results.rds")

# Extract results for each condition and cell type

# For KS vs Control in iPSCs
res_KS_iPSCs <- results(dds_combat_seq, name = "Condition_KS_vs_Control", alpha=0.05)

# For HGA vs Control in iPSCs
res_HGA_iPSCs <- results(dds_combat_seq, name = "Condition_HGA_vs_Control", alpha=0.05)

# For KS vs Control in NSCs
res_KS_NSCs <- results(dds_combat_seq, contrast = list(c("Condition_KS_vs_Control","Cell_TypeNSC.ConditionKS")), alpha=0.05)

# For HGA vs Control in NSCs
res_HGA_NSCs <- results(dds_combat_seq, contrast = list(c("Condition_HGA_vs_Control","Cell_TypeNSC.ConditionHGA")), alpha=0.05)

# For KS vs Control in Neurons
res_KS_Neurons <- results(dds_combat_seq, contrast = list(c("Condition_KS_vs_Control","Cell_TypeNeurons.ConditionKS")), alpha=0.05)

# For HGA vs Control in Neurons
res_HGA_Neurons <- results(dds_combat_seq, contrast = list(c("Condition_HGA_vs_Control","Cell_TypeNeurons.ConditionHGA")), alpha=0.05)

#________________________SUMMARIZE RESULTS________________________
# Summarize results for each condition and cell type
summary(res_KS_iPSCs)
summary(res_HGA_iPSCs)
summary(res_KS_NSCs)
summary(res_HGA_NSCs)
summary(res_KS_Neurons)
summary(res_HGA_Neurons)

#________________________FILTER RESULTS________________________
# Filter results to keep only non significant transcripts with minimal expression differences (no more than 15%)
res_KS_iPSCs <- res_KS_iPSCs[which(res_KS_iPSCs$padj > 0.05 & abs(res_KS_iPSCs$log2FoldChange) <= 0.05), ]
res_HGA_iPSCs <- res_HGA_iPSCs[which(res_HGA_iPSCs$padj > 0.05 & abs(res_HGA_iPSCs$log2FoldChange) <= 0.05), ]
res_KS_NSCs <- res_KS_NSCs[which(res_KS_NSCs$padj > 0.05 & abs(res_KS_NSCs$log2FoldChange) <= 0.05), ]
res_HGA_NSCs <- res_HGA_NSCs[which(res_HGA_NSCs$padj > 0.05 & abs(res_HGA_NSCs$log2FoldChange) <= 0.05), ]
res_KS_Neurons <- res_KS_Neurons[which(res_KS_Neurons$padj > 0.05 & abs(res_KS_Neurons$log2FoldChange) <= 0.05), ]
res_HGA_Neurons <- res_HGA_Neurons[which(res_HGA_Neurons$padj > 0.05 & abs(res_HGA_Neurons$log2FoldChange) <= 0.05), ]

# Extract the transcript list form the filtered results and merge those from same cell type
transcript_list_iPSCs <- unique(c(rownames(res_KS_iPSCs), rownames(res_HGA_iPSCs)))
length(transcript_list_iPSCs)
transcript_list_NSCs <- unique(c(rownames(res_KS_NSCs), rownames(res_HGA_NSCs)))
length(transcript_list_NSCs)
transcript_list_Neurons <- unique(c(rownames(res_KS_Neurons), rownames(res_HGA_Neurons)))
length(transcript_list_Neurons)

# Map the transcript id to their corresponding gene id with ensembl
# Connect to the Ensembl database
edb <- EnsDb.Hsapiens.v86

# Fetch data using ensembldb
annot_ipscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_list_iPSCs), columns = "gene_id")
# Convert to data frame the annotation
annot_ipscs_ensembldb_df <- as.data.frame(annot_ipscs_ensembldb)
# Extract unique gene ids from the annotation
gene_list_ipscs <- unique(annot_ipscs_ensembldb_df$gene_id)
write.table(gene_list_ipscs, file = "../results_DTE/forTFBS/gene_list_background_ipscs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Repeat with other cell types
annot_nscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_list_NSCs), columns = "gene_id")
annot_nscs_ensembldb_df <- as.data.frame(annot_nscs_ensembldb)
gene_list_nscs <- unique(annot_nscs_ensembldb_df$gene_id)
write.table(gene_list_nscs, file = "../results_DTE/forTFBS/gene_list_background_nscs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

annot_neurons_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_list_Neurons), columns = "gene_id")
annot_neurons_ensembldb_df <- as.data.frame(annot_neurons_ensembldb)
gene_list_neurons <- unique(annot_neurons_ensembldb_df$gene_id)
write.table(gene_list_neurons, file = "../results_DTE/forTFBS/gene_list_background_neurons.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




