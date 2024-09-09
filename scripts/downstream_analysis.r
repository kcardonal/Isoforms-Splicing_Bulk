#________________________DIFFERENTIAL TRANSCRIPT EXPRESSION (DTE) ANALYSIS________________________

#DTE analysis focuses on identifying differentially expressed transcripts between conditions.
# This script performs differential transcript expression (DTE) analysis using the DESeq2 package.
# DTE provides a general overview of transcript-level expression changes across conditions.
# It can help identify key transcripts of interest and provide context for DTU and splicing analyses.

#________________________LOAD REQUIRED LIBRARIES________________________
library(DESeq2)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(biomaRt)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(UpSetR)


#________________________LOAD DATA________________________
# Use combat-seq corrected data for DTE analysis With DESeq2
# This model includes an interaction term, allowing to test how the effect of aneuploidy level across different cell types. 
#This is useful to test the hypothesis that karyotype effects are different in iPSCs versus NSCs or Neurons.

# Load the combat-seq corrected data
load("../Large_Files_No_repo/combat_seq_counts_isoforms.RData")
# Load the sample information
load("../Large_Files_No_repo/isoforms_data.RData")

dds_combat_seq <- DESeqDataSetFromMatrix(countData = iso_corrected,
                                      colData = sample_info,
                                      design = ~ Cell_Type + Condition + Cell_Type:Condition)

#________________________RUN DTE ANALYSIS________________________
# Run DESeq2 analysis
dds_combat_seq <- DESeq(dds_combat_seq)
# Extract results
#list available results names to see all interaction terms
resultsNames(dds_combat_seq)
#[1] "Intercept"                     "Cell_Type_Neurons_vs_iPSCs" - difference between Neurons and iPSCs  
#[3] "Cell_Type_NSC_vs_iPSCs" - difference between NSCs and iPSCs         "Condition_HGA_vs_Control" -  difference between HGA and Control  
#[5] "Condition_KS_vs_Control" -difference between KS And control      "Cell_TypeNeurons.ConditionHGA" - interation effect of HGA in neurons compared to iPSCs
#[7] "Cell_TypeNSC.ConditionHGA" - interaction effect of HGA in NSCs compared to control in iPSCs    "Cell_TypeNeurons.ConditionKS" - interaction effect of KS in neurons compared to iPSCs
#[9] "Cell_TypeNSC.ConditionKS" - interaction effect of KS in NSCs compared to iPSCs control

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
# Filter results to keep only significant transcripts
res_KS_iPSCs <- res_KS_iPSCs[which(res_KS_iPSCs$padj < 0.05), ]
res_HGA_iPSCs <- res_HGA_iPSCs[which(res_HGA_iPSCs$padj < 0.05), ]
res_KS_NSCs <- res_KS_NSCs[which(res_KS_NSCs$padj < 0.05), ]
res_HGA_NSCs <- res_HGA_NSCs[which(res_HGA_NSCs$padj < 0.05), ]
res_KS_Neurons <- res_KS_Neurons[which(res_KS_Neurons$padj < 0.05), ]
res_HGA_Neurons <- res_HGA_Neurons[which(res_HGA_Neurons$padj < 0.05), ]

#________________________ANNOTATE RESULTS USING ensembldb________________________
# Annotate results with gene names, chromosome, transcript biotype, and gene biotype
# Get the transcript IDs from the results
transcript_ids_ks_ipscs <- rownames(res_KS_iPSCs)
# Connect to the Ensembl database
edb <- EnsDb.Hsapiens.v86

# Fetch data using ensembldb
annot_ks_ipscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_ks_ipscs), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
# Convert to data frame the annotation
annot_ks_ipscs_ensembldb_df <- as.data.frame(annot_ks_ipscs_ensembldb)
# Convert to data frame the results
res_KS_iPSCs_df <- as.data.frame(res_KS_iPSCs)
# Merge annotations with results
res_KS_iPSCs_annotated <- merge(res_KS_iPSCs_df, annot_ks_ipscs_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
# Rename the row.names column to transcript_id for clarity
colnames(res_KS_iPSCs_annotated)[1] <- "transcript_id"
# Rename the seqnames column to chromosome for clarity
colnames(res_KS_iPSCs_annotated)[8] <- "chromosome"

# Repeat the same process for the other conditions and cell types
# For HGA in iPSCs
transcript_ids_hga_ipscs <- rownames(res_HGA_iPSCs)
annot_hga_ipscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_hga_ipscs), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
annot_hga_ipscs_ensembldb_df <- as.data.frame(annot_hga_ipscs_ensembldb)
res_HGA_iPSCs_df <- as.data.frame(res_HGA_iPSCs)
res_HGA_iPSCs_annotated <- merge(res_HGA_iPSCs_df, annot_hga_ipscs_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
colnames(res_HGA_iPSCs_annotated)[1] <- "transcript_id"
colnames(res_HGA_iPSCs_annotated)[8] <- "chromosome"

# For KS in NSCs
transcript_ids_ks_nscs <- rownames(res_KS_NSCs)
annot_ks_nscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_ks_nscs), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
annot_ks_nscs_ensembldb_df <- as.data.frame(annot_ks_nscs_ensembldb)
res_KS_NSCs_df <- as.data.frame(res_KS_NSCs)
res_KS_NSCs_annotated <- merge(res_KS_NSCs_df, annot_ks_nscs_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
colnames(res_KS_NSCs_annotated)[1] <- "transcript_id"
colnames(res_KS_NSCs_annotated)[8] <- "chromosome"

# For HGA in NSCs
transcript_ids_hga_nscs <- rownames(res_HGA_NSCs)
annot_hga_nscs_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_hga_nscs), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
annot_hga_nscs_ensembldb_df <- as.data.frame(annot_hga_nscs_ensembldb)
res_HGA_NSCs_df <- as.data.frame(res_HGA_NSCs)
res_HGA_NSCs_annotated <- merge(res_HGA_NSCs_df, annot_hga_nscs_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
colnames(res_HGA_NSCs_annotated)[1] <- "transcript_id"
colnames(res_HGA_NSCs_annotated)[8] <- "chromosome"

# For KS in Neurons
transcript_ids_ks_neurons <- rownames(res_KS_Neurons)
annot_ks_neurons_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_ks_neurons), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
annot_ks_neurons_ensembldb_df <- as.data.frame(annot_ks_neurons_ensembldb)
res_KS_Neurons_df <- as.data.frame(res_KS_Neurons)
res_KS_Neurons_annotated <- merge(res_KS_Neurons_df, annot_ks_neurons_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
colnames(res_KS_Neurons_annotated)[1] <- "transcript_id"
colnames(res_KS_Neurons_annotated)[8] <- "chromosome"

# For HGA in Neurons
transcript_ids_hga_neurons <- rownames(res_HGA_Neurons)
annot_hga_neurons_ensembldb <- transcripts(edb, filter = TxIdFilter(transcript_ids_hga_neurons), columns = c("tx_id", "gene_id", "gene_name", "seq_name", "tx_biotype", "gene_biotype"))
annot_hga_neurons_ensembldb_df <- as.data.frame(annot_hga_neurons_ensembldb)
res_HGA_Neurons_df <- as.data.frame(res_HGA_Neurons)
res_HGA_Neurons_annotated <- merge(res_HGA_Neurons_df, annot_hga_neurons_ensembldb_df, by.x = "row.names", by.y = "tx_id", all.x = TRUE)
colnames(res_HGA_Neurons_annotated)[1] <- "transcript_id"
colnames(res_HGA_Neurons_annotated)[8] <- "chromosome"

#________________________SAVE RESULTS________________________
#SAVE RESULTS
# Save results to the same excel file with multiple sheets
write.xlsx(list("KS_iPSCs" = res_KS_iPSCs_annotated, "HGA_iPSCs" = res_HGA_iPSCs_annotated,
                "KS_NSCs" = res_KS_NSCs_annotated, "HGA_NSCs" = res_HGA_NSCs_annotated,
                "KS_Neurons" = res_KS_Neurons_annotated, "HGA_Neurons" = res_HGA_Neurons_annotated),
           file = "../results/DTE_results.xlsx")

#________________________VISUALIZE RESULTS________________________

#________________________PLOT NUMBER OF SIGNIFICANT TRANSCRIPTS:BARPLOT________________________
# Funtion to count number of significant upregulated and downregulated transcripts
count_sig <- function(res){
  up <- sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)
  down <- sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)
  return(c(Upregulated = up, Downregulated = down))
}
# Count number of significant upregulated and downregulated transcripts for each condition and cell type
counts_ks_ipscs <- count_sig(res_KS_iPSCs)
counts_hga_ipscs <- count_sig(res_HGA_iPSCs)
counts_ks_nscs <- count_sig(res_KS_NSCs)
counts_hga_nscs <- count_sig(res_HGA_NSCs)
counts_ks_neurons <- count_sig(res_KS_Neurons)
counts_hga_neurons <- count_sig(res_HGA_Neurons)

# Create a data frame with the counts for plotting
df_counts <- data.frame(
    Celltype = rep(c("iPSCs", "NSCs", "Neurons"), each = 2),
    condition = rep(c("KS", "HGA"), 3),
    upreg = c(counts_ks_ipscs["Upregulated"], counts_hga_ipscs["Upregulated"],
              counts_ks_nscs["Upregulated"], counts_hga_nscs["Upregulated"],
              counts_ks_neurons["Upregulated"], counts_hga_neurons["Upregulated"]),
    downreg = c(counts_ks_ipscs["Downregulated"], counts_hga_ipscs["Downregulated"],
                counts_ks_nscs["Downregulated"], counts_hga_nscs["Downregulated"],
                counts_ks_neurons["Downregulated"], counts_hga_neurons["Downregulated"])
)
#reshape the dataframe for plotting with ggplot2
df_counts_melt <- df_counts %>% tidyr::pivot_longer(cols = c(upreg, downreg), names_to = "Direction", values_to = "Count")

#Define the order of the cell types for plotting
df_counts_melt$Celltype <- factor(df_counts_melt$Celltype, levels = c("iPSCs", "NSCs", "Neurons"))

#Create barplot
ggplot(df_counts_melt, aes(x = Celltype, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ condition, scales = "free_y") +
    labs(title = "Number of Differentially Expressed Transcripts by condition and cell type",
         x = "Cell Type", y = "Count") +
    theme_minimal() + scale_fill_manual(values = c("upreg" = "gold", "downreg" = "purple"))

#Save the plot
ggsave("./plots/DTE_summary.png", width = 8, height = 6, units = "in")

#________________________PLOT UPSET PLOTS________________________

#Load the DTE results from an excel file
res_KS_iPSCs_annotated <- read.xlsx("../results/DTE_results.xlsx", sheet = "KS_iPSCs", detectDates = FALSE)
res_HGA_iPSCs_annotated <- read.xlsx("../results/DTE_results.xlsx", sheet = "HGA_iPSCs", detectDates = FALSE)
res_KS_NSCs_annotated <- read.xlsx("../results/DTE_results.xlsx", sheet = "KS_NSCs", detectDates = FALSE)
res_HGA_NSCs_annotated <- read.xlsx("../results/DTE_results.xlsx", sheet = "HGA_NSCs", detectDates = FALSE)
res_KS_Neurons_annotated <- read.xlsx("../results/DTE_results.xlsx", sheet = "KS_Neurons", detectDates = FALSE)
res_HGA_Neurons_annotated <- read.xlsx("../results/DTE_results.xlsx", sheet = "HGA_Neurons", detectDates = FALSE)

# Fix rownames
rownames(res_KS_iPSCs_annotated) <- res_KS_iPSCs_annotated$transcript_id
res_KS_iPSCs_annotated <- res_KS_iPSCs_annotated[, -1]
rownames(res_HGA_iPSCs_annotated) <- res_HGA_iPSCs_annotated$transcript_id
res_HGA_iPSCs_annotated <- res_HGA_iPSCs_annotated[, -1]
rownames(res_KS_NSCs_annotated) <- res_KS_NSCs_annotated$transcript_id
res_KS_NSCs_annotated <- res_KS_NSCs_annotated[, -1]
rownames(res_HGA_NSCs_annotated) <- res_HGA_NSCs_annotated$transcript_id
res_HGA_NSCs_annotated <- res_HGA_NSCs_annotated[, -1]
rownames(res_KS_Neurons_annotated) <- res_KS_Neurons_annotated$transcript_id
res_KS_Neurons_annotated <- res_KS_Neurons_annotated[, -1]
rownames(res_HGA_Neurons_annotated) <- res_HGA_Neurons_annotated$transcript_id
res_HGA_Neurons_annotated <- res_HGA_Neurons_annotated[, -1]

# Prepare data for UpSet plots.
# Example function for fetching up/downregulated transcripts
get_up_down <- function(res) {
  upregulated <- rownames(res[which(res$log2FoldChange > 0 & res$padj < 0.05), ])
  downregulated <- rownames(res[which(res$log2FoldChange < 0 & res$padj < 0.05), ])
  return(list(up = upregulated, down = downregulated))
}

# process_common_transcripts(res_KS_iPSCs_up_down, res_HGA_iPSCs_up_down, res_KS_iPSCs_annotated, res_HGA_iPSCs_annotated, "../results/common_upregulated_iPSCs.xlsx")
# iPSCs
res_KS_iPSCs_up_down <- get_up_down(res_KS_iPSCs_annotated)
res_HGA_iPSCs_up_down <- get_up_down(res_HGA_iPSCs_annotated)

# NSCs
res_KS_NSCs_up_down <- get_up_down(res_KS_NSCs_annotated)
res_HGA_NSCs_up_down <- get_up_down(res_HGA_NSCs_annotated)

# Neurons
res_KS_Neurons_up_down <- get_up_down(res_KS_Neurons_annotated)
res_HGA_Neurons_up_down <- get_up_down(res_HGA_Neurons_annotated)

# Combine the up/downregulated lists for each condition and cell type

# iPSCs
upset_ks_ipscs <- list(KS_iPSCs_Up = res_KS_iPSCs_up_down$up,
                       KS_iPSCs_Down = res_KS_iPSCs_up_down$down,
                       HGA_iPSCs_Up = res_HGA_iPSCs_up_down$up,
                       HGA_iPSCs_Down = res_HGA_iPSCs_up_down$down)
# convert the list to a format suitable for UpSetR
upset_data_ipscs <- fromList(upset_ks_ipscs)
#Generate the UpSet plot
pdf(file = "./plots/DTE_upset_ipscs.pdf", width = 10, height = 10)
upset(upset_data_ipscs, sets = names(upset_ks_ipscs),
      keep.order = TRUE,
      order.by = "freq",
      text.scale = 1.5)
dev.off()


# NSCs
upset_ks_nscs <- list(KS_NSCs_Up = res_KS_NSCs_up_down$up,
                      KS_NSCs_Down = res_KS_NSCs_up_down$down,
                      HGA_NSCs_Up = res_HGA_NSCs_up_down$up,
                      HGA_NSCs_Down = res_HGA_NSCs_up_down$down)
# convert the list to a format suitable for UpSetR
upset_data_nscs <- fromList(upset_ks_nscs)
#Generate the UpSet plot
pdf(file = "./plots/DTE_upset_nscs.pdf", width = 10, height = 10)
upset(upset_data_nscs, sets = names(upset_ks_nscs),
      keep.order = TRUE,
      order.by = "freq",
      text.scale = 1.5)
dev.off()

# Neurons
upset_ks_neurons <- list(KS_Neurons_Up = res_KS_Neurons_up_down$up,
                         KS_Neurons_Down = res_KS_Neurons_up_down$down,
                         HGA_Neurons_Up = res_HGA_Neurons_up_down$up,
                         HGA_Neurons_Down = res_HGA_Neurons_up_down$down)
# convert the list to a format suitable for UpSetR
upset_data_neurons <- fromList(upset_ks_neurons)
#Generate the UpSet plot
pdf(file = "./plots/DTE_upset_neurons.pdf", width = 10, height = 10)
upset(upset_data_neurons, sets = names(upset_ks_neurons),
      keep.order = TRUE,
      order.by = "freq",
      text.scale = 1.5)
dev.off()

#Upset plots by karyotype in all cell types
# KS
upset_ks_all <- list(KS_iPSCs_Up = res_KS_iPSCs_up_down$up,
                     KS_iPSCs_Down = res_KS_iPSCs_up_down$down,
                     KS_NSCs_Up = res_KS_NSCs_up_down$up,
                     KS_NSCs_Down = res_KS_NSCs_up_down$down,
                     KS_Neurons_Up = res_KS_Neurons_up_down$up,
                     KS_Neurons_Down = res_KS_Neurons_up_down$down)
# convert the list to a format suitable for UpSetR
upset_data_ks_all <- fromList(upset_ks_all)
#Generate the UpSet plot
pdf(file = "./plots/DTE_upset_ks_all.pdf", width = 10, height = 10)
upset(upset_data_ks_all, sets = names(upset_ks_all),
      keep.order = TRUE,
      order.by = "freq",
      text.scale = 1.5)
dev.off()

# HGA
upset_hga_all <- list(HGA_iPSCs_Up = res_HGA_iPSCs_up_down$up,
                      HGA_iPSCs_Down = res_HGA_iPSCs_up_down$down,
                      HGA_NSCs_Up = res_HGA_NSCs_up_down$up,
                      HGA_NSCs_Down = res_HGA_NSCs_up_down$down,
                      HGA_Neurons_Up = res_HGA_Neurons_up_down$up,
                      HGA_Neurons_Down = res_HGA_Neurons_up_down$down)
# convert the list to a format suitable for UpSetR
upset_data_hga_all <- fromList(upset_hga_all)
#Generate the UpSet plot
pdf(file = "./plots/DTE_upset_hga_all.pdf", width = 10, height = 10)
upset(upset_data_hga_all, sets = names(upset_hga_all),
      keep.order = TRUE,
      order.by = "freq",
      text.scale = 1.5)
dev.off()

#________________________GET THE TRANSCRIPS IDENTITY FROM THE UPSET PLOTS________________________
# Get the transcripts identity from the UpSet plots

process_common_transcripts <- function(res_1, res_2, res_df1, res_df2, dir1, dir2, name1, name2, output_file) {
  # function to process common transcripts between two conditions

  # res_1: list of up/downregulated transcripts for condition 1
  # res_2: list of up/downregulated transcripts for condition 2
  # res_df1: data frame with the results for condition 1
  # res_df2: data frame with the results for condition 2
  # dir1: direction of the transcripts to be compared in res_1 (up or down)
  # dir2: direction of the transcripts to be compared in res_2 (up or down)
  # name1: name to discriminate results that cames from the results in condition 1 (e.g., "_KS")
  # name2: name to discriminate results that cames from the results in condition 2 (e.g., "_HGA")
  # output_file: path to save the output file as excel


  # Calculate the intersection of upregulated transcripts
  common <- intersect(res_1[[dir1]], res_2[[dir2]])
  
  # Create a function to map the common transcripts on the original results, and extract in a new data frame
  map_common_transcripts <- function(common, res) {
    # Check if res is a data frame
    if (!is.data.frame(res)) {
      stop("The input 'res' is not a data frame.")
    }
    
    # Check if common contains valid row names
    valid_rows <- rownames(res) %in% common
    if (sum(valid_rows) == 0) {
      stop("No common transcripts found in the row names of 'res'.")
    }
    
    common_res <- res[valid_rows, ]
    return(common_res)
  }
  
  # Map the common upregulated transcripts to the original results
  common_1 <- map_common_transcripts(common, res_df1)
  common_2 <- map_common_transcripts(common, res_df2)
  
  # Merge by rowname the results in a single data frame
  common_merged <- merge(common_1, common_2, by = "row.names", all = TRUE)
  
  # Rename the row.names column to transcript_id
  colnames(common_merged)[1] <- "transcript_id"
  # Rename .x and .y columns to indicate the source
  colnames(common_merged) <- gsub("\\.x", name1, colnames(common_merged))
  colnames(common_merged) <- gsub("\\.y", name2, colnames(common_merged))
  # Remove the redundant columns
  common_merged <- common_merged[, -c(23:31)]
  # Save the results
  write.xlsx(common_merged, file = output_file)
}

# Process common transcripts in iPSCs:
# KS vs HGA in Upregulated transcripts
process_common_transcripts(res_KS_iPSCs_up_down, res_HGA_iPSCs_up_down, res_KS_iPSCs_annotated, res_HGA_iPSCs_annotated, "up", "up", "_KS", "_HGA", "../results/upset_upregulated_iPSCs.xlsx")
# KS vs HGA in Downregulated transcripts
process_common_transcripts(res_KS_iPSCs_up_down, res_HGA_iPSCs_up_down, res_KS_iPSCs_annotated, res_HGA_iPSCs_annotated, "down", "down", "_KS", "_HGA", "../results/upset_downregulated_iPSCs.xlsx")
# KS vs HGA upregulated to downregulated transcripts
process_common_transcripts(res_KS_iPSCs_up_down, res_HGA_iPSCs_up_down, res_KS_iPSCs_annotated, res_HGA_iPSCs_annotated, "up", "down", "_KS", "_HGA", "../results/upset_up_to_down_iPSCs.xlsx")
# KS vs HGA downregulated to upregulated transcripts
process_common_transcripts(res_KS_iPSCs_up_down, res_HGA_iPSCs_up_down, res_KS_iPSCs_annotated, res_HGA_iPSCs_annotated, "down", "up", "_KS", "_HGA", "../results/upset_down_to_up_iPSCs.xlsx")

# Process common transcripts in NSCs:
# KS vs HGA in Upregulated transcripts
process_common_transcripts(res_KS_NSCs_up_down, res_HGA_NSCs_up_down, res_KS_NSCs_annotated, res_HGA_NSCs_annotated, "up", "up", "_KS", "_HGA", "../results/upset_upregulated_NSCs.xlsx")
# KS vs HGA in Downregulated transcripts
process_common_transcripts(res_KS_NSCs_up_down, res_HGA_NSCs_up_down, res_KS_NSCs_annotated, res_HGA_NSCs_annotated, "down", "down", "_KS", "_HGA", "../results/upset_downregulated_NSCs.xlsx")
# KS vs HGA upregulated to downregulated transcripts
process_common_transcripts(res_KS_NSCs_up_down, res_HGA_NSCs_up_down, res_KS_NSCs_annotated, res_HGA_NSCs_annotated, "up", "down", "_KS", "_HGA", "../results/upset_up_to_down_NSCs.xlsx")
# KS vs HGA downregulated to upregulated transcripts
process_common_transcripts(res_KS_NSCs_up_down, res_HGA_NSCs_up_down, res_KS_NSCs_annotated, res_HGA_NSCs_annotated, "down", "up", "_KS", "_HGA", "../results/upset_down_to_up_NSCs.xlsx")

# Process common transcripts in Neurons:
# KS vs HGA in Upregulated transcripts
process_common_transcripts(res_KS_Neurons_up_down, res_HGA_Neurons_up_down, res_KS_Neurons_annotated, res_HGA_Neurons_annotated, "up", "up", "_KS", "_HGA", "../results/upset_upregulated_Neurons.xlsx")
# KS vs HGA in Downregulated transcripts
process_common_transcripts(res_KS_Neurons_up_down, res_HGA_Neurons_up_down, res_KS_Neurons_annotated, res_HGA_Neurons_annotated, "down", "down", "_KS", "_HGA", "../results/upset_downregulated_Neurons.xlsx")
# KS vs HGA upregulated to downregulated transcripts
process_common_transcripts(res_KS_Neurons_up_down, res_HGA_Neurons_up_down, res_KS_Neurons_annotated, res_HGA_Neurons_annotated, "up", "down", "_KS", "_HGA", "../results/upset_up_to_down_Neurons.xlsx")
# KS vs HGA downregulated to upregulated transcripts
process_common_transcripts(res_KS_Neurons_up_down, res_HGA_Neurons_up_down, res_KS_Neurons_annotated, res_HGA_Neurons_annotated, "down", "up", "_KS", "_HGA", "../results/upset_down_to_up_Neurons.xlsx")

# Independent processing to get the common transcripts between all cell types by karyotype
# KS upregulated transcripts
common_ks_up <- Reduce(intersect, list(res_KS_iPSCs_up_down$up, res_KS_NSCs_up_down$up, res_KS_Neurons_up_down$up))
# KS downregulated transcripts
common_ks_down <- Reduce(intersect, list(res_KS_iPSCs_up_down$down, res_KS_NSCs_up_down$down, res_KS_Neurons_up_down$down))
# HGA upregulated transcripts
common_hga_up <- Reduce(intersect, list(res_HGA_iPSCs_up_down$up, res_HGA_NSCs_up_down$up, res_HGA_Neurons_up_down$up))
# HGA downregulated transcripts
common_hga_down <- Reduce(intersect, list(res_HGA_iPSCs_up_down$down, res_HGA_NSCs_up_down$down, res_HGA_Neurons_up_down$down))

# Map the common transcripts to the original results

#upregulated KS
common_ks_up_i  <- map_common_transcripts(common_ks_up, res_KS_iPSCs_annotated)
common_ks_up_ns  <- map_common_transcripts(common_ks_up, res_KS_NSCs_annotated)
common_ks_up_n <- map_common_transcripts(common_ks_up, res_KS_Neurons_annotated)
# Merge the three sets in a single data frame
common_ks_up_merged <- merge(common_ks_up_i, common_ks_up_ns, by = "row.names", all = TRUE)
# Rename the row.names column to transcript_id
colnames(common_ks_up_merged)[1] <- "transcript_id"
# Make transcript_id the rownames
rownames(common_ks_up_merged) <- common_ks_up_merged$transcript_id
# do the last merge
common_ks_up_merged <- merge(common_ks_up_merged, common_ks_up_n, by = "row.names", all = TRUE)
# Remove the Row.names column
common_ks_up_merged <- common_ks_up_merged[, -1]
# make transcript_id the rownames
rownames(common_ks_up_merged) <- common_ks_up_merged$transcript_id
# Remove the transcript_id column
common_ks_up_merged <- common_ks_up_merged[, -1]
# Remove redundant columns
common_ks_up_merged <- common_ks_up_merged[, -c(7:15, 22:30)]
# Rename the columns to indicate the source
colnames(common_ks_up_merged) <- gsub("\\.x", "_iPSCs", colnames(common_ks_up_merged))
colnames(common_ks_up_merged) <- gsub("\\.y", "_NSCs", colnames(common_ks_up_merged))
# Save the results
write.xlsx(common_ks_up_merged, file = "../results/common_upregulated_KS_all.xlsx")

#downregulated KS
common_ks_down_i <- map_common_transcripts(common_ks_down, res_KS_iPSCs_annotated)
common_ks_down_ns <- map_common_transcripts(common_ks_down, res_KS_NSCs_annotated)
common_ks_down_n <- map_common_transcripts(common_ks_down, res_KS_Neurons_annotated)
# Merge the three sets in a single data frame
common_ks_down_merged <- merge(common_ks_down_i, common_ks_down_ns, by = "row.names", all = TRUE)
# Rename the row.names column to transcript_id
colnames(common_ks_down_merged)[1] <- "transcript_id"
# Make transcript_id the rownames
rownames(common_ks_down_merged) <- common_ks_down_merged$transcript_id
# do the last merge
common_ks_down_merged <- merge(common_ks_down_merged, common_ks_down_n, by = "row.names", all = TRUE)
# Remove the Row.names column
common_ks_down_merged <- common_ks_down_merged[, -1]
# make transcript_id the rownames
rownames(common_ks_down_merged) <- common_ks_down_merged$transcript_id
# Remove the transcript_id column
common_ks_down_merged <- common_ks_down_merged[, -1]
# Remove redundant columns
common_ks_down_merged <- common_ks_down_merged[, -c(7:15, 22:30)]
# Rename the columns to indicate the source
colnames(common_ks_down_merged) <- gsub("\\.x", "_iPSCs", colnames(common_ks_down_merged))
colnames(common_ks_down_merged) <- gsub("\\.y", "_NSCs", colnames(common_ks_down_merged))
# Save the results
write.xlsx(common_ks_down_merged, file = "../results/common_downregulated_KS_all.xlsx")

#upregulated HGA
common_hga_up_i <- map_common_transcripts(common_hga_up, res_HGA_iPSCs_annotated)
common_hga_up_ns <- map_common_transcripts(common_hga_up, res_HGA_NSCs_annotated)
common_hga_up_n <- map_common_transcripts(common_hga_up, res_HGA_Neurons_annotated)
# Merge the three sets in a single data frame
common_hga_up_merged <- merge(common_hga_up_i, common_hga_up_ns, by = "row.names", all = TRUE)
# Rename the row.names column to transcript_id
colnames(common_hga_up_merged)[1] <- "transcript_id"
# Make transcript_id the rownames
rownames(common_hga_up_merged) <- common_hga_up_merged$transcript_id
# do the last merge
common_hga_up_merged <- merge(common_hga_up_merged, common_hga_up_n, by = "row.names", all = TRUE)
# Remove the Row.names column
common_hga_up_merged <- common_hga_up_merged[, -1]
# make transcript_id the rownames
rownames(common_hga_up_merged) <- common_hga_up_merged$transcript_id
# Remove the transcript_id column
common_hga_up_merged <- common_hga_up_merged[, -1]
# Remove redundant columns
common_hga_up_merged <- common_hga_up_merged[, -c(7:15, 22:30)]
# Rename the columns to indicate the source
colnames(common_hga_up_merged) <- gsub("\\.x", "_iPSCs", colnames(common_hga_up_merged))
colnames(common_hga_up_merged) <- gsub("\\.y", "_NSCs", colnames(common_hga_up_merged))
# Save the results
write.xlsx(common_hga_up_merged, file = "../results/common_upregulated_HGA_all.xlsx")

#downregulated HGA
common_hga_down_i <- map_common_transcripts(common_hga_down, res_HGA_iPSCs_annotated)
common_hga_down_ns <- map_common_transcripts(common_hga_down, res_HGA_NSCs_annotated)
common_hga_down_n <- map_common_transcripts(common_hga_down, res_HGA_Neurons_annotated)
# Merge the three sets in a single data frame
common_hga_down_merged <- merge(common_hga_down_i, common_hga_down_ns, by = "row.names", all = TRUE)
# Rename the row.names column to transcript_id
colnames(common_hga_down_merged)[1] <- "transcript_id"
# Make transcript_id the rownames
rownames(common_hga_down_merged) <- common_hga_down_merged$transcript_id
# do the last merge
common_hga_down_merged <- merge(common_hga_down_merged, common_hga_down_n, by = "row.names", all = TRUE)
# Remove the Row.names column
common_hga_down_merged <- common_hga_down_merged[, -1]
# make transcript_id the rownames
rownames(common_hga_down_merged) <- common_hga_down_merged$transcript_id
# Remove the transcript_id column
common_hga_down_merged <- common_hga_down_merged[, -1]
# Remove redundant columns
common_hga_down_merged <- common_hga_down_merged[, -c(7:15, 22:30)]
# Rename the columns to indicate the source
colnames(common_hga_down_merged) <- gsub("\\.x", "_iPSCs", colnames(common_hga_down_merged))
colnames(common_hga_down_merged) <- gsub("\\.y", "_NSCs", colnames(common_hga_down_merged))
# Save the results
write.xlsx(common_hga_down_merged, file = "../results/common_downregulated_HGA_all.xlsx")
