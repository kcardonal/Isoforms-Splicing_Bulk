#August 07, 2024
#Kelly J. Cardona L.

#PCA analysis of isoform level data to check for batch effects and other sources of variation

#Importing libaries
library(tximport)
library(DESeq2)
library(readxl)
library(ggplot2)
library(patchwork)
library(sva)
library(pheatmap)
library(cluster)

#_______________LOADING DATA_____________________

#load Isoform level data. Since i have already run the RSEM pipeline, i will use the isoforms.results files 
#to load the data. The tximport function will load the data and create a SummarizedExperiment object.
#The tximport function will also calculate the counts from the abundance data.

files_iso <- list.files(path = "../counts", pattern = ".isoforms.results", full.names = TRUE)
sample_names_iso <- gsub(".isoforms.results", "", basename(files_iso))
names(files_iso) <- sample_names_iso #assign the sample names to the files

txi_iso <- tximport(files_iso, type = "rsem", txIn = TRUE, txOut = TRUE, countsFromAbundance = "lengthScaledTPM")
#get the samples names in the order they appear on the txt_iso object
samples_in_txi_iso <- colnames(txi_iso$counts)

#Load metadata and match it with the samples in the txi_iso object
sample_info <- read_excel("../metadata_Isoforms_01_07_24.xlsx", sheet = "Metadata")
#convert the sample_info object to a DataFrame object
sample_info <- as.data.frame(sample_info)
rownames(sample_info) <- sample_info$Sample_ID
#remove the Sample_ID column
sample_info <- sample_info[, -1]

#Reorder the sample_info object to match the order of the samples in the txi_iso object
sample_info <- sample_info[samples_in_txi_iso, ]
#check if the sample names are the same in the sample_info and txi_iso objects
all(colnames(txi_iso$counts) == rownames(sample_info)) 
#TRUE
#Save the object with the counts and the sample information
save(txi_iso, sample_info, file = "./data/isoforms_data.RData")

#_______________PCA ANALYSIS BEFORE BATCH CORRECTION_____________________

#Create a DESeqDataSet object
#The design formula will include:
#Batch: controls for variability due to different sequencing runs
#karyotype: tests for effect of different karotypes in isoform expression
#cell type: tests for effect of different cell types in isoform expression
#karyotype:cell type: test wehather the effect of karyotype on isoform expression is different in different cell types

design_formula <- ~ Sequencing_Batch + Cell_Type + Karyotype + Karyotype:Cell_Type
dds_iso <- DESeqDataSetFromTximport(txi_iso, colData = sample_info, design = design_formula)

#Apply variance stabilizing transformation to stabilize the variance across the mean for PCA and other plots
vsd_iso <- vst(dds_iso, blind = FALSE)
#Calculate PCA
pcaData_iso <- plotPCA(vsd_iso, intgroup=c("Sequencing_Batch", "Karyotype", "Cell_Type"), returnData=TRUE)
percentVar_iso <- round(100 * attr(pcaData_iso, "percentVar"))
#Save PCA data and variance explained
save(pcaData_iso, percentVar_iso, file = "./data/pca_isoforms.RData")

#Plot PCA
#PCA colored by sequencing batch
pca_batch <- ggplot(pcaData_iso, aes(x = PC1, y = PC2, color = Sequencing_Batch)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_iso[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_iso[2], "% variance")) +
  ggtitle("PCA by Batch") +
  theme_minimal()

#PCA colored by Karyotype
pca_karyotype <- ggplot(pcaData_iso, aes(x = PC1, y = PC2, color = Karyotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_iso[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_iso[2], "% variance")) +
  ggtitle("PCA by Karyotype") +
  theme_minimal()

#PCA colored by Cell Type
pca_cell_type <- ggplot(pcaData_iso, aes(x = PC1, y = PC2, color = Cell_Type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_iso[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_iso[2], "% variance")) +
  ggtitle("PCA by Cell Type") +
  theme_minimal()

#Combine PCA plots in a single plot 
pca_iso <- pca_batch + pca_karyotype + pca_cell_type + plot_layout(ncol = 3)
#Display the PCA plot
print(pca_iso)
#Save the PCA plot
ggsave("./plots/PCA_isoforms.png", pca_iso, width = 12, height = 6)

#_______________BATCH CORRECTION WITH COMBAT-SEQ_____________________

#applying ComBat-Seq to correct for batch effects in the isoform level data.
#Combat-Seq works with raw RNA-seq counts and correct for batch effects without transforming the data. 
#This method is particularly useful to maintain the original count distribution.

iso_corrected <- ComBat_seq(counts = counts(dds_iso), batch = sample_info$Sequencing_Batch)
#Create a DESeqDataSet object with the corrected data
dds_iso_corrected <- DESeqDataSetFromMatrix(countData = iso_corrected, colData = sample_info, design = design_formula)

#_______________PCA ANALYSIS AFTER BATCH CORRECTION_____________________
#vst transformation for PCA plotting 
vsd_iso_corrected <- vst(dds_iso_corrected, blind = FALSE)
#Calculate PCA
pcaData_iso_corrected <- plotPCA(vsd_iso_corrected, intgroup=c("Sequencing_Batch", "Karyotype", "Cell_Type"), returnData=TRUE)
percentVar_iso_corrected <- round(100 * attr(pcaData_iso_corrected, "percentVar"))
#Save PCA data and variance explained
save(pcaData_iso_corrected, percentVar_iso_corrected, file = "./data/pca_iso_combat_seq.RData")

#Plot PCA by batch
pca_batch_corrected <- ggplot(pcaData_iso_corrected, aes(x = PC1, y = PC2, color = Sequencing_Batch)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_iso_corrected[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_iso_corrected[2], "% variance")) +
  ggtitle("PCA by Batch (Corrected)") +
  theme_minimal()
#Plot PCA by Karyotype
pca_karyotype_corrected <- ggplot(pcaData_iso_corrected, aes(x = PC1, y = PC2, color = Karyotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_iso_corrected[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_iso_corrected[2], "% variance")) +
  ggtitle("PCA by Karyotype (Corrected)") +
  theme_minimal()
 #Plot PCA by Cell Type
pca_cell_type_corrected <- ggplot(pcaData_iso_corrected, aes(x = PC1, y = PC2, color = Cell_Type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_iso_corrected[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_iso_corrected[2], "% variance")) +
  ggtitle("PCA by Cell Type (Corrected)") +
  theme_minimal()
 #Combine PCA plots in a single plot
pca_iso_corrected <- pca_batch_corrected + pca_karyotype_corrected + pca_cell_type_corrected + plot_layout(ncol = 3)
#Display the PCA plot
print(pca_iso_corrected)
#Save the PCA plot
ggsave("./plots/PCA_isoforms_corrected.png", pca_iso_corrected, width = 12, height = 6)  

#_______________BATCH CORRECTION COMBAT_____________________

#Combat is used in already transformed RNA-seq data (in my case VST/rlog from DESeq2) 
#and need to remove batch effects on this transformed data.

#Extract vst counts to use in ComBat
vst_counts <- assay(vsd_iso)
#Create a model matrix to include additional covariates in the batch correction process.
#In this case, we will include the karyotype and cell type as covariates.
model_matrix <- model.matrix(~ Karyotype + Cell_Type + Karyotype:Cell_Type, data = sample_info)
#Apply ComBat to the vst counts
vst_corrected <- ComBat(dat = vst_counts, batch = sample_info$Sequencing_Batch, mod = model_matrix)
#Calculate PCA
pcaData_combat <- prcomp(t(vst_corrected), center=TRUE, scale.=TRUE)
#Extract PCA results for plotting
pcaData_combat <- as.data.frame(pcaData_combat$x)
#Add sample information to the PCA results
pcaData_combat$Sequencing_Batch <- sample_info$Sequencing_Batch
pcaData_combat$Karyotype <- sample_info$Karyotype
pcaData_combat$Cell_Type <- sample_info$Cell_Type
#Calculate the percentage of variance explained by each PC
percentVar_combat <- (pcaData_combat$sdev^2) / sum(pcaData_combat$sdev^2) * 100
#Save PCA data and variance explained
save(pcaData_combat, percentVar_combat, file = "./data/pca_iso_combat.RData")
#Plot PCA
#PCA colored by sequencing batch
pca_batch_combat <- ggplot(pcaData_combat, aes(x = PC1, y = PC2, color = Sequencing_Batch)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar_combat[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_combat[2],2), "% variance")) +
  ggtitle("PCA by Batch (Combat)") +
  theme_minimal()
#PCA colored by Karyotype
pca_karyotype_combat <- ggplot(pcaData_combat, aes(x = PC1, y = PC2, color = Karyotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar_combat[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_combat[2], 2), "% variance")) +
  ggtitle("PCA by Karyotype (Combat)") +
  theme_minimal()
#PCA colored by Cell Type
pca_cell_type_combat <- ggplot(pcaData_combat, aes(x = PC1, y = PC2, color = Cell_Type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar_combat[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_combat[2], 2), "% variance")) +
  ggtitle("PCA by Cell Type (Combat)") +
  theme_minimal()
#Combine PCA plots in a single plot
pca_combat <- pca_batch_combat + pca_karyotype_combat + pca_cell_type_combat + plot_layout(ncol = 3)
#Display the PCA plot
print(pca_combat)
#Save the PCA plot
ggsave("./plots/PCA_corrected_combat.png", pca_combat, width = 12, height = 6)


#_______________ASSESS BATCH EFFECT REMOVAL_____________________

#Create heatmaps showing the distance or correlation between samples before and after batch correction.
#To visualize the similarity between replicates and see if batch correction has effectively reduced 
#unwanted batch variation.

#Define colors for annotations
ann_colors <- list(
  Sequencing_Batch = c("April_2019" = "red", "April_2023" = "brown", "August_2023" = "yellow", "December_2023" = "darkgreen", "January_2020" = "salmon", "July_2018" = "orange", "July_2023" = "black", "March_2021" = "violet", "March_2022" = "cyan", "November_2018" ="magenta"),
  Karyotype = c("46,XY" = "green", "47,XXY" = "lightblue", "48,XXXY" = "blue", "49,XXXXY" = "purple"),
  Cell_Type = c("iPSCs" = "orange", "NSC" = "pink", "Neurons" = "gray")
)

#Extracting relevant information from the sample_info object
annotation_df <- sample_info[, c("Sequencing_Batch", "Karyotype", "Cell_Type")]

#Sample distance heatmap before batch correction
sample_dist_before <- dist(t(assay(vsd_iso)))
sample_dist_before <- as.matrix(sample_dist_before)

#Ensure the annotation dataframe row names match the sample names
rownames(annotation_df) <- colnames(sample_dist_before)
#Create the heatmap
hm_raw <- pheatmap(sample_dist_before,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = "Sample Distance Before Batch Correction",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize = 10,
  #fontsize_row = 8,
  #fontsize_col = 8,
  cellwidth = 10,
  cellheight = 10,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA
  #filename = "./plots/sample_dist_raw_no_labels.png"
)

#Sample distance heatmap after batch correction with ComBat-Seq
sample_dist_combatseq <- dist(t(iso_corrected))
sample_dist_combatseq <- as.matrix(sample_dist_combatseq)
#Create the heatmap
hm_combatseq <- pheatmap(sample_dist_combatseq,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = "Sample Distance After Batch Correction (ComBat-Seq)",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  #fontsize_row  = 8,
  #fontsize_col = 8,
  fontsize = 10,
  cellwidth = 10,
  cellheight = 10,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA
  #filename = "./plots/sample_dist_combatSeq_no_labels.png"
)

#Sample distance heatmap after batch correction with ComBat
sample_dist_combat <- dist(t(vst_corrected))
sample_dist_combat <- as.matrix(sample_dist_combat)
#Create the heatmap
hm_combat <- pheatmap(sample_dist_combat,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = "Sample Distance After Batch Correction (ComBat)",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row  = 8,
  fontsize_col = 8,
  #fontsize = 10,
  cellwidth = 10,
  cellheight = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA
  #filename = "./plots/sample_dist_combat.png"
)

#_________________BATCHEFFECT ASSESMENT METRICS_____________________

#Silhouette width analysis: a measure of how similar an object is to its own cluster compared to other clusters.
#This metric ranges from -1 to 1, where a value close to 1 indicates that the object is well matched to its own cluster.
#A value close to -1 indicates that the object is poorly matched to its own cluster.

# Assuming 'dist_raw' is your pheatmap object
# Extract the order of samples
sample_order <- hm_raw$tree_row$order
# Extract cluster assignments; this uses cutree to determine cluster groups
# Adjust 'k' based on expected number of biological clusters or what you used
clusters <- cutree(dist_raw$tree_row, k = 3)  # For example, 3 clusters
# Map the cluster assignments back to the sample names
cluster_assignments <- data.frame(
  Sample = names(clusters),
  Cluster = as.factor(clusters)
)
# Check the head of cluster_assignments to confirm
head(cluster_assignments)
ordered_clusters <- cluster_assignments[order(cluster_assignments$Sample), "Cluster"]
# Calculate silhouette widths
sil_widths <- silhouette(ordered_clusters, sample_dist_before)


















#Silhouette width before batch correction
sil_width_before <- silhouette(dist(sample_dist_before), as.factor(annotation_df$Sequencing_Batch))
mean_sil_width_before <- mean(sil_width_before[, "sil_width"])
cat("Mean Silhouette Width Before Batch Correction: ", mean_sil_width_before, "\n")

#Silhouette width after batch correction with ComBat-Seq
sil_width_combatseq <- silhouette(as.dist(sample_dist_combatseq), as.factor(annotation_df$Sequencing_Batch))
mean_sil_width_combatseq <- mean(sil_width_combatseq[, "sil_width"])
cat("Mean Silhouette Width After Batch Correction (ComBat-Seq): ", mean_sil_width_combatseq, "\n")

#Silhouette width after batch correction with ComBat
sil_width_combat <- silhouette(as.dist(sample_dist_combat), as.factor(annotation_df$Sequencing_Batch))
mean_sil_width_combat <- mean(sil_width_combat[, "sil_width"])
cat("Mean Silhouette Width After Batch Correction (ComBat): ", mean_sil_width_combat, "\n")