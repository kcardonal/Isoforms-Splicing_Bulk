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

#_______________PCA ANALYSIS BEFORE BATCH CORRECTION_____________________

#Create a DESeqDataSet object
#The design formula will include:
#Batch: controls for variability due to different sequencing runs
#karyotype: tests for effect of different karotypes in isoform expression
#cell type: tests for effect of different cell types in isoform expression
#karyotype:cell type: test wehather the effect of karyotype on isoform expression is different in different cell types

design_formula <- ~ Sequencing_Batch + Cell_Type + Karyotype + Karyotype:Cell_Type
dds_iso <- DESeqDataSetFromTximport(txi_iso, colData = sample_info, design = design_formula)

#Apply varaince stabilizing transformation to stabilize the variance across the mean for PCA and other plots
vsd_iso <- vst(dds_iso, blind = FALSE)
#Calculate PCA
pcaData_iso <- plotPCA(vsd_iso, intgroup=c("Sequencing_Batch", "Karyotype", "Cell_Type"), returnData=TRUE)
percentVar_iso <- round(100 * attr(pcaData_iso, "percentVar"))
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

#_______________BATCH CORRECTION_____________________
#applying ComBat-Seq using the mod argument is more advanced and can be effective for my complex experimental design.
#The mod matrix allows me to include other variables (like karyotype and cell type) in the batch correction process, 
#helping to adjust for batch effects while preserving the biological variability associated with these factors.
# this is to ensure that batch correction is done without diminishing the biological signal associated with my primary variables.
iso_corrected <- ComBat_seq(counts = counts(dds_iso), batch = sample_info$Sequencing_Batch)



#Create a DESeqDataSet object with the corrected data
dds_iso_corrected <- DESeqDataSetFromMatrix(countData = iso_corrected, colData = sample_info, design = design_formula)

#_______________PCA ANALYSIS AFTER BATCH CORRECTION_____________________
#vst transformation for PCA plotting 
vsd_iso_corrected <- vst(dds_iso_corrected, blind = FALSE)
#Calculate PCA
pcaData_iso_corrected <- plotPCA(vsd_iso_corrected, intgroup=c("Sequencing_Batch", "Karyotype", "Cell_Type"), returnData=TRUE)
percentVar_iso_corrected <- round(100 * attr(pcaData_iso_corrected, "percentVar"))
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