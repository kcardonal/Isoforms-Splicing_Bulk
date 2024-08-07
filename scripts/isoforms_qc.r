
#Importing libaries
library(tximport)
library(DESeq2)
library(readxl)
library(ggplot2)
library(patchwork)


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
