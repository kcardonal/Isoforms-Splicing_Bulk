#Script to draw volcano plots for differentially expressed transcripts in common between the different contrasts

#load packages needed
library(ggplot2)
library(reshape2)

# Load the results of the common differentially expressed transcripts
#ipscs up
print(det_ipscs_up)

melted_data <- melt(det_ipscs_up,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes

# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated DETs Across iPSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "blue", "log2FoldChange_HGA" = "darkblue"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_ipscs_up.png", width = 8, height = 6, dpi = 300)

#ipscs down
melted_data <- melt(det_ipscs_down,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes

# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated DETs Across iPSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "blue", "log2FoldChange_HGA" = "darkblue"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_ipscs_down.png", width = 8, height = 6, dpi = 300)

#ipscs up to down
melted_data <- melt(det_ipscs_up_to_down,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated to Downregulated DETs Across iPSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "blue", "log2FoldChange_HGA" = "darkblue"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_ipscs_up_to_down.png", width = 8, height = 6, dpi = 300)

#ipscs down to up
melted_data <- melt(det_ipscs_down_to_up,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated to Upregulated DETs Across iPSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "blue", "log2FoldChange_HGA" = "darkblue"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_ipscs_up_to_down.png", width = 8, height = 6, dpi = 300)

#nscs up
melted_data <- melt(det_nscs_up,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated DETs Across NSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "green", "log2FoldChange_HGA" = "darkgreen"))  # Customize colors
 # Save the plot
ggsave("./plots/boxplot_det_common_nscs_up.png", width = 8, height = 6, dpi = 300) 

#nscs down
melted_data <- melt(det_nscs_down,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated DETs Across NSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "green", "log2FoldChange_HGA" = "darkgreen"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_nscs_down.png", width = 8, height = 6, dpi = 300)

#nscs up to down
melted_data <- melt(det_nscs_up_to_down,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated to Downregulated DETs Across NSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "green", "log2FoldChange_HGA" = "darkgreen"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_nscs_up_to_down.png", width = 8, height = 6, dpi = 300)

#nscs down to up
melted_data <- melt(det_nscs_down_to_up,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated to Upregulated DETs Across NSCs",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "green", "log2FoldChange_HGA" = "darkgreen"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_nscs_down_to_up.png", width = 8, height = 6, dpi = 300)

#neurons up
melted_data <- melt(det_neurons_up,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated DETs Across Neurons",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "violet", "log2FoldChange_HGA" = "purple"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_neurons_up.png", width = 8, height = 6, dpi = 300)

#neurons down
melted_data <- melt(det_neurons_down,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated DETs Across Neurons",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "violet", "log2FoldChange_HGA" = "purple"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_neurons_down.png", width = 8, height = 6, dpi = 300)

#neurons up to down
melted_data <- melt(det_neurons_up_to_down,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated DETs Across Neurons",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "violet", "log2FoldChange_HGA" = "purple"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_neurons_down.png", width = 8, height = 6, dpi = 300)

#neurons down to up
melted_data <- melt(det_neurons_down_to_up,
                    id.vars = "transcript_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_KS", "log2FoldChange_HGA"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated to Upregulated DETs Across Neurons",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_KS" = "violet", "log2FoldChange_HGA" = "purple"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_neurons_down_to_up.png", width = 8, height = 6, dpi = 300)

#ks up
melted_data <- melt(det_ks_up,
                    id.vars = "gene_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_iPSCs", "log2FoldChange_NSCs", "log2FoldChange_Neurons"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated DETs Across KS",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_iPSCs" = "red", "log2FoldChange_NSCs" = "darkred", "log2FoldChange_Neurons" = "brown"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_ks_up.png", width = 8, height = 6, dpi = 300)

#ks down
melted_data <- melt(det_ks_down,
                    id.vars = "gene_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_iPSCs", "log2FoldChange_NSCs", "log2FoldChange_Neurons"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated DETs Across KS",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_iPSCs" = "red", "log2FoldChange_NSCs" = "darkred", "log2FoldChange_Neurons" = "brown"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_ks_down.png", width = 8, height = 6, dpi = 300)

#hga up
melted_data <- melt(det_hga_up,
                    id.vars = "gene_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_iPSCs", "log2FoldChange_NSCs", "log2FoldChange_Neurons"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Upregulated DETs Across HGA",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_iPSCs" = "yellow", "log2FoldChange_NSCs" = "orange", "log2FoldChange_Neurons" = "darkorange"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_hga_up.png", width = 8, height = 6, dpi = 300)

#hga down
melted_data <- melt(det_hga_down,
                    id.vars = "gene_id",  # The unique transcript identifier
                    measure.vars = c("log2FoldChange_iPSCs", "log2FoldChange_NSCs", "log2FoldChange_Neurons"),  # logFC columns
                    variable.name = "Contrast",  # New column to distinguish contrasts
                    value.name = "Log2FoldChange")  # Column to hold log2 fold changes
# Check the reshaped data
head(melted_data)
# Create the boxplot to compare log2FC for common transcripts in both contrasts
ggplot(melted_data, aes(x = Contrast, y = Log2FoldChange, fill = Contrast)) +
  geom_boxplot(alpha = 0.7) +  # Create boxplot
  labs(title = "Common Downregulated DETs Across HGA",
       x = "Contrast",
       y = "Log2 Fold Change") +
  theme_minimal() +
  scale_fill_manual(values = c("log2FoldChange_iPSCs" = "yellow", "log2FoldChange_NSCs" = "orange", "log2FoldChange_Neurons" = "darkorange"))  # Customize colors
# Save the plot
ggsave("./plots/boxplot_det_common_hga_down.png", width = 8, height = 6, dpi = 300)