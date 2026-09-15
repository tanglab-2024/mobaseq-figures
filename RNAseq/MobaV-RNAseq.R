###TCR##RP48 RNAseq code

###TCR### load packages
##or if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
  #BiocManager::install("rWikiPathways")

#BiocManager::install("DESeq2")

library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(tidyr)
library(ggrepel)
library(ggrastr)
library(clusterProfiler)
library(stringr)
library(corrplot)
library(heatmaply)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(apeglm)
library(grid)
library(lazyeval)
#library(RVA)

###TCR### setwd and load featurecount data
setwd("/Users/tanglab/Library/CloudStorage/Box-Box/AX_RNASeq_Data/MobaV")
setwd("/Users/tanglab/Library/CloudStorage/Box-Box/AX_data_analysis/AX_RNASeq_Data/MobaV")
setwd("C:/Users/Andy/Box/AX_data_analysis/AX_RNASeq_Data/MobaV")
#setwd("C:/Users/Andy/Box/AX_RNASeq_Data/CXCL_mCherry")
file_path <- "5_Counts/2025_MobaV_RNASeq_counts.txt"
countdata <- read.table(file_path, sep="\t", header = T, row.names = 1)
summary_path <- "5_Counts/2025_MobaV_RNASeq_counts.txt.summary"
summary_data <- read.table(summary_path, sep = "\t", header = T, row.names = 1)
project_code <- sub("_counts.*$", "", basename(file_path))

###AX### generate QC plot from summary data
colnames(summary_data) <- gsub(paste0("^.*", project_code, "_"), "", colnames(summary_data))
colnames(summary_data) <- gsub("\\.sort\\.bam$", "", colnames(summary_data))
# Transpose the summary data so that samples are rows, and categories are columns
summary_t <- as.data.frame(t(summary_data))
summary_t$Sample <- rownames(summary_t)
summary_t <- summary_t[summary_t$Sample != "Sample", ]
# Make a long df for plotting
summary_long <- pivot_longer(
  summary_t,
  cols = -Sample,
  names_to = "Category",
  values_to = "Count"
)
# Change counts to numeric and filter zeroes
summary_long$Count <- as.numeric(summary_long$Count)
summary_long_filtered <- summary_long %>%
  filter(!is.na(Count) & !is.na(Category)) %>%
  filter(Count > 0)
# Calculate the percent of assigned reads
assigned_percent <- summary_long_filtered %>%
  group_by(Sample) %>%
  mutate(Percent = ifelse(Category == "Assigned", Count / sum(Count) * 100, NA))
assigned_percent_filtered <- assigned_percent %>%
  filter(!is.na(Percent) & Category == "Assigned")
category_colors <- c("Assigned" = "#1F77B4",  # Dark blue for Assigned
                     "Unassigned_Ambiguity" = "#A6CEE3",
                     "Unassigned_Chimera" = "#ADD8E6",
                     "Unassigned_Duplicate" = "#B0E0E6",
                     "Unassigned_FragmentLength" = "#B0C4DE",
                     "Unassigned_MappingQuality" = "#C0D9E9",
                     "Unassigned_MultiMapping" = "#C8D8E6",
                     "Unassigned_NoFeatures" = "#D1E2F0",
                     "Unassigned_NonSplit" = "#D6E8F5",
                     "Unassigned_Overlapping_Length" = "#DAF0F9",
                     "Unassigned_Read_Type" = "#E0F5FB",
                     "Unassigned_Secondary" = "#E6FBFF",
                     "Unassigned_Singleton" = "#F0FEFF",
                     "Unassigned_Unmapped" = "lightgray")

# Calculate the percent of assigned reads in assigned_percent
assigned_percent <- summary_long_filtered %>%
  group_by(Sample) %>%
  mutate(Percent = ifelse(Category == "Assigned", Count / sum(Count) * 100, NA))

# Filter for rows where Percent is not NA and Category is "Assigned"
assigned_percent_filtered <- assigned_percent %>%
  filter(!is.na(Percent) & Category == "Assigned")

# Merge assigned_percent with summary_long_filtered
summary_long_filtered <- left_join(summary_long_filtered,
                                   assigned_percent %>% dplyr::select(Sample, Category, Percent),
                                   by = c("Sample", "Category"))

# Create a new column to indicate if assigned percent is less than 60
summary_long_filtered <- summary_long_filtered %>%
  mutate(RedOutline = ifelse(Category == "Assigned" & Percent < 60, "red", NA))

# Plot the data with conditional red outline for bars where Assigned % < 60
ggplot(summary_long_filtered, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_stack(),
           color = summary_long_filtered$RedOutline,
           size = 1) +  # Add size for the border
  scale_fill_manual(values = category_colors) +  # Apply custom fill colors
  geom_text(data = assigned_percent_filtered,
            aes(x = Sample, y = Count, label = paste0(round(Percent, 1), "%")),
            position = position_stack(vjust = 1.1),
            inherit.aes = FALSE, size = 4, color = "black") +
  labs(title = "Read Counts per Sample and Category",
       x = "Sample",
       y = "Read Count",
       fill = "Category") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(project_code, "mappingBar.pdf"), width = 10, height = 8)

###AX### extract columns based on if they are .bam data
col_names <- colnames(countdata)
print(col_names)
# Identify columns with .bam and without .bam
genedata_cols <- grep("\\.bam$", col_names, invert = TRUE)  # Columns without .bam
countdata_cols <- grep("\\.bam$", col_names)  # Columns with .bam
# Split the data into genedata and countdata
genedata <- countdata[, genedata_cols]
countdata <- countdata[, countdata_cols]

###AX### extract colnames and auto-rename
col_names_samples <- colnames(countdata)
# Get all column names (except the first which is likely a gene annotation or length info)
clean_col_names <- sub(paste0(".*", project_code, "_"), "", col_names_samples)
clean_col_names <- sub("\\.sort\\.bam$", "", clean_col_names)
print(clean_col_names)
colnames(countdata) <- clean_col_names
write.csv(countdata, paste0(project_code, "_CountData_Clean.csv"), row.names = TRUE)

###TCR### cal TPM and CPM
cpm <- t(t(countdata)/colSums(countdata))
avg_cpm <- data.frame(avg_cpm=rowMeans(cpm))

#-----TPM Calculation------
kb <- genedata$Length / 1000
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
avg_tpm <- data.frame(avg_tpm=rowMeans(tpm))

# write and reload
write.csv(avg_tpm, paste0(project_code, "_avg_tpm.csv"))
write.csv(avg_cpm, paste0(project_code, "_avg_cpm.csv"))
write.csv(tpm, paste0(project_code, "_tpm.csv"))
write.csv(cpm, paste0(project_code, "_cpm.csv"))

# heatmap of TPM
df <- read.csv(file = paste0(project_code, "_tpm.csv"), row.names = 1)
# set cutoff
df <- df[which(rowSums(df)>10),]
df <- na.omit(df)
J_map <- df
J_map <- J_map[rowSums(df[, -1] != 0) > 0, ]
pheatmap(J_map, scale = "row", show_rownames = FALSE, clustering_method = "ward.D")

# set TPM>1 filter
#tpmfilter <- read.csv('RP48_RNAseq_tpm.csv', sep= ',', header=T, row.names = 1)
#colnames(tpmfilter) <- c("preTran_1","preTran_2", "preTran_3", "IV_1","IV_2", "IV_3", "SubQ_1", "SubQ_2", "SubQ_3" )
#tpmfilter <- tpmfilter %>% filter((tpmfilter$preTran_1+tpmfilter$preTran_2+tpmfilter$preTran_3)>3)
#write.csv(tpmfilter,"RP48_RNAseq_tpmfilter.csv")


#FilterCount <- merge(countdata, tpmfilter, by = 0)
#FinalInput <- FilterCount[,1:10]
#colnames(FinalInput) <- c("ENSEMBL","preTran_1","preTran_2", "preTran_3", "IV_1","IV_2", "IV_3", "SubQ_1", "SubQ_2", "SubQ_3" )

#write.csv(FinalInput,"RP48_RNAseq_count_afterfilterTPM.csv", row.names=FALSE)

##change col names and reload
#Data <- read.csv('RP48_RNAseq_count_afterfilterTPM.csv', sep= ',', header=T, row.names = 1)

#remove NA,
Data <- na.omit(countdata)

#chose Data set
Mat <- Data

###AX### generate metadata file from counts data
sample_ids <- colnames(Data)

# Extract sample type (based on sample name before the number)
prefixes <- sub("_[0-9]+$", "", sample_ids)

# Prompt user to name each sample type based on unique samples and generate the sample column
unique_prefixes <- unique(prefixes)
sample_names <- character(length(unique_prefixes))
for (i in seq_along(unique_prefixes)) {
  sample_names[i] <- readline(prompt = paste0("→ What is the sample name for '", unique_prefixes[i], "'? "))
}
names(sample_names) <- unique_prefixes
sample_column <- sample_names[prefixes]

# Assign replicate numbers within each sample group
replicate_column <- ave(sample_ids, sample_column, FUN = function(x) seq_along(x))

# Create the metadata dataframe
Metadata <- data.frame(
  sample_id = sample_ids,
  sample = sample_column,
  Replicate = replicate_column,
  stringsAsFactors = FALSE
)
print(Metadata)
write.csv(Metadata, paste0(project_code, "_meta.csv"), row.names = FALSE)

###TCR### generate dds matrix
Metadata$sample <- as.factor(Metadata$sample)
dds <-DESeqDataSetFromMatrix(countData=Mat,
                             colData=Metadata,
                             design=~sample)

###TCR### filter dds >1
dds <- dds[rowSums(counts(dds))>1]

###TCR### sample clustering
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
hc <- hclust(sampleDists, method = "ward.D2")
plot(hc, hang = -1)

## Run DESeq2
ds2 <- DESeq(dds, betaPrior = F) # run DESeq using experiment file
resultsNames(ds2)

###AX### filter metadata
# Specify the samples you want to exclude
exclude_samples <- c("mCh_505", "GFP_505")
# Subset the count data to exclude the unwanted samples
Mat_filtered <- Mat[, !colnames(Mat) %in% exclude_samples]
# Subset the metadata to exclude the corresponding rows
Metadata_filtered <- Metadata[!Metadata$sample_id %in% exclude_samples, ]
# Now create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = Mat_filtered,
                              colData = Metadata_filtered,
                              design = ~sample)
dds <- dds[rowSums(counts(dds)) > 1]
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
hc <- hclust(sampleDists, method = "ward.D2")
plot(hc, hang = -1)
ds2 <- DESeq(dds, betaPrior = F) # run DESeq using experiment file
resultsNames(ds2)

###TCR### Normalized counts
normalized_counts <- as.data.frame(counts(ds2, normalized=TRUE))
write.csv(normalized_counts, 'MobaV_normalizedcounts.csv')

###TCR### PCA plot
vsdB <- varianceStabilizingTransformation(ds2)
p <- plotPCA(vsdB, intgroup = "sample")
# Add points and labels using rownames from colData
p +
  geom_point(shape = 21, size = 3, stroke = 1, color = "black", aes(fill = group)) +
  geom_text_repel(aes(label = rownames(vsdB@colData)), color = "black", size = 3, vjust = -0.5)  # Adjust label position
ggsave("PCA_plot_f.pdf", plot = p+
         geom_point(shape = 21, size = 3, stroke = 1, color = "black", aes(fill = group))  +
         geom_text_repel(aes(label = rownames(vsdB@colData)), color = "black", size = 3, vjust = -0.5) , width = 8, height = 8, units = "in")

###TCR### Corplot
CorMatrix <- cor(normalized_counts)
heatmaply_cor(CorMatrix, limits = c(min(CorMatrix),1))

###AX### generate fold change results from user input
# Get the levels of the 'sample' factor
levels_sample <- levels(ds2$sample)

# Display numbered options
cat("Available sample groups:\n", paste(sprintf("%d: %s", seq_along(levels_sample), levels_sample), collapse = "\n"), "\n")

# Prompt user to enter two numbers
user_input <- readline(prompt = "Enter the numbers for the test and reference group in order, separated by a comma: ")

# Split and parse the input
selected_indices <- as.numeric(strsplit(user_input, ",")[[1]])

# Validate input
if (length(selected_indices) == 2 &&
    all(!is.na(selected_indices)) &&
    all(selected_indices %in% seq_along(levels_sample))) {

  test <- levels_sample[selected_indices[1]]
  reference <- levels_sample[selected_indices[2]]

  # Confirm selection
  cat("\nYou selected:\n")
  cat("Test:", test, "\n")
  cat("Reference:", reference, "\n")

  # Run DESeq2 results
  ds2_results <- results(ds2, contrast = c("sample", test, reference))

} else {
  cat("Invalid input. Please enter two valid numbers corresponding to the sample groups.\n")
}

## output
write.csv(ds2_results, paste0(test, "_vs_", reference, "_result.csv"))

###AX### generate ENSEMBL results based on user input
# List all result files in the current directory
result_files <- list.files(pattern = "result.*\\.csv$")

# Display choices to user
cat("Available result files:\n", paste(sprintf("%d: %s", seq_along(result_files), result_files), collapse = "\n"), "\n")

# Prompt user to select a file by number
selection <- as.integer(readline(prompt = "Select a result file by number: "))
selected_file <- result_files[selection]

# Read selected file
res <- read.csv(selected_file, row.names = NULL)

res$symbol <- mapIds(org.Mm.eg.db,
                     keys = res[, 1],
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

res$entrez <- mapIds(org.Mm.eg.db,
                     keys = res[, 1],
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

res <- na.omit(res)

# Save updated results
write.csv(res, sub("\\.csv$", "_final.csv", selected_file))

###AX### final data input for volcano plot
final_files <- list.files(pattern = "final.*\\.csv$")
if (length(final_files) == 0) {
  stop("No files ending in 'final.csv' found.")
}
cat("Available final result files:\n",
    paste(sprintf("%d: %s", seq_along(final_files), final_files), collapse = "\n"), "\n")

# Prompt user to select a file
file_index <- as.integer(readline(prompt = "Enter the number corresponding to the file: "))
selected_file <- final_files[file_index]

# Read the selected file
voldata <- read.csv(file = selected_file, header = TRUE, row.names = 1)

###TCR## volplot
voldata$label=ifelse(((voldata$padj < 0.001)&(abs(voldata$log2FoldChange) > 2.5)&(voldata$baseMean >10)),'Yes','No')
sig_data <- voldata %>% filter(label == "Yes")
# Find top 20 by lowest padj
top_padj <- sig_data %>%
  arrange(padj) %>%
  slice_head(n = 10)

# Find top 20 by largest |log2FoldChange|
top_lfc <- sig_data %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)

# Combine (remove duplicates if overlap)
top_labels <- bind_rows(top_padj, top_lfc) %>%
  distinct(symbol, .keep_all = TRUE)
# Find the range of the x and y values
x_range <- range(voldata$log2FoldChange, na.rm = TRUE)
x_range <- c(-4.5, 4)
y_range <- range(-1 * log10(voldata$padj), na.rm = TRUE)
y_range <- c(0, 50)
plot_title <- paste("Volcano Plot", gsub("_", " ", sub("_result_final$", "", tools::file_path_sans_ext(basename(selected_file)))))

# Create the plot with dynamic coordinate limits
p <- ggplot(data = voldata, aes(x = log2FoldChange, y = -1 * log10(padj))) +
  geom_point(aes(color = label)) +
  scale_color_manual(values = c("#d2dae2", "#b50719")) +
  geom_text_repel(data = top_labels, aes(x = log2FoldChange, y = -1 * log10(padj), label = symbol), size = 3,
                  max.overlaps = 20) +
  labs(title = plot_title, x = expression(log[2](FC)), y = expression(-log[10](padj)), subtitle = "log2FoldChange>2.5; padj<0.001; MeanConts >10") +
  geom_hline(yintercept = 1.3, linetype = 4) +  # set line adj = 0.05
  geom_vline(xintercept = c(-1, 1), linetype = 4) +
  coord_cartesian(ylim = c(min(y_range), max(y_range)), xlim = c(min(x_range), max(x_range))) +
  theme_bw() + theme(panel.grid = element_blank())
p <- p + theme(legend.position = "none")
print(p)

# Count upregulated and downregulated genes
up_count <- sum(voldata$log2FoldChange >= 1 & voldata$padj <= 0.05, na.rm = TRUE)
down_count <- sum(voldata$log2FoldChange <= -1 & voldata$padj <= 0.05, na.rm = TRUE)

# Add annotation to the plot
p <- p +
  annotate("text", x = max(x_range) - 0.5, y = max(y_range) - 5,
           label = paste0("Up: ", up_count), color = "red", hjust = 1, size = 4) +
  annotate("text", x = min(x_range) + 0.5, y = max(y_range) - 5,
           label = paste0("Down: ", down_count), color = "blue", hjust = 0, size = 4)

print(p)

# Export downregulated genes for HOMER
downregulated_genes <- voldata %>%
  dplyr::filter(padj < 0.001, log2FoldChange < -2.5, baseMean > 10) %>%
  dplyr::pull(symbol)

# Generate list of background genes
background_genes <- voldata %>%
  dplyr::filter(baseMean > 10) %>%
  dplyr::pull(symbol)

# Save to a text file for HOMER
writeLines(downregulated_genes, "cdx2_downregulated_genes.txt")
writeLines(background_genes, "background_genes.txt")

##### Replot volcano plot with HOMER targets #####
target_hits <- read.table("cdx2_promoter_genes.csv", header = TRUE, sep = ",")
voldata$highlight <- voldata$symbol %in% merged_df$Gene.Name

# Plot
p <- ggplot(data = voldata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = label), alpha = 0.8) +   # color by significance
  scale_color_manual(values = c("No" = "#d2dae2", "Yes" = "#b50719")) +

  # Add an outline or different shape for highlight genes
  geom_point(data = filter(voldata, highlight),
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, color = "steelblue", size = 2, stroke = 1) +  # black circle outline

  # Label only highlighted genes
  geom_text_repel(data = filter(voldata, highlight),
                  aes(label = symbol),
                  size = 3,
                  max.overlaps = 20) +

  labs(title = plot_title,
       x = expression(log[2](FC)),
       y = expression(-log[10](padj)),
       subtitle = "log2FoldChange>2; padj<0.005; MeanConts >10 colored red; Highlighted genes labeled with outline") +

  geom_hline(yintercept = 1.3, linetype = 4) +
  geom_vline(xintercept = c(-1, 1), linetype = 4) +

  coord_cartesian(ylim = c(min(y_range), max(y_range)),
                  xlim = c(min(x_range), max(x_range))) +

  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none")

print(p)

# Add in non-DE genes and DE gene list
# findMotifs.pl cdx2_downregulated_genes.txt mouse homer_output/ -bg background_genes.txt -start -1000 -end 100 -p 4
# List of genes and promoter and p-value of if its related

# Save the plot using the generated filename
ggsave(filename = file.path(getwd(), paste0("VolcPlot", strsplit(selected_file, "_")[[1]][1], "filt.pdf")), width = 11, height = 11)

##### AX generate safe vs KO count data for heatmap #####
safeKOcountdata <- countdata %>%
  mutate(GeneName = mapIds(org.Mm.eg.db,
                           keys = rownames(countdata),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")) %>%
  filter(!is.na(GeneName)) %>%  # Remove rows where GeneName is NA
  dplyr::select(GeneName, everything()) %>%  # Move GeneName to the first column
  as.data.frame()

# Remove row names
rownames(safeKOcountdata) <- NULL
# Remove rows where all gene expression is 0
safeKOcountdata <- safeKOcountdata[rowSums(safeKOcountdata[, -1]) > 0, ]
# Sum up rows with same gene name
safeKOcountdata <- safeKOcountdata %>%
  group_by(GeneName) %>%
  summarise(across(everything(), sum, na.rm = TRUE))  # Sum expression values across columns
safeKOcountdata <- as.data.frame(safeKOcountdata)
# Ensure GeneName is set as row names and remove it from the dataframe
rownames(safeKOcountdata) <- safeKOcountdata$GeneName  # Assign GeneName as row names
safeKOcountdataR <- safeKOcountdata[, -1]
#gene_variances <- apply(safeKOcountdata, 1, var)

# Get the indices of the top 50 most variable genes
#top50_genes <- order(gene_variances, decreasing = TRUE)[1:50]

# Subset the data for the top 50 most variable genes
#top50_matrix <- safeKOcountdata[top50_genes, ]

# Convert to a numeric matrix
#top50_matrix <- as.matrix(top50_matrix)
# Convert to a numeric matrix
safeKO_matrix <- as.matrix(safeKOcountdataR)

# Ensure all values are numeric
safeKO_matrix <- apply(safeKO_matrix, 2, as.numeric)
rownames(safeKO_matrix) <- rownames(safeKOcountdataR)  # Retain gene names
# Log 2 transform
safeKO_matrixlog2 <- log2(safeKO_matrix)
safeKO_matrixlog2[safeKO_matrixlog2 == -Inf] <- -1

# Generate heatmap with gene names
pheatmap(safeKO_matrixlog2,
         color = colorRampPalette(c("magenta", "black", "yellow"))(50),
         scale = "row",
         show_rownames = F,  # Display gene names
         clustering_method = "ward.D"
         )

###231206 TCR## find shared genes

LMet <-read.csv(file = "LiverMet_vs_SubQ_Final.csv",header = TRUE, row.names =1)
LMet <- LMet %>% filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)

SMet <-read.csv(file = "SubQLiver_vs_SubQ_Final.csv",header = TRUE, row.names =1)
SMet <- SMet %>% filter((log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)

Expression <- read.csv(file = "RP48_LiverMet_tpm.csv",header = TRUE, row.names =1)
Exp <- Expression %>% filter(rowMeans(select(., CellCulture_1, CellCulture_2, CellCulture_3)) > 1)


merged_Met <- merge(LMet, SMet, by = "symbol")
Exp_Met <- merge(merged_Met, Exp, by.x = "ENSEMEL", by.y = "ENSEMEL.1")



p <- plot(merged_Met$log2FoldChange.x, merged_Met$log2FoldChange.y, main = "Scatter Plot", xlab = "FoldChange in LiverMet", ylab = "FoldChange in SubQMet", pch = 16, col = "blue")
print(p)

write.csv(merged_Met, 'merged_Met_genes.csv')


###TCR##GSEA analysis############################################
res <- read.csv("WT_vs_CrebbpKO_result_final_filtered.csv", row.names = NULL)

#select 3 cols
gene_df <- res %>%
  dplyr::select(symbol, log2FoldChange, entrez) %>%
  filter(!is.na(entrez)) %>%
  distinct(entrez, .keep_all = TRUE)

##logFC
geneList <- gene_df$log2FoldChange
##gene name
names(geneList) = gene_df$entrez
##sort
geneList = sort(geneList, decreasing = TRUE)

## GSEA
gseaKEGG <- gseKEGG(geneList     = geneList,
                    organism     = 'mmu',
                    minGSSize    = 10,
                    pvalueCutoff = 0.05,
                    verbose      = F)
##Plot pathway
dotplot(gseaKEGG)
dotplot(gseaKEGG,split=".sign")+facet_grid(~.sign)
ggsave(filename = paste0(getwd(), "/", "gseaKEGG.pdf"), width = 11, height = 11)

##GSEA pathway report
gseaKEGG_results <- gseaKEGG@result
write.csv(gseaKEGG_results, 'crebbpKO_vs_wildtype_GSEA.csv')

##GSEA enrichment plot
library(enrichplot)
pathway.id = "hsa04261"
gseaplot2(gseaKEGG,
          color = "blue",
          geneSetID = pathway.id,
          pvalue_table = T)
ggsave(filename = paste0(getwd(), "/", "gseaEnrichplot.pdf"), width = 11, height = 11)

#### TCR### GSEA analysis with file selection and dynamic naming AX ####
set.seed(1234)
# Select result file
final_files <- list.files(pattern = "_result_final\\.csv$")
if (length(final_files) == 0) {
  stop("No files ending in '_result_final.csv' found.")
}
cat("Available final result files:\n",
    paste(sprintf("%d: %s", seq_along(final_files), final_files), collapse = "\n"), "\n")

file_index <- as.integer(readline(prompt = "Enter the number corresponding to the file: "))
selected_file <- final_files[file_index]

# Extract comparison name from filename
comparison_name <- sub("_result_final\\.csv$", "", selected_file)
cat("Analysis for:", comparison_name, "\n")

# Read the selected file
res <- read.csv(selected_file, row.names = NULL)

# Prepare gene list
gene_df <- res %>%
  dplyr::select(symbol, log2FoldChange, entrez) %>%
  filter(!is.na(entrez) & !is.na(symbol)) %>%
  distinct(entrez, .keep_all = TRUE)

geneList <- gene_df$log2FoldChange
names(geneList) = gene_df$entrez
geneList = sort(geneList, decreasing = TRUE)

# Gene list with symbols for plotting
geneList_symbols <- gene_df$log2FoldChange
names(geneList_symbols) = gene_df$symbol
geneList_symbols = sort(geneList_symbols, decreasing = TRUE)

## GSEA KEGG
gseaKEGG <- gseKEGG(geneList     = geneList,
                    organism     = 'mmu',
                    minGSSize    = 10,
                    pvalueCutoff = 0.05,
                    verbose      = F)

# Plot pathway
dotplot(gseaKEGG, split=".sign") + facet_grid(~.sign)
ggsave(filename = paste0(comparison_name, "_gseaKEGG.pdf"), width = 11, height = 11)

# Save GSEA pathway report
gseaKEGG_results <- gseaKEGG@result
write.csv(gseaKEGG_results, paste0(comparison_name, '_GSEA_KEGG.csv'))

## GSE GO analysis
res_go <- gseGO(geneList,
                ont = "BP",
                OrgDb = org.Mm.eg.db,
                keyType = "ENTREZID",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH")

dotplot(res_go, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave(filename = paste0(comparison_name, "_gseGO_split.pdf"), width = 11, height = 11)

#### GSEA concept network plots ####
gseaKEGG <- setReadable(gseaKEGG, 'org.Mm.eg.db', 'ENTREZID')

cnetplot(gseaKEGG,
         categorySize = "pvalue",
         foldChange = geneList_symbols,
         showCategory = 5,
         colorEdge = TRUE,
         node_label = "all")
ggsave(filename = paste0(comparison_name, "_KEGG_cnetplot2.pdf"),
       width = 14, height = 12)

res_go <- setReadable(res_go, 'org.Mm.eg.db', 'ENTREZID')

cnetplot(res_go,
         categorySize = "pvalue",
         foldChange = geneList_symbols,
         showCategory = 5,
         colorEdge = TRUE,
         node_label = "all")
ggsave(filename = paste0(comparison_name, "_GO_cnetplot2.pdf"),
       width = 14, height = 12)

#### GSEA enrichment maps ####
library(enrichplot)
gseaKEGG_pairwise <- pairwise_termsim(gseaKEGG)
emapplot(gseaKEGG_pairwise, showCategory = 30)
ggsave(filename = paste0(comparison_name, "_KEGG_emapplot.pdf"),
       width = 12, height = 12)

res_go_pairwise <- pairwise_termsim(res_go)
emapplot(res_go_pairwise, showCategory = 30)
ggsave(filename = paste0(comparison_name, "_GO_emapplot2.pdf"),
       width = 12, height = 12)

##GSEA enrichment plot
pathway.id = "hsa04261"
gseaplot2(gseaKEGG,
          color = "blue",
          geneSetID = pathway.id,
          pvalue_table = T)
ggsave(filename = paste0(getwd(), "/", "gseaEnrichplot.pdf"), width = 11, height = 11)

#GSEA pathway mapping
library(pathview)
pathway.id = "hsa04261"
pv.out <- pathview(gene.data  = geneList,
                   pathway.id = pathway.id,
                   species    = "hsa")



##GSE GO analysis
res <- gseGO(
  geneList,
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

dotplot(res)
ggsave(filename = paste0(getwd(), "/", "gseGO.pdf"), width = 11, height = 11)

dotplot(
  res,
  showCategory=10,
  split=".sign") + facet_grid(.~.sign)
ggsave(filename = paste0(getwd(), "/", "gseGOsplit.pdf"), width = 11, height = 11)

#############################
##TCR## plot VennDiagram



library(VennDiagram)

dat <- read.table('LiverMet_DownGenes.csv', header = TRUE, sep = ',')

venn_list <- list(LiverMet = dat$LiverMet, SubQMet = dat$SubQLiver)

venn.diagram(venn_list,
             resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),
             fill=c("red","blue"),
             main="Down-regulated genes in LiverMet and SubQMet",
             cat.col = rep('black', 2),

             col = 'black', cex = 1.5, fontfamily = 'serif',

             cat.cex = 1.5, cat.fontfamily = 'serif',
             filename = "VennDiagram_Downgene.tif")


inter <- get.venn.partitions(venn_list)
###TCR### QC of Fold change in known genes, such as FN1 (Ensembl:ENSG00000115414)
plotCounts(ds2, gene = "ENSMUSG00000026180", intgroup=c("sample"))

###TCR### Fold change, clean those with small counts, and MAplot
contrast <- c("sample", "Bone_Chip_EC", "MSC_Chip_EC")
dd1 <- results(ds2, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-4,4))

resLFC <- lfcShrink(ds2, coef="sample_MSC_Chip_EC_vs_Bone_Chip_EC", type="apeglm")
plotMA(resLFC, ylim=c(-2,2), main="MSC_vs_Bone", alpha = 0.05)

##### Heat map #####
# Define samples to exclude
exclude_samples <- c("mCh_505","GFP_505")
countdata <- read.table("/Users/tanglab/Library/CloudStorage/Box-Box/AX_data_analysis/AX_RNASeq_Data/MobaV/MobaV_normalizedcounts.csv", sep=",", header = T )
rownames(countdata) <- countdata$X
countdata <- countdata %>% dplyr::select(-X)
# Process count data (your existing code with exclusion)
safeKOcountdata <- countdata %>%
  mutate(GeneName = mapIds(org.Mm.eg.db,
                           keys = rownames(countdata),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")) %>%
  filter(!is.na(GeneName)) %>%
  dplyr::select(GeneName, everything()) %>%
  as.data.frame()

# Remove row names
rownames(safeKOcountdata) <- NULL

# Exclude specified samples
safeKOcountdata <- safeKOcountdata %>%
  dplyr::select(-any_of(exclude_samples))

# Remove rows where all gene expression is 0
safeKOcountdata <- safeKOcountdata[rowSums(safeKOcountdata[, -1]) > 0, ]

# Sum up rows with same gene name
safeKOcountdata <- safeKOcountdata %>%
  group_by(GeneName) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

safeKOcountdata <- as.data.frame(safeKOcountdata)

# Set GeneName as row names
rownames(safeKOcountdata) <- safeKOcountdata$GeneName
safeKOcountdataR <- safeKOcountdata[, -1]

# Convert to numeric matrix
safeKO_matrix <- as.matrix(safeKOcountdataR)
safeKO_matrix <- apply(safeKO_matrix, 2, as.numeric)
rownames(safeKO_matrix) <- rownames(safeKOcountdataR)

# Log2 transform with pseudocount
safeKO_matrixlog2 <- log2(safeKO_matrix + 1)

# Identify sample groups based on column names
sample_names <- colnames(safeKO_matrixlog2)
mCh_samples <- grep("mCh", sample_names, value = TRUE)
GFP_samples <- grep("GFP", sample_names, value = TRUE)

# Reorder columns: mCh (crebbpKO) first, then GFP (Wild-type)
column_order <- c(mCh_samples, GFP_samples)
safeKO_matrixlog2 <- safeKO_matrixlog2[, column_order]

# Calculate mean expression for each group
mCh_mean <- rowMeans(safeKO_matrixlog2[, mCh_samples, drop = FALSE])
GFP_mean <- rowMeans(safeKO_matrixlog2[, GFP_samples, drop = FALSE])

# Calculate fold change (mCh vs GFP)
fold_change <- mCh_mean - GFP_mean

# Filter genes with meaningful expression
mean_expression <- rowMeans(safeKO_matrixlog2)
expressed_genes <- names(mean_expression)[mean_expression > 1]
fold_change_filtered <- fold_change[expressed_genes]

# Get top 20 up-regulated and down-regulated genes
top20_up <- names(sort(fold_change_filtered, decreasing = TRUE)[1:20])
top20_down <- names(sort(fold_change_filtered, decreasing = FALSE)[1:20])
top_genes <- c(top20_up, top20_down)

# Create subset matrix with top regulated genes
heatmap_matrix <- safeKO_matrixlog2[top_genes, ]

# Create annotation for samples
sample_annotation <- data.frame(
  Group = ifelse(colnames(heatmap_matrix) %in% mCh_samples, "crebbpKO", "Wild-type"),
  row.names = colnames(heatmap_matrix)
)

# Define colors for annotations
ann_colors <- list(
  Group = c("crebbpKO" = "#FF6B6B", "Wild-type" = "#4ECDC4")
)

# Generate heatmap with clustering disabled for columns
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("magenta", "black", "yellow"))(100),
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,  # Disable column clustering to maintain order
         clustering_method = "ward.D2",
         annotation_col = sample_annotation,
         annotation_colors = ann_colors,
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Top 20 Up- and Down-regulated Genes\n(crebbpKO vs Wild-type)",
         border_color = NA)

# Top 500 genes heatmap
top600_up <- names(sort(fold_change_filtered, decreasing = TRUE)[1:500])
top600_down <- names(sort(fold_change_filtered, decreasing = FALSE)[1:500])
selected_genes <- c(top600_up, top600_down)

# Subset expression matrix
heatmap_matrix <- safeKO_matrixlog2[selected_genes, ]

# Z-score normalization
z <- t(scale(t(heatmap_matrix)))

# Sample annotation
sample_annotation <- data.frame(
  Group = ifelse(colnames(z) %in% mCh_samples, "crebbpKO", "Wild-type"),
  row.names = colnames(z)
)

# Color annotations
ann_colors <- list(
  Group = c("crebbpKO" = "#FF6B6B", "Wild-type" = "#4ECDC4")
)

# Draw heatmap with column clustering
heat <- pheatmap(z,
                 color = colorRampPalette(c("magenta", "black", "yellow"))(100),
                 scale = "none",
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 cluster_cols = FALSE,  # Disable column clustering
                 clustering_method = "ward.D2",
                 annotation_col = sample_annotation,
                 annotation_colors = ann_colors,
                 fontsize_row = 6,
                 fontsize_col = 10,
                 main = "Top 500 Up & Down Genes",
                 border_color = NA)

# Add flags to top 20 genes
flagged_heatmap <- add.flag(heat,
                            kept.labels = c("Cdx2"),
                            repel.degree = 0)
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels,
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)),
                    length.out = sum(d.select)),
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.15, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4,
                                     l = 4
  )

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid.newpage()
  grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

#### Heatmap revision ####
# reload data and plot heat map
countdata <- read.table("/Users/tanglab/Library/CloudStorage/Box-Box/AX_data_analysis/AX_RNASeq_Data/MobaV/MobaV_normalizedcounts.csv", sep=",", header = T)

# Set rownames and remove X column
rownames(countdata) <- countdata$X
countdata <- countdata %>% dplyr::select(-X)

# Map Ensembl IDs to gene symbols
countdata_with_symbols <- countdata %>%
  mutate(GeneName = mapIds(org.Mm.eg.db,
                           keys = rownames(countdata),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")) %>%
  filter(!is.na(GeneName))

# Sum rows with the same gene name
countdata_with_symbols <- countdata_with_symbols %>%
  group_by(GeneName) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  as.data.frame()

# Now set gene symbols as rownames
rownames(countdata_with_symbols) <- countdata_with_symbols$GeneName
countdata <- countdata_with_symbols %>% dplyr::select(-GeneName)

##### Heatmap AX - FOLD CHANGE APPROACH #####
## set cutoff
countdata <- countdata[which(rowSums(countdata) > 10), ]
countdata <- na.omit(countdata)
row_sums <- rowSums(countdata)

# Get the 10 genes with lowest rowSums
lowest_genes <- sort(row_sums)[1:10]
print(lowest_genes)

# Log2 transform
countdata_log2 <- log2(countdata + 1)

# Identify sample groups
sample_names <- colnames(countdata_log2)
mCh_samples <- grep("mCh", sample_names, value = TRUE)
GFP_samples <- grep("GFP", sample_names, value = TRUE)

# Calculate mean expression for each group
mCh_mean <- rowMeans(countdata_log2[, mCh_samples, drop = FALSE])
GFP_mean <- rowMeans(countdata_log2[, GFP_samples, drop = FALSE])

# Calculate fold change (mCh vs GFP)
fold_change <- mCh_mean - GFP_mean

# Filter genes with meaningful expression
mean_expression <- rowMeans(countdata_log2)
expressed_genes <- names(mean_expression)[mean_expression > 1]
fold_change_filtered <- fold_change[expressed_genes]

# Get top 1000 up-regulated and 1000 down-regulated genes
top_up <- names(sort(fold_change_filtered, decreasing = TRUE)[1:500])
top_down <- names(sort(fold_change_filtered, decreasing = FALSE)[1:500])
selected_genes <- c(top_up, top_down)

# Create matrix for heatmap
mat_heatmap <- as.matrix(countdata_log2[selected_genes, ])

# Reorder columns: mCh first, then GFP
column_order <- c(mCh_samples, GFP_samples)
mat_heatmap <- mat_heatmap[, column_order]

# Label only Cdx2
rownames(mat_heatmap) <- ifelse(grepl("Cdx2", rownames(mat_heatmap)),
                                rownames(mat_heatmap), "")

# Z-score normalization
z <- t(scale(t(mat_heatmap)))

# Sample annotation
sample_annotation <- data.frame(
  Group = ifelse(colnames(z) %in% mCh_samples, "crebbpKO", "Wild-type"),
  row.names = colnames(z)
)

# Color annotations
ann_colors <- list(
  Group = c("crebbpKO" = "#FF6B6B", "Wild-type" = "#4ECDC4")
)

# Generate heatmap
heat_lab <- pheatmap(z,
                     color = colorRampPalette(c("magenta", "black", "yellow"))(50),
                     legend_breaks = c(-2, -1, 0, 1, 2, max(z)),
                     main = "Top 500 Up & Down Regulated Genes",
                     legend_labels = c("-2", "-1", "0", "1", "2", "z-score\n"),
                     show_rownames = TRUE,
                     cluster_cols = FALSE,  # Keep mCh and GFP grouped
                     clustering_method = "ward.D",
                     annotation_col = sample_annotation,
                     annotation_colors = ann_colors,
                     fontsize_row = 6)

add.flag(heat_lab,
         kept.labels = "Cdx2",
         repel.degree = 0)
