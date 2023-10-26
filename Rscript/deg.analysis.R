# Install and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

#######################################################
## Copy First Half Before Plotting for Gene Ontology ##
#######################################################

# Load data
# file <- "C:/biomedical_data_mining/teamproject/BS6202.count (1).mx"

setwd("C:/Irfan/Personal/Masters/BS6202/Project")
file <- "BS6202.count.mx"
data <- read.table(file, header = TRUE, row.names = 1)
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
data <- data+1

# Define sample conditions
# Create DESeq2 object
data <- data[, c('hoxd9.1',	'hoxd9.2',	'hoxd9.3', 'vector.1',	'vector.2',	'vector.3')]

conditions <- as.factor(c('HOXD9', 'HOXD9', 'HOXD9','Vector', 'Vector', 'Vector'))
# conditions <- as.factor(c('Vector', 'Vector', 'Vector', 'HOXD9', 'HOXD9', 'HOXD9'))
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=DataFrame(conditions), 
                              design=~conditions)

# Perform differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res)
sorted_genes <- res[order(res$padj),]
write.csv(sorted_genes, "DEG_genelist_P-adj.csv")
# write.csv(sorted_genes, "C:/biomedical_data_mining/teamproject/DEG_genelist_results.csv")

gene_type_file <- 'gene_info.csv'
gene_type <- read.table(gene_type_file, header = TRUE, row.names = 1, sep=',')

library(dplyr)
library(tibble)

sorted_genes <- as.data.frame(sorted_genes) %>%
  mutate(gene_id = rownames(.))

gene_type <- gene_type %>%
  mutate(gene_id = rownames(.))

merged_df <- inner_join(sorted_genes, gene_type, by = "gene_id")
merged_df$flippedlog2FoldChange <- -(merged_df$log2FoldChange)
merged_df$FoldChange <- 2^merged_df$flippedlog2FoldChange
merged_df$upregulated <- ifelse(merged_df$flippedlog2FoldChange > 0, 1, 0)

write.csv(merged_df, "DEG_genelist_P-adj.gene_info.csv")

merged_df <- read.csv("DEG_genelist_P-adj.gene_info.csv", sep=",", row.names=1, header=1)

significant_genes <- subset(merged_df, padj <= 0.05)
significant_genes <- subset(significant_genes, abs(log2FoldChange) >.75)

#Drop HOXD9
significant_genes <- subset(significant_genes, rownames(significant_genes) != "ENSG00000128709.11")

selected_rows <- significant_genes[significant_genes$gene_type == "lincRNA" & significant_genes$upregulated
                                   == 0, ]
selected_genes <- selected_rows$gene_name
print (selected_genes)

unique_counts1 <- table(significant_genes$gene_type)
unique_counts2 <- table(significant_genes$upregulated)

# Create a new column that represents the pairs
df <- paste(significant_genes$gene_type, "-", significant_genes$upregulated)

# Calculate unique counts for the pairs
unique_counts_pairs <- table(df)
print(unique_counts_pairs)

compare_genes <- subset(merged_df, padj <= 0.05)
compare_genes <- subset(compare_genes, abs(log2FoldChange) >.27)

li_file <- 'Li.TopGenes.fixed.csv'
li_data <- read.csv(li_file, sep=",", header=T, row.names=1)
selected_li_rows <- li_data[li_data$Gene_Type == "lincRNA" & li_data$Upregulated
 == 0, ]
selected_li_genes <- selected_li_rows$Gene_Name


selected_our_rows <- compare_genes[compare_genes$gene_type == "lincRNA" & compare_genes$upregulated
                            == 0, ]
selected_our_rows <- selected_our_rows$gene_name
overlap <- intersect(selected_our_rows, selected_li_genes)
overlap_length <- length(overlap)
print (overlap)
print (overlap_length)


#######################################################
## Copy First Half Before Plotting for Gene Ontology ##
#######################################################

count_data <- read.csv(file, sep="\t", row.names=1)
deg_data <- significant_genes

# Prepare the data for volcano plot
volcano_data <- as.data.frame(deg_data)

# Set the number of top genes you want to label
num_genes_label <- 12

# Select top genes based on adjusted p-value
top_genes <- volcano_data[order(volcano_data$padj), ][1:num_genes_label, ]

library(ggplot2)
library(pheatmap)
threshold_padj <- 0.05
# Create the enhanced volcano plot with labels
base_plot <- ggplot(merged_df, aes(x=-log2FoldChange, y=-log10(padj), color="black")) + 
  geom_point(alpha=0.6) +
  geom_text(data=top_genes, aes(label=top_genes$gene_name), vjust=-1, hjust=0.5, size=3) +
  labs(title="Enhanced Volcano Plot", x="Log2 Fold Change", y="-Log10(FDR)") +
  geom_hline(yintercept=-log10(threshold_padj), color="red", linetype="dashed") +
  geom_vline(xintercept=.75, color="red", linetype="dashed") +
  geom_vline(xintercept=-.75, color="red", linetype="dashed") +
  scale_color_identity()

x_limits <- c(-5, 5)  # Set your desired limits here

final_plot <- base_plot +
  geom_point(data=volcano_data, aes(x=-log2FoldChange, y=-log10(padj), color="red"), alpha=0.6) +
  scale_color_identity() +
  coord_cartesian(xlim = x_limits)  # Set the x-axis limits

print(final_plot)

# Heatmap of Top 50 Most Significant DEGs
row.names(deg_data) <- deg_data$gene_id
# top_genes <- rownames(deg_data[order(deg_data$padj)[1:50],]) #Top 50
# top_genes <- deg_data[deg_data$gene_type == "protein_coding", ] #Protein Coding
# top_genes <- rownames(top_genes[order(top_genes$padj)[1:50],])

top_genes <- deg_data[deg_data$gene_type == "lincRNA", ] #lincRNA
top_genes <- rownames(top_genes)
heatmap_data <- count_data[top_genes,]

heatmap_data <- as.data.frame(heatmap_data) %>%
  mutate(gene_id = rownames(.))

heatmap_data_merged <- inner_join(deg_data, heatmap_data, by = "gene_id")
row.names(heatmap_data_merged) <- heatmap_data_merged$gene_name

heatmap_data_merged <- heatmap_data_merged %>% select(vector.1, vector.2, vector.3, hoxd9.1, hoxd9.2, hoxd9.3)


# Drawing heatmap with adjusted font sizes
pheatmap(heatmap_data_merged, 
         show_rownames=T, 
         show_colnames=T, 
         scale="row", 
         fontsize=8,       
         fontsize_row=6) 

# MA Plot
DESeq2::plotMA(res, ylim=c(-2,2))


# GO Enrichment
# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

# Load the necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Assuming you have already run the previous DESeq2 code and obtained 'significant_genes'
# Extract gene symbols from DESeq2 results
gene_list <- rownames(significant_genes)

# Perform GO enrichment analysis
ego <- enrichGO(gene         = gene_list,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP", 
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05)

# Visualization
dotplot(ego, showCategory=20)

