
# Load necessary libraries
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

setwd("C:/Irfan/Personal/Masters/BS6202/Project")

file <- 'test_results_C2_intron_retention.tsv'
data <- read.table(file, header = TRUE, row.names = 1)
significant <- subset(data, p_val_adj <= 0.05)
genes_to_test <- significant$gene_id
print(length(genes_to_test))

library(stringr)
newList <- vector()
for (i in genes_to_test) {
  newList <- c(newList, strsplit(i, '.', fixed = TRUE)[[1]][1])
}

gene_type_file <- 'gene_info.csv'
gene_type <- read.table(gene_type_file, header = TRUE, row.names = 1, sep=',')

ID <- bitr(newList, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene = as.vector(ID$ENTREZID), OrgDb=org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1, readable = TRUE)
head(summary(ego))

plot <- dotplot(ego, showCategory=15) + 
  theme(axis.text.y = element_text(size=8)) + 
  labs(title = "GO Enrichment Analysis of DSEAs")

deg_file <- "DEG_genelist_P-adj.gene_info.csv"
deg_data <- read.table(deg_file, header = TRUE, row.names = 1, sep=',')

significant_genes <- subset(deg_data, padj <= 0.05)
significant_genes <- subset(significant_genes, abs(log2FoldChange) >.75)

#Drop HOXD9
significant_genes <- subset(significant_genes, rownames(significant_genes) != "ENSG00000128709.11")

deg_gene_list <- significant_genes$gene_id

new_deg_gene_list <- vector()
for (i in deg_gene_list) {
  new_deg_gene_list <- c(new_deg_gene_list, strsplit(i, '.', fixed = TRUE)[[1]][1])
}

GO_results <- enrichGO(gene = new_deg_gene_list, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")

# ID <- bitr(new_deg_gene_list, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

selected_rows <- deg_data[deg_data$gene_id %in% genes_to_test, ]
print (selected_rows)

combine_list <- c(newList, new_deg_gene_list)
combine_ID <- bitr(combine_list, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

GO_results1 <- enrichGO(gene = as.vector(combine_ID$ENTREZID), OrgDb=org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", minGSSize = 1, pvalueCutoff = 0.05, qvalueCutoff = 1, readable = TRUE)

GO_results1 <- enrichGO(gene = combine_list, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")


