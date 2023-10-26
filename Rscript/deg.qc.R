library(data.table)
library(ggplot2)
library(tximport)
library(jsonlite)
library(edgeR)
library(ggrepel)

setwd("C:/Irfan/Personal/Masters/BS6202/Project")

############################################
## RNA-seq (QC) Project by Muhammad Irfan ##
############################################

count_file <- "BS6202.count.mx"
data <- read.table(count_file, header = TRUE, row.names = 1)
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
cts <- data+1

group <- factor(c("vector", "vector", "vector", "hoxd9", "hoxd9", "hoxd9"))
y <- DGEList(counts = cts, group = group)
y <- calcNormFactors(y)

library(reshape2)
df_long <- melt(as.data.frame(y$counts))

p <- ggplot(df_long, aes(x = variable, y = value)) + geom_boxplot() + theme_bw() + ggtitle("Box plot of read counts") + xlab("sample") + ylab("read counts") + scale_y_log10()

y <- t(y$counts)
constant_columns <- apply(y, 2, var) == 0
y <- y[, !constant_columns]

fit <- prcomp(y, scale. = TRUE, center = TRUE)
pca_pcs <- cbind(group, as.data.frame(fit$x))
pca_pcs$sample <- rownames(pca_pcs)
pca_var <- data.frame(pc=1:length(fit$sdev), variance=(fit$sdev^2)/sum(fit$sdev^2))
pca_loadings <- cbind(parameter=rownames(fit$rotation), as.data.frame(fit$rotation))

p <- ggplot(pca_var, aes(x=pc, y=variance*100)) + geom_col() + theme_bw() + ggtitle("Variance accounted for")
p <- ggplot(pca_pcs, aes(x=PC1, y=PC2, col=group, label=sample)) + geom_point() + geom_text_repel(col="black") + theme_bw() + ggtitle("PCA plot") + xlab(sprintf("PC1 (%0.3f%%)", pca_var[1, 2]*100)) + ylab(sprintf("PC2 (%0.3f%%)", pca_var[2, 2]*100))

##################
## Test in RPKM ##
##################

e <- DGEList(cts)

# Calculate library sizes
e <- calcNormFactors(e)

length_file <- "BS6202.GeneLength.tsv"
gene_lengths <- as.data.frame(read.table(length_file, header = TRUE, row.names = 1))

common_rows <- intersect(rownames(cts), rownames(gene_lengths))
common_genes <- gene_lengths[common_rows, , drop = FALSE]
e$genes <- common_genes

e_rpkm <- rpkm(e, log=T)
e_rpkm_long <- melt(as.data.frame(e_rpkm))

p <- ggplot(e_rpkm_long, aes(x = variable, y = value)) + geom_boxplot() + theme_bw() + ggtitle("Box plot of read counts") + xlab("sample") + ylab("read counts") + scale_y_log10()

y <- t(e_rpkm)
constant_columns <- apply(y, 2, var) == 0
y <- y[, !constant_columns]

fit <- prcomp(y, scale. = TRUE, center = TRUE)
pca_pcs <- cbind(group, as.data.frame(fit$x))
pca_pcs$sample <- rownames(pca_pcs)
pca_var <- data.frame(pc=1:length(fit$sdev), variance=(fit$sdev^2)/sum(fit$sdev^2))
pca_loadings <- cbind(parameter=rownames(fit$rotation), as.data.frame(fit$rotation))

p <- ggplot(pca_var, aes(x=pc, y=variance*100)) + geom_col() + theme_bw() + ggtitle("Variance accounted for")
p <- ggplot(pca_pcs, aes(x=PC1, y=PC2, col=group, label=sample)) + geom_point() + geom_text_repel(col="black") + theme_bw() + ggtitle("PCA plot") + xlab(sprintf("PC1 (%0.3f%%)", pca_var[1, 2]*100)) + ylab(sprintf("PC2 (%0.3f%%)", pca_var[2, 2]*100))
