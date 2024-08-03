count_data <- read.csv(gzfile("~/Downloads/GSE187420_Raw_Read_count.csv.gz"), row.names = 1)
tpm_data <- read.csv(gzfile("~/Downloads/GSE187420_Normalized_TPM_count.csv.gz"), row.names = 1)

######convert gene ids
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = rownames(count_data),
                 mart = ensembl)
count_data = merge( mapping, count_data, by.x = "ensembl_gene_id", by.y = 0)
count_data = count_data[, -1]
count_data = aggregate(. ~ hgnc_symbol, count_data, median)
count_data = count_data[count_data$hgnc_symbol != "", ]
rownames(count_data) = count_data$hgnc_symbol
count_data = count_data[, -1]

tpm_data = merge( mapping, tpm_data, by.x = "ensembl_gene_id", by.y = 0)
tpm_data = tpm_data[, -1]
tpm_data = aggregate(. ~ hgnc_symbol, tpm_data, median)
tpm_data = tpm_data[tpm_data$hgnc_symbol != "", ]
rownames(tpm_data) = tpm_data$hgnc_symbol
tpm_data = tpm_data[, -1]

################
#PCA
library(ggplot2)
library(readr)
library(FactoMineR)
library(factoextra)
pca_results <- PCA(t(tpm_data), scale.unit = TRUE, graph = FALSE)

# Extract variance explained by each principal component
variance_explained <- pca_results$eig

# Extract PCA coordinates
pca_coords <- pca_results$ind$coord

# Create a data frame for PCA coordinates
pca_df <- data.frame(pca_coords)
pca_df$sample <- rownames(pca_df)

# Plot PCA
ggplot(pca_df, aes(x = Dim.1, y = Dim.2, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1) +
  xlab(paste0("PC1: ", round(variance_explained[1, 2], 2), "% variance")) +
  ylab(paste0("PC2: ", round(variance_explained[2, 2], 2), "% variance")) +
  theme_minimal() +
  ggtitle("PCA of TPM Data")

####################
#DE analysis
library(DESeq2)

sample_info = data.frame(condition = c("Control", "Control", "Control", "Infection", "Infection", "Infection", "Treatment", "Treatment", "Treatment"))
sample_info$condition = factor(sample_info$condition, levels = c("Control", "Infection", "Treatment"))

rownames(sample_info) = colnames(count_data)
count_data_raw = round(2^(count_data+0.1))

hist(count_data_raw$Control.1)
hist(count_data_raw$Control.1)

dds <- DESeqDataSetFromMatrix(countData = count_data_raw,
                              colData = sample_info,
                              design = ~ condition)  
# Run the DESeq pipeline
dds <- DESeq(dds)

# Get the results
res <- results(dds, contrast=c("condition","Infection","Control"))

plotMA(res, ylim=c(-2,2))

# View a summary of the results
summary(res, alpha = 0.05)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj),]

# View the top differentially expressed genes
head(res_ordered)

# Convert results to a data frame for ggplot2
res_df <- as.data.frame(res_ordered)

# Create a volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value") +
  ggtitle("Volcano Plot of Differential Expression Analysis") +
  theme_minimal()

res_df = res_df[order(res_df$log2FoldChange), ]

# Perform GSEA using KEGG pathways
library("clusterProfiler")
library(org.Hs.eg.db)

# Example gene expression data (replace with your actual data)
gene_list <- data.frame(gene = rownames(res_df), logFC = res_df$log2FoldChange)
gene_list <- gene_list[order(gene_list$logFC, decreasing = TRUE),]  # Sort by logFC
gene_list <- setNames(gene_list$logFC, gene_list$gene)  # Named vector

# Run GSEA using GO terms
gsea_go_results <- gseGO(geneList = gene_list,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "ALL",  # You can also specify "BP", "MF", or "CC"
                         nPerm = 1000,
                         minGSSize = 10,
                         pvalueCutoff = 0.1,
                         verbose = FALSE)

# Print summary of results
summary(gsea_go_results)

# Plot enrichment for the top pathway
gseaplot(gsea_go_results, geneSetID = gsea_go_results@result$ID[1])

# Dot plot for top pathways
dotplot(gsea_go_results)

#################
