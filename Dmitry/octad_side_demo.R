library('octad')
library('ggrepel')

df = read.csv('/Users/leshchi4/colo_data/colo_expression.csv', row.names = 1)
df = log2(df + 1)
df = as.matrix(df)

coldata <- read.csv('/Users/leshchi4/colo_data/colo_coldata.csv', row.names = 1)

all(rownames(coldata) %in% colnames(df))

all(rownames(coldata) == colnames(df))

df <- df[, rownames(coldata)]
all(rownames(coldata) == colnames(df))

ctrl <- coldata[coldata$type == 'Solid Tissue Normal', ]
control_id <- rownames(ctrl)

case <- coldata[coldata$type == 'Primary Tumor', ]
case_id <- rownames(case)

res=diffExp(case_id,control_id,source='side',output=TRUE,n_topGenes=10000,expSet=df,annotate=FALSE,DE_method="limma")

res$gene_symbol <- rownames(res)
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 1 & res$pvalue < 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange < -1 & res$pvalue < 0.05] <- "DOWN"

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
res$delabel <- NA
res$delabel[res$diffexpressed != "NO"] <- res$gene_symbol[res$diffexpressed != "NO"]

ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

res$Symbol <- res$identifier

sRGES = runsRGES(res,max_gene_size=100,permutations=10000)

sRGESf = sRGES[sRGES$sRGES < -0.2,]
sRGESf = sRGESf[sRGESf$n > 3,]
sRGESf = sRGESf[!grepl("BRD-", sRGESf$pert_iname),]

folderEnrich = '/Users/leshchi4/colo_enrich'
octadDrugEnrichment(sRGES=sRGES, target = c('chembl_targets','mesh','ChemCluster'), outputFolder = folderEnrich)

path = '/Users/leshchi4/colo_data/res_colo.csv'
write.csv(res, path)

# modify this code to compare male and female samples
