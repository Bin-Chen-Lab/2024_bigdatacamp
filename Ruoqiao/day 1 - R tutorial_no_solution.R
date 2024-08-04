#Author: Ruoqiao Chen

#HPCC OnDemand platform: https://ondemand.hpcc.msu.edu/

#R Script: where you write your codes. A R script is saved as a ".R" file.

#Shortcut key of running a R code line: Use "command + return" on your Mac keyboard. 

#Today we'll use the data "OCTAD_cell_line.RData" to introduce basic R coding and how to use R to analyze your data.

#---------------------------------------------------------------------------------------------------------------------
#read and save data:
#read a RData file:
load("/mnt/research/IPSTP_2024/Unix_and_R/OCTAD_cell_line.RData")

#Data viewer:
View(octad_cell_line_matrix)
View(octad_cell_line_features)
View(octad_cell_line_meta)

#---------------------------------------------------------------------------------------------------------------------
#About the "OCTAD_cell_line.RData" data:
#"octad_cell_line_matrix": Rows are cell lines, columns consist of multiple types of measurements including gene mutation, CNV, gene expression, protein level, metabolite level, gene effect score, and drug sensitivity.
#"octad_cell_line_features": Information on each column (i.e., measurement name) of "octad_cell_line_matrix". 
#"octad_cell_line_meta": Information on each row (i.e., cell line) of "octad_cell_line_matrix". 
#---------------------------------------------------------------------------------------------------------------------

#save an object as a csv file:
write.csv(octad_cell_line_meta, "octad_cell_line_meta.csv") 
#save an object as a txt file:
write.table(octad_cell_line_meta, "octad_cell_line_meta.txt") 
#read a csv file:
octad_cell_line_meta_new = read.csv("octad_cell_line_meta.csv", stringsAsFactors = FALSE, check.names=FALSE) 
#read a txt file:
octad_cell_line_meta_new = read.table("octad_cell_line_meta.txt", stringsAsFactors = FALSE) 
#save objects as a RData file:
save(octad_cell_line_meta, octad_cell_line_meta_new, file = "octad_cell_line_meta.RData")

#---------------------------------------------------------------------------------------------------------------------
#check the user help: 
?read.table

#---------------------------------------------------------------------------------------------------------------------
#R packages:

#load a package:
library(dplyr)
library(ggplot2)

#If the package is not there, you need to install it first before you load it:
#install R package from CRAN:
install.packages("dplyr")
install.packages("ggplot2")
#install bioconductor R packages:
install.packages("BiocManager")
BiocManager::install("GEOquery")

#---------------------------------------------------------------------------------------------------------------------
#data type:
#tips: To name an object: 1) no space, 2) no special character such as #, &, 3) don't start with number, and 4) understandable 

#numeric:
age = 60
#show the class of an object:
class(age)
#rounding of numbers:
round(0.267, digits = 2)

#Basic Operators and Calculations:
age_new = age + 20
age_new = age - 20
age_new = age / 3
age_new = age * 3

#character:
cancer = "liver cancer"
class(cancer)
age_2 = "60"
class(age_2)
#concatenate two objects after converting to characters:
paste(cancer, age_2, sep = '')
paste(cancer, age, sep = '')
paste(cancer, age, sep = ',')
tmp = strsplit(paste(cancer, age, sep = ','), ',')
unlist(tmp)[1]
unlist(tmp)[2]

#We can't do calculations on character objects! 
age_new = age_2 + 20
#convert one class to another class (e.g., character to numeric before you do calculations):
age_new = as.numeric(age_2) + 20
as.character(age_new)

#logical (TRUE/FALSE):
age == 60
age != 60
age < 80
age > 80
#
#"if" statement:
if (age == 40){
  print("age == 40")
}else{
  print("age != 40")
}

#vector:
cancers = c("liver cancer", "lung cancer")
ages = c(60, 30:34)
ages[1]
ages[c(1,3)]
ages[-c(1,2)]
sort(ages, decreasing = T)
rev(ages)
#Generate a numeric sequence:
seq(1, 10, by = 2)
#repeat an element N times:
rep(2, 5)
rep(cancers, 5)
#show unique elements:
unique(rep(cancers, 5))
unique(rep(2, 5))
#
#Data stats summary:
ages
mean(ages)
median(ages)
sum(ages)
log(ages)
min(ages)
max(ages)
#
#"for" loop:
for (x in ages){
  print(paste("this is", x))
}
#
for (x in cancers){
  print(paste("this is", x))
}

#matrix / data frame:
class(octad_cell_line_matrix)
#
#show dimensions of a matrix:
dim(octad_cell_line_matrix)
#
#show the number of rows / columns:
nrow(octad_cell_line_matrix)
ncol(octad_cell_line_matrix)
#
#show column names / row names:
colnames(octad_cell_line_matrix)
rownames(octad_cell_line_matrix)
#
#subset a matrix:
octad_cell_line_matrix[c(2:3, 7), 1:2]
octad_cell_line_matrix[c('ACH-000002', 'ACH-000003', 'ACH-000007'), c('mutation_VPS13D', 'mutation_AADACL4')]
octad_cell_line_matrix[, 1:2]
octad_cell_line_matrix[ , c('mutation_VPS13D', 'mutation_AADACL4')]
octad_cell_line_matrix[c(2:3, 7), ]
octad_cell_line_matrix[c('ACH-000002', 'ACH-000003', 'ACH-000007'), ]
#
#Data frame:
octad_cell_line_matrix = as.data.frame(octad_cell_line_matrix)
class(octad_cell_line_matrix)
octad_cell_line_matrix[c(2:3, 7), 1:2]
octad_cell_line_matrix[c('ACH-000002', 'ACH-000003', 'ACH-000007'), c('mutation_VPS13D', 'mutation_AADACL4')]
octad_cell_line_matrix$mutation_VPS13D
#
#rename rows of a data frame:
rownames(octad_cell_line_meta) = octad_cell_line_meta$DepMap_ID
#
#transpose a data frame / matrix:
t(octad_cell_line_meta)
View(t(octad_cell_line_meta))
#
#combine two data frames by columns / rows:
new_data = cbind(octad_cell_line_meta, octad_cell_line_matrix[, 1:5])
dim(new_data)
new_data2 = rbind(t(octad_cell_line_meta), t(octad_cell_line_matrix[, 1:5]))
dim(new_data2)
#
#use the "dplyr" package to assist with subsetting a data frame:
library(dplyr)
#e.g., obtain a subset containing only breast cancer, pancreatic cancer, and liver cancer cell lines:
df = dplyr::filter(octad_cell_line_meta, disease == 'Breast Cancer' | disease == 'Pancreatic Cancer' | disease == 'Liver Cancer')
unique(df$disease)
unique(df$lineage)
#e.g., remove rows containing undefined lineages:
df = df[df$lineage != '', ]
unique(df$lineage)

#Factor:
factor(df$disease) # 3 disease factor levels
factor(df$lineage) # 3 cell line lineage factor levels

#summarize the counts for each factor (or each combination of factors):
table(df$disease)
table(df$disease, df$lineage)

#---------------------------------------------------------------------------------------------------------------------
#statistical analysis:

#Fisher's exact test:
fisher.test(factor(df$disease), factor(df$lineage))

#t-tests & Wilcoxon tests:
#compare MYC expression in liver cancer and breast cancer
MYC_liver = octad_cell_line_matrix[rownames(octad_cell_line_meta[octad_cell_line_meta$disease == "Liver Cancer", ]), "expression_MYC"]
MYC_liver
#remove NA values from a vector:
MYC_liver = MYC_liver[!is.na(MYC_liver)]
MYC_breast = octad_cell_line_matrix[rownames(octad_cell_line_meta[octad_cell_line_meta$disease == "Breast Cancer", ]), "expression_MYC"]
MYC_breast
MYC_breast = MYC_breast[!is.na(MYC_breast)]
boxplot(MYC_liver, MYC_breast)
#
#t-test (unpaired, two-sided by default):
t.test(MYC_liver, MYC_breast)
?t.test
#
#Wilcoxon test (unpaired, two-sided by default):
wilcox.test(MYC_liver, MYC_breast)
?wilcox.test

#correlation analysis:
#compute correlation between MYC and MAX (a.k.a., myc-associated factor X) gene expression:
#Pearson correlation:
cor.test(octad_cell_line_matrix[, "expression_MYC"], octad_cell_line_matrix[, "expression_MAX"], method = "pearson")
#Spearman correlation:
cor.test(octad_cell_line_matrix[, "expression_MYC"], octad_cell_line_matrix[, "expression_MAX"], method = "spearman")

#---------------------------------------------------------------------------------------------------------------------
#Advanced: Use functions to simplify your code writing when encountering repetitive tasks:
#e.g., concatenate any possible input character strings with another character string "cancer":
#write the function:
my_function <- function(fname) {
  paste(fname, "cancer")
}
#
#Execute the function on any possible new inputs:
my_function("Liver")
my_function("Lung")
my_function("Pancreas")
#---------------------------------------------------------------------------------------------------------------------
#visualization

#A simple scatter plot:
plot(octad_cell_line_matrix[, "expression_MYC"], octad_cell_line_matrix[, "expression_MAX"])

#A simple boxplot:
boxplot(MYC_liver, MYC_breast)

#Advanced: Use "ggplot2" package to refine your plots:
#e.g., Draw a scatter plot for any two interested genes' expression:
plot_gene_coexpression <- function(gene1, gene2){
  library(ggplot2)
  plot = ggplot(octad_cell_line_matrix, aes(octad_cell_line_matrix[, gene1], octad_cell_line_matrix[, gene2])) + 
    geom_point(shape=1, color = 'blue') +
    xlab(gene1) +
    ylab(gene2) +
    ggtitle(paste0('Cor: ', round(cor.test(octad_cell_line_matrix[, gene1], octad_cell_line_matrix[, gene2], method = 'pearson')$estimate, 2),
                   ' (p=', round(cor.test(octad_cell_line_matrix[, gene1], octad_cell_line_matrix[, gene2], method = 'pearson')$p.value, 2), ')')) +
    theme_bw()
  print(plot)
}
#
plot_gene_coexpression('expression_MYC', 'expression_MAX')
plot_gene_coexpression('expression_PTEN', 'expression_AKT1')
plot_gene_coexpression('expression_MTOR', 'expression_AKT1')

#PCA:
#e.g., here we want to do PCA on the gene expression data in breast cancer, pancreatic cancer, and liver cancer cell lines, and then visualize the PCA result:
df = dplyr::filter(octad_cell_line_meta, disease == 'Breast Cancer' | disease == 'Pancreatic Cancer' | disease == 'Liver Cancer')
cell_lines = rownames(df)
exp_features = as.character(octad_cell_line_features$id[octad_cell_line_features$type == "expression"])
data = octad_cell_line_matrix[cell_lines, exp_features]
#
#Apply a function to each row of a matrix / data frame:
#e.g., here we want to summarize the counts of NA values for each row:
na_samples = apply(data, 1, function(x){sum(is.na(x))})
#only use rows with zero NA values:
data = data[na_samples == 0, ]
dim(data)
#
#Run PCA. The data's feature dimension reduces after running PCA:
pca = prcomp(data)
dim(pca$x)
#
#Visualize the PCA results:
library(ggplot2)
tmp = cbind(pca$x, df[rownames(data), ])
plot_pca = ggplot(tmp, aes(PC1, PC2, color = disease)) + 
      geom_point(shape=1, size = 3) +
      theme_bw()
print(plot_pca)

#Q: What does the PCA plot indicate?

#---------------------------------------------------------------------------------------------------------------------
#Practice:

#Q1: Is there co-expression between the ERBB2 (a.k.a., HER2) gene and the AKT1 gene in breast cancer cell lines? How about in pancreatic cancer cell lines?

#Q2: Perform PCA on the drug sensitivity data in lung cancer, colon/colorectal Cancer, and pancreatic cancer cell lines. Visualize the PCA results.

