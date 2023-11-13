# This is the script to run DESeq2 on the count data
# The input is the count data from gene_count.csv

# Load the DESeq2 library
# library("DESeq2")

# Read the count data
# Install R and DESeq2. Upon installing R, install DESeq2 on R:
# source("https://bioconductor.org/biocLite.R")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

BiocManager::install("DESeq2")

biocLite("DESeq2")
# Import DESeq2 library in R
library("DESeq2")
# Load gene(/transcript) count matrix and labels
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <- read.csv(PHENO_DATA, sep="\t", row.names=1)
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
# Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
# [1] TRUE
countData <- countData[, rownames(colData)]
# all(rownames(colData) == colnames(countData))
print(all(rownames(colData) == colnames(countData)))
# [1] TRUE
# Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData,
colData = colData, design = ~ CHOOSE_FEATURE)
# Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
res <- results(dds)
# Sort by adjusted p-value and display
print(res[order(res$padj), ])
# (resOrdered <- res[order(res$padj), ])