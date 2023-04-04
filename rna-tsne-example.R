#Download rnaseq_count_matrix_norm.tsv in advance
readCount <- read.csv("./Downloads/rnaseq_count_matrix_norm.tsv", sep = "\t")
#install.packages("ggplot2")
library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)

#Seurat prep
phgr <- CreateSeuratObject(counts = readCount, project = "phgr")
phgr
phgr <- NormalizeData(phgr)
phgr <- FindVariableFeatures(phgr, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(phgr), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(phgr)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(phgr)
phgr <- ScaleData(phgr, features = all.genes)
phgr <- RunPCA(phgr, features = VariableFeatures(object = phgr))
print(phgr[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(phgr, dims = 1:2, reduction = "pca")
DimPlot(phgr, reduction = "pca")
phgr <- FindNeighbors(phgr, dims = 1:10)
phgr <- FindClusters(phgr, resolution = 0.5)
phgr <- RunUMAP(phgr, dims = 1:10)
# S3 method for Seurat
tsne <- RunTSNE(
  phgr,
  cells = NULL,
  dims = 1:30,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
)
tsne
DimPlot(tsne,reduction = "tsne", pt.size = 2) + ggtitle("tSNE")
