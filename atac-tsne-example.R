#Download atacseq_count_matrix_norm.tsv in advance
readCount <- read.csv("./Downloads/atacseq_count_matrix_norm.tsv", sep = "\t")
head(readCount)
#install.packages("ggplot2")
library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)

#Seurat prep
phga <- CreateSeuratObject(counts = readCount, project = "phga", assay = "ATAC")
phga
phga <- NormalizeData(phga)
phga <- FindVariableFeatures(phga, selection.method = "vst")
all.genes <- rownames(phga)
phga <- ScaleData(phga, features = all.genes)
phga <- RunPCA(phga, features = VariableFeatures(object = phga))
print(phga[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(phga, dims = 1:2, reduction = "pca")
DimPlot(phga, reduction = "pca")
phga <- FindNeighbors(phga, dims = 1:10)
phga <- FindClusters(phga, resolution = 0.5)
phga <- RunUMAP(phga, dims = 1:10)

#Plot
plot <- RunTSNE(
  phga,
  cells = NULL,
  dims = 1:10,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = NULL,
  reduction.name = "tsne",
  reduction.key = "tSNE_",
)
plot
DimPlot(plot,reduction = "tsne", pt.size = 2) + ggtitle("tSNE")

