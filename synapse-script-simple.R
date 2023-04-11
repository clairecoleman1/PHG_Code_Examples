#Contents:

# RNA-Seq: QC plots, sex checks and tsne
# ATAC-Seq: QC plots, sex checks and tsne

################################################################################

# RNA-Seq

# QC plots 
# Please download RNA-seq metadata from this repository. 
rnaqc <- read.csv("./Downloads/metadata-phg-RNA.csv")
celllabels <- c(astroAndMicroglia = "Astrocytes and Microglia", neuron = "Neurons", oligodendrocytes = "Oligodendrocytes")
library(ggplot2)
library(ggh4x)
require(gridExtra)
class(rnaqc)
rplot1 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=RIN, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("RIN") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) + theme(strip.background = element_blank(), strip.text = element_blank())
rplot1
rplot2 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$rnaseqc2_PCT_INTERGENIC_BASES, fill = cell_subtype)) + geom_boxplot() + theme_classic()  + guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Intergenic Rate") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
rplot2
rplot3 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$rnaseqc2_PCT_INTRONIC_BASES, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Intronic Rate") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
rplot3
rplot4 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$picard_MEDIAN_INSERT_SIZE, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Median Insert Size") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
rplot4
rplot5 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$picard_READ_PAIRS_EXAMINED, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Mapped Read Pairs") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
rplot5
rplot6 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$star_Uniquely_mapped_reads_pct, fill = cell_subtype)) + geom_boxplot() + theme_classic()   + guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("% Uniquely Mapped Reads") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
rplot6
rplot7 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$picard_meanGcContent, fill = cell_subtype)) + geom_boxplot() + theme_classic()+ guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Mean GC Content") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
rplot7
rplot8 <- ggplot(rnaqc, aes(x=Dx_asFactor, y=rnaqc$rnaseqAlign_PCT_RIBOSOMAL_BASES, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("% Ribosomal Bases") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086"))  +theme(strip.background = element_blank(), strip.text = element_blank())
rplot8
rgc <- grid.arrange(rplot1, rplot2, rplot3, rplot4, rplot5, rplot6, rplot7, rplot8, ncol=3)
rgc
#ggsave("dec1-rna-qc-8plot-pink-green-purple.pdf", agt)
#ggsave("dec1-rna-qc-8plot-pink-green-purple.png", agt)

################################################################################

# RNA-Seq

#Sex check 
ggplot(rnaqc, aes(x = XIST_expression, y = RPS4Y1_expression, color = Sex)) + geom_point() + theme_classic()
#ggsave("rna-sex-check.pdf", last_plot())


################################################################################

# RNA-Seq

#tsne

# Please download RNA-seq count matrix at https://www.synapse.org/#!Synapse:syn51138026 
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

plottsne <- DimPlot(tsne,reduction = "tsne", pt.size = 2) + ggtitle("tSNE")
plottsne

#ggsave("rna-tsne.pdf", last_plot())


################################################################################
# ATAC-Seq

# QC plots 
# Please download ATAC-seq metadata from this repository.
atacqc <- read.csv("./Downloads/metadata-phg-ATAC.csv")
celllabels <- c(microgliaAndAstrocytes = "Astrocytes and Microglia", neuron = "Neurons", oligodendrocytes = "Oligodendrocytes")
library(ggh4x)
require(gridExtra)
library(ggh4x)
class(atacqc)
aplot1 <- ggplot(atacqc, aes(x=Dx_asFactor, y=peakNarrowFDR1pctCount, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Number of Narrow Peaks") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) + theme(strip.background = element_blank(), strip.text = element_blank())
aplot1
aplot2 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$fracReadsInNonBlacklistedPeaks, fill = cell_subtype)) + geom_boxplot() + theme_classic()  + guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Fraction of reads in peaks") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot2
aplot3 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$chrMFrac, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("% mitDNA Reads") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot3
aplot4 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$insertMetrics_MEAN_INSERT_SIZE, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Median Insert Size") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot4
aplot5 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$star_Uniquely_mapped_reads_pct, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Fraction of uniquely mapped reads") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot5
aplot6 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$picard_meanGcContent, fill = cell_subtype)) + geom_boxplot() + theme_classic()+ guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Mean GC Content") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot6
aplot7 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$finalReadCount, fill = cell_subtype)) + geom_boxplot() + theme_classic() + guides(fill="none")+ facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Number of uniquely mapped reads") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot7
aplot8 <- ggplot(atacqc, aes(x=Dx_asFactor, y=atacqc$picard_PERCENT_DUPLICATION, fill = cell_subtype)) + geom_boxplot() + theme_classic()   + guides(fill="none") + facet_wrap2(vars(cell_subtype), labeller = as_labeller(celllabels)) + xlab("") + ylab("Fraction of duplicated reads") + labs(fill = "Diagnosis") + scale_fill_manual(values=c("#8434fb", "#0d8405","#fb3086")) +theme(strip.background = element_blank(), strip.text = element_blank())
aplot8


agt <- grid.arrange(aplot1, aplot2, aplot3, aplot4, aplot5, aplot6, aplot7, aplot8, ncol=3)
#ggsave("atac-qc-8plot-pink-green-purple.pdf", agt)
#ggsave("atac-qc-8plot-pink-green-purple.png", agt)


################################################################################

# ATAC-Seq

#Sex check 
ggplot(atacqc, aes(x = readCountInfo_Total, y = chryCounts, color = sex)) + geom_point() + theme_classic()
#ggsave("atac-sex-check.pdf", last_plot())

################################################################################
# Please download ATAC-seq count matrix at https://www.synapse.org/#!Synapse:syn51137300
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

################################################################################
