atacqc <- read.csv("./Downloads/metadata-phg-ATAC.csv")
celllabels <- c(microgliaAndAstrocytes = "Astrocytes and Microglia", neuron = "Neurons", oligodendrocytes = "Oligodendrocytes")
library(ggh4x)
require(gridExtra)
library(ggh4x)

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
