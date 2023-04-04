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