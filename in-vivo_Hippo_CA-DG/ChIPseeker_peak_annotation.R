### peak annotation with chipseeker ###

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(diffloop)
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

setwd("/data/manke/group/ferrari/PhD_project/reproducible_code_paper/Foxg1_code/in-vivo_Hippo_CA-DG/")

ca_peaks = readPeakFile("INTERMEDIATE_FILES/CA1-1_FOXG1_mm_i69.filtered.BAM_peaks.filtered.bed")
dg_peaks = readPeakFile("INTERMEDIATE_FILES/DG_FOXG1_mm_i68.filtered.BAM_peaks.filtered.bed")
consensus_peaks = readPeakFile("INTERMEDIATE_FILES/consensus_DG-CA_filtered.bed")
intersect_peaks = readPeakFile("INTERMEDIATE_FILES/intersect_CA-DG.filtered.bed")

ca_peaks = addchr(ca_peaks)
dg_peaks = addchr(dg_peaks)
consensus_peaks = addchr(consensus_peaks)
intersect_peaks = addchr(intersect_peaks)

peakAnno_ca <- annotatePeak(ca_peaks, tssRegion=c(-1000, 0),
                        TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_dg <- annotatePeak(dg_peaks, tssRegion=c(-1000, 0),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_consensus <- annotatePeak(consensus_peaks, tssRegion=c(-1000, 0),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_intersect <- annotatePeak(intersect_peaks, tssRegion=c(-1000, 0),
                                   TxDb=txdb, annoDb="org.Mm.eg.db")

# pdf("FIGURES/Annotation_CA.pdf")
# plotAnnoPie(peakAnno_ca) 
# dev.off()
# 
# pdf("FIGURES/Annotation_DG.pdf")
# plotAnnoPie(peakAnno_dg)
# dev.off()

# pdf("FIGURES/Annotation_consensus.pdf")
# plotAnnoPie(peakAnno_consensus)
# dev.off()
# 
# pdf("FIGURES/Annotation_intersect.pdf")
# plotAnnoPie(peakAnno_intersect)
# dev.off()

pdf("FIGURES/upsetplot_CA.pdf",width=10,height = 5, onefile=FALSE)
upsetplot(peakAnno_ca,vennpie=T)
dev.off()

pdf("FIGURES/upsetplot_DG.pdf", width=10,height = 5, onefile=FALSE)
upsetplot(peakAnno_dg,vennpie=T)
dev.off()

annotation_ca = as.data.frame(peakAnno_ca)
annotation_dg = as.data.frame(peakAnno_dg)
annotation_consensus = as.data.frame(peakAnno_consensus)
annotation_intersect = as.data.frame(peakAnno_intersect)

write.table(annotation_ca,"INTERMEDIATE_FILES/ChIPseeker_Annotation_CApeaks.tsv",sep="\t",quote=F,row.names = F)
write.table(annotation_dg,"INTERMEDIATE_FILES/ChIPseeker_Annotation_DGpeaks.tsv",sep="\t",quote=F,row.names = F)
write.table(annotation_consensus,"INTERMEDIATE_FILES/ChIPseeker_Annotation_consensus.tsv",sep="\t",quote=F,row.names = F)
write.table(annotation_intersect,"INTERMEDIATE_FILES/ChIPseeker_Annotation_intersect.tsv",sep="\t",quote=F,row.names = F)

