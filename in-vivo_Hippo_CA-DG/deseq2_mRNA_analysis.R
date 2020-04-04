############################
### Deseq2 analysis mRNA ###
############################
libr = "/home/ferrari/R/x86_64-redhat-linux-gnu-library/3.5/"
suppressMessages(library(argparse, lib=libr))

# create parser object
parser <- ArgumentParser()

parser$add_argument("-c", "--count_matrix", type="character", 
                    help="count matrix")
parser$add_argument("-meta", "--meta_data", type="character", 
                    help="matadata table")
parser$add_argument('-d','--design', type="character", 
                    help="design for comparision")
parser$add_argument('-ref','--reference_condition', type="character", default=NULL,
                    help="reference condition in a comparison")
parser$add_argument('-lfcThr','--lfc_threshold', type="character", default=NULL,
                    help="lfc threshold for test")
parser$add_argument('-o','--outdir', type="character", 
                    help="output directory")
parser$add_argument("-go_pca", "--analyze_PCA", action="store_true", default=FALSE,
                    help="perform GO enrichment analysis of loadings correlated with first 3 principal components")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

### LOAD PACKAGES ###

suppressMessages(library(stringr))
suppressMessages(library(DESeq2, lib=libr))
suppressMessages(library(vsn, lib=libr))
suppressMessages(library(ggplot2, lib=libr))
library(FactoMineR, lib=libr)
library(org.Mm.eg.db, lib="/data/manke/group/ferrari/3.5/")
library(clusterProfiler, lib="/data/manke/group/ferrari/3.5/")
# library(AnnotationDbi, lib=libr)
#library(MSigDB)
# library(gskb, lib=libr)
# library(gage, lib=libr)
# library(ggpubr, lib=libr)
suppressMessages(library(RColorBrewer, lib=libr))
suppressMessages(library(pheatmap, lib=libr))
suppressMessages(library(plyr, lib=libr))


### auxiliary functions start ###

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

load_housekeeping = function(){
  dk = read.csv("/data/manke/group/ferrari/my_repository/housekeeping_genes.txt",sep="\t", header=F, comment.char = "#")
  colnames(dk) = c("symbol","id")
  dk$symbol = trim(dk$symbol)
  dk$capSymbol = str_to_title(dk$symbol)
  return(dk$capSymbol)
}

### auxiliary functions end ###

#setwd("/data/manke/group/ferrari/PhD_project/reference_datasets/Ferrari_mESC_DMSOvsEPZ_RNA-Seq/downstream_analysis/DESeq2/")
out_dir = args$outdir
print(out_dir)
ifelse(!dir.exists(out_dir), dir.create(out_dir), FALSE)

file.copy(args$meta_data, paste(out_dir,"metadata.tsv",sep='/'))
file.copy(args$count_matrix, paste(out_dir,"counts.tsv",sep='/'))
cat(args$design,file=paste(out_dir,"design.txt",sep='/'))
          

# ### prepare design table ###
print(args$meta_data)
design_table = read.csv(args$meta_data, sep="\t", header = T, comment.char = '#')
#design_table = read.csv("metadata_TSA_shFoxg1VSCtr.tsv", sep="\t", header = T, comment.char = '#')
rownames(design_table) = design_table$sample
design_table$sample = NULL



# ### LOAD COUNT TABLE ###
print(args$count_matrix)
count=read.csv(args$count_matrix, sep='\t', header = T, comment.char = '#')
#count=read.csv("../output_snakePipes_hippoShFOXG1KD_TSA/featureCounts/counts.tsv", sep='\t', header = T, comment.char = '#')

rownames(count) = count$X
count$X = NULL
print(colnames(count))
print(rownames(design_table))
select_cols = as.vector(row.names(design_table))
count = count[,select_cols]
print(head(count))
# 
# 
# 
#                             #######################
#                             ### mRNA processing ###
#                             #######################
# 
# 

print(args$design)

# 
print(all(rownames(design_table) == colnames(count)))
# 
# 
# ### create dds object ###
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = design_table,
                              design = as.formula(args$design))
                              #design = ~ condition)
                              
dds
# 
# ### add feature data ### 
suppressMessages(library(rtracklayer))

anno18 <- import.gff2("/data/manke/group/ferrari/my_repository/annotations_gencode/mouse/M18/gencode.vM18.annotation.sorted.gtf", feature.type = "gene")

featureData = data.frame(row.names = anno18$gene_id,
                         gene_name = anno18$gene_name,
                         gene_type = anno18$gene_type,
                         level = anno18$level,
                         chr = seqnames(anno18))
print(head(row.names(dds)))
print(head(row.names(featureData)))
index=match(row.names(dds),row.names(featureData))
mcols(dds) = DataFrame(mcols(dds), featureData[index,])
mcols(dds)


### PREFILTERING ###
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

### SET REFERENCE LEVELS ###
if (is.null(args$reference_condition)){
  print(dds$condition)
} else {
  dds$condition = relevel(dds$condition, ref = args$reference_condition)
  print(dds$condition)
  #
}
 

# 
# 
### SAMPLE CLUSTERING ###
#########################

rld <- rlog(dds, blind=FALSE)

### Mean/sd plot ###
pdf(paste(out_dir,'Mean_SDplot.pdf',sep="/"), height = 5)
meanSdPlot(assay(rld))
dev.off()

### Compute eucledian distance between samples ###
sampleDists = dist(t(assay(rld)))
### hierarchical clustering of samples based on distance ###
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = colnames(dds)
colnames(sampleDistMatrix) = NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Reds")) )(255)

pdf(paste(out_dir,"eucledian_distance_clustering.pdf",sep="/"), height = 5 , onefile = F)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

### PCA ###
pdf(paste(out_dir,'PCA.pdf',sep="/"))
plotPCA(rld, intgroup=c('condition'))
dev.off()

A=plotPCA(rld, ntop=500, returnData=T, intgroup=c("condition"))

write.table(A, paste(out_dir,"PCA_coords.tsv",sep="/"), sep="\t", quote = F)

percentVar <- round(100 * attr(A, "percentVar"))
pca_plt=ggplot(A, aes(x=PC1, y=PC2, col=name))+
  geom_point(size=4)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

pdf(paste(out_dir,"PCA_detailed.pdf",sep="/"))
pca_plt
dev.off()
#
# # ### extract features that mostly correlate with PC1 and PC2
matrix_rlog=as.data.frame(t(assay(rld)))
variances=apply(matrix_rlog, 2, var)
variances=variances[order(variances, decreasing = T)]
matrix_rlog_1 = matrix_rlog[,names(matrix_rlog) %in% names(variances)[1:500]]



# # 
# # pca=prcomp(matrix_rlog)
# # pca_1=prcomp(matrix_rlog_1)
# # summary_pca=(summary(pca))
# # summary_pca_1=summary(pca_1)
# # plot(pca$rotation[,2]~pca$rotation[,1])
# # pdf('output_DESeq2_analysis/PCA_2vs3.pdf', height= 4, width = 8)
# # par(mfrow=c(1,2))
# # plot(pca_1$x[,2]~pca_1$x[,1], col=c(rep('blue',3),rep('red',3)), pch=20, xlab='PC1', ylab='PC2'); #legend("topleft", legend=c("EPZ", "DMSO"), col=c(2,4), pch=20)
# # plot(pca_1$x[,3]~pca_1$x[,2], col=c(rep('blue',3),rep('red',3)), pch=20, xlab='PC2', ylab='PC3'); legend("topright", legend=c("EPZ", "DMSO"), col=c(2,4), pch=20)
# # dev.off()
# # 
# # 
# # loadings <- data.frame(gene_ens=rownames(pca$rotation), pca$rotation[,1:2]^2, stringsAsFactors=F)
# # loadings <- loadings[order(loadings$PC1, decreasing=T),]
# # loadings_def = merge(loadings, featureData, by='row.names')
# # loadings_def_1 = loadings_def[order(loadings_def$PC1, decreasing = T),]
# # loadings_def_2 = loadings_def[order(loadings_def$PC2, decreasing = T),]
# # top_var_genes_1=as.character(head(unique(loadings_def_1$gene_name),n=300))
# # top_var_genes_2=as.character(head(unique(loadings_def_2$gene_name),n=300))
# # 
# # rownames(loadings) <- NULL
# # head(loadings)
# # png("screeplot.png", w=600, h=600)
# # screeplot(pca_1, type="lines", pch=16, main="Scree plot")
# # dev.off()
# # png("loadings.png", w=800, h=500)
# # par(mfrow=c(1,2))
# # plot(sort(loadings$PC1, decr=T), pch=20, main="PC1 loadings", ylab='loadings', xlab = 'gene index'); points(sort(loadings$PC1, decr=T)[1:50], col=2, pch=20); legend("topright", legend=c("Top 500 Genes", "Other Genes"), col=c(2,1), pch=20)
# # plot(sort(loadings$PC2, decr=T), pch=20, main="PC2 loadings", ylab='loadings', xlab = 'gene index'); points(sort(loadings$PC2, decr=T)[1:50], col=2, pch=20); legend("topright", legend=c("Top 500 Genes", "Other Genes"), col=c(2,1), pch=20)
# # dev.off()
# # par(mfrow=c(1,1))
# # 
# # biplot(pca)
# # summary(pca)
# # 
# # 
# ### PCA with FactoMineR ###

if (args$analyze_PCA){

pca_facto = FactoMineR::PCA(matrix_rlog_1)
result_pca_facto = as.list(summary(pca_facto,nbelements=Inf))
correlation_with_PC = FactoMineR::dimdesc(pca_facto, proba = 0.1)

sig_cor_PC1 = as.vector(featureData[rownames(as.data.frame(correlation_with_PC[[1]]$quanti)),]$gene_name)
sig_cor_PC2 = as.vector(featureData[rownames(as.data.frame(correlation_with_PC[[2]]$quanti)),]$gene_name)
sig_cor_PC3 = as.vector(featureData[rownames(as.data.frame(correlation_with_PC[[3]]$quanti)),]$gene_name)

eg_1 = bitr(sig_cor_PC1, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
universe = bitr(featureData$gene_name, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")

ego <- enrichGO(gene          = eg_1[["ENTREZID"]],
                universe      = universe[["ENTREZID"]],
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

barplot(ego, drop=TRUE, showCategory=12)

eg_2 = bitr(sig_cor_PC2, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")

ego_2 <- enrichGO(gene          = eg_2[["ENTREZID"]],
                universe      = universe[["ENTREZID"]],
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

barplot(ego_2)

eg_3 = bitr(sig_cor_PC3, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")

list_genes = list(PC1_corrGenes = eg_1$ENTREZID, PC2_corrGenes = eg_2$ENTREZID, PC3_corrGenes = eg_3$ENTREZID)

ck <- compareCluster(geneCluster = list_genes, fun = "enrichGO", OrgDb=org.Mm.eg.db, ont="BP",pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
head(as.data.frame(ck))
write.csv(as.data.frame(ck), paste(out_dir,"PCA_GOEnrichment.tsv",sep='/'), quote=F, sep="\t")

pdf(paste(out_dir,"Enrichment_PCcorrGenes.pdf",sep="/"))
dotplot(ck)
dev.off()
}

### DIFFERENTIAL EXPRESSION ANALYSIS ###

dds = DESeq(dds)

print(resultsNames(dds))

### without lfc threshold
res=results(dds)
resLFCShrink_normal = lfcShrink(dds, coef=tail(resultsNames(dds),n=1), type="normal")
resLFCShrink_apeglm = lfcShrink(dds, coef=tail(resultsNames(dds),n=1), type="apeglm")

res$symbol=mcols(dds)$gene_name
resLFCShrink_normal$symbol=mcols(dds)$gene_name
resLFCShrink_apeglm$symbol=mcols(dds)$gene_name

### with lfc threshold
if (is.null(args$lfc_threshold)){
  hk_genes = load_housekeeping()
  hk_thr = res[resLFCShrink_apeglm$symbol %in% hk_genes,"log2FoldChange"]
  lfc_thr = round(2*sd(hk_thr),2)
} else {
  lfc_thr = round(as.numeric(args$lfc_threshold),2)
}
cat(paste("\nlog2 Fold Change threshold:",lfc_thr,sep="\t"), file=paste(out_dir,"design.txt",sep="/"), append=TRUE)

res_lfc_cutoff = results(dds, lfcThreshold=lfc_thr, altHypothesis="greaterAbs")
resLFCShrink_cutoff_normal = lfcShrink(dds, lfcThreshold=lfc_thr, coef=tail(resultsNames(dds),n=1), type="normal")
resLFCShrink_cutoff_apeglm = lfcShrink(dds, lfcThreshold=lfc_thr, coef=tail(resultsNames(dds),n=1), type="apeglm")

print("res done")


res_lfc_cutoff$symbol=mcols(dds)$gene_name
resLFCShrink_cutoff_normal$symbol=mcols(dds)$gene_name
resLFCShrink_cutoff_apeglm$symbol=mcols(dds)$gene_name

print("add metadata done")

### write results to file ###
resOrdered = res[order(res$padj),]
write.table(as.data.frame(resOrdered),
            file = paste(out_dir,"DE_genes_basic.tsv",sep="/"),
            quote = F,
            sep = '\t')
print("done 1")

resOrdered_lfc = res_lfc_cutoff[order(res_lfc_cutoff$padj),]
write.table(as.data.frame(resOrdered_lfc),
            file = paste(out_dir,paste("DE_genes_lfcThr_",lfc_thr,".tsv",sep=""),sep="/"),
            quote = F,
            sep = '\t')
print("done 2")

resLFCShrink_normal_ordered = as.data.frame(resLFCShrink_normal[order(resLFCShrink_normal$padj),])
write.table(resLFCShrink_normal_ordered,
            file = paste(out_dir,"DE_genes_shrinked_normal.tsv",sep="/"),
            quote = F,
            sep = '\t')
print("done 3")

resLFCShrink_apeglm_ordered = as.data.frame(resLFCShrink_apeglm[order(resLFCShrink_apeglm$padj),])
write.table(resLFCShrink_apeglm_ordered,
            file = paste(out_dir,"DE_genes_shrinked_apeglm.tsv", sep="/"),
            quote = F,
            sep = '\t')
print("done 4")

resLFCShrink_normal_lfcThr_ordered = as.data.frame(resLFCShrink_cutoff_normal[order(resLFCShrink_cutoff_normal$padj),])
write.table(resLFCShrink_normal_lfcThr_ordered,
            file = paste(out_dir,paste("DE_genes_shrinked_lfcThr_",lfc_thr,"_normal.tsv",sep=""),sep="/"),
            quote = F,
            sep = '\t')
print("done 5")

resLFCShrink_apeglm_lfcThr_ordered = as.data.frame(resLFCShrink_cutoff_apeglm[order(resLFCShrink_cutoff_apeglm$svalue),])
#
write.table(resLFCShrink_apeglm_lfcThr_ordered,
            file = paste(out_dir,paste("DE_genes_shrinked_lfcThr_",lfc_thr,"_apeglm.tsv",sep=""), sep="/"),
            quote = F,
            sep = '\t')
print("done 6")

print("write to file done")

#combined_sig_res=na.omit(combined_sig_res)


### create df combining results and genes info ###
#combined_sig_res=as.data.frame(cbind(res_lfc_cutoff,mcols(dds)))
#combined_sig_res=combined_sig_res[!is.na(combined_sig_res$padj),]
#row.names(combined_sig_res) = combined_sig_res$gene_ens



### VISUALIZATION ###
#####################

### Volcano Plot ###
# drawLines_v <- function() abline(v=c(-0.585,0.585),col="dodgerblue",lwd=2)
#
# pdf("output_DESeq2_analysis/volcano_plot_withLFCcutoff.pdf")
# plot(res_lfc_cutoff$log2FoldChange, -log(res_lfc_cutoff$pvalue),pch=20, ylab='-Log10(pvalue)', xlab='Log2(FC)', main='Volcano Plot with LFC cut-off=1.5\nEPZ vs DMSO\nsignificant: padj<0.05')
# points(res_lfc_cutoff$log2FoldChange[res_lfc_cutoff$padj<0.05], -log(res_lfc_cutoff$pvalue[res_lfc_cutoff$padj<0.05]), col='red', pch=21)
# drawLines_v()
# dev.off()
#
# # pdf("output_DESeq2_analysis/volcano_plot_noLFCcutoff.pdf")
# # plot(res$log2FoldChange, -log(res$pvalue),pch=20, ylab='-Log10(pvalue)', xlab='Log2(FC)', main='Volcano Plot\nEPZ vs DMSO\nsignificant: padj<0.05')
# # points(res$log2FoldChange[res$padj<0.05], -log(res$pvalue[res$padj<0.05]), col='red', pch=21)
# # dev.off()
# # 
# # 
### MA-plots ###
drawLines <- function() abline(h=c(-lfc_thr,lfc_thr),col="dodgerblue",lwd=2)
pdf(paste(out_dir,"MAplot_NoLFCcutoff.pdf",sep="/"))
DESeq2::plotMA(res, alpha= 0.05)
#drawLines()
dev.off()

pdf(paste(out_dir,"MAplot_apeglm.pdf",sep="/"))
DESeq2::plotMA(resLFCShrink_apeglm, alpha= 0.05)
#drawLines()
dev.off()

pdf(paste(out_dir,"MA_plot_apeglm_LCFcutoff.pdf",sep="/"))
DESeq2::plotMA(resLFCShrink_cutoff_apeglm, alpha= 0.05)
drawLines()
dev.off()
# # 
# # pdf("output_DESeq2_analysis/MA_PlotMINE_NoLFCthr.pdf")
# # plot(log10(res$baseMean), res$log2FoldChange,
# #      pch=20, 
# #      main='MA plot\nEPZ vs DMSO\nsignificant: padj<0.05',
# #      ylab='Log2FC',
# #      xlab='Mean of normalized counts (Log10 scale)',
# #      ylim=c(-2,2))
# # points(log10(res$baseMean[res$padj < 0.05]), res$log2FoldChange[res$padj < 0.05], 
# #        col='red',
# #        pch=21,
# #        ylim=c(-2,2))
# # #drawLines()
# # abline(h=0, col='red', lwd=3)
# # dev.off()
# # 
# # pdf("output_DESeq2_analysis/MA_PlotMINE_LFCthr.pdf")
# # plot(log10(res_lfc_cutoff$baseMean), res_lfc_cutoff$log2FoldChange,
# #      pch=20, 
# #      main='MA plot\nEPZ vs DMSO\nsignificant: padj<0.05',
# #      ylab='Log2FC',
# #      xlab='Mean of normalized counts (Log10 scale)',
# #      ylim=c(-11,11))
# # points(log10(res_lfc_cutoff$baseMean[res_lfc_cutoff$padj < 0.05]), res_lfc_cutoff$log2FoldChange[res_lfc_cutoff$padj < 0.05], 
# #        col='red',
# #        pch=21,
# #        ylim=c(-11,11))
# # drawLines()
# # abline(h=0, col='red', lwd=3)
# # dev.off()
# # 
# # 
# # 
# # ####################################
# # ## FUNCTIONAL ENRICHMENT ANALYSIS ##
# # ####################################


res_no.na = na.omit(resLFCShrink_apeglm)
res_lfc_cutoff_no.na = na.omit(res_lfc_cutoff)

target_genes_noLFCthr = res_no.na$symbol[res_no.na$padj < 0.05]
target_genes_withLFCthr = res_lfc_cutoff_no.na$symbol[res_lfc_cutoff_no.na$padj < 0.05]
background_genes = res_no.na[res_no.na$baseMean>30,]$symbol

keytypes(org.Mm.eg.db)

### without LFC cutoff ###
eg = bitr(target_genes_noLFCthr, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
universe = bitr(background_genes, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
lfc = as.vector(res_no.na$log2FoldChange[match(eg$SYMBOL, res_no.na$symbol)])
names(lfc) = eg$ENTREZID

ego <- enrichGO(gene          = eg[["ENTREZID"]],
                universe      = universe[["ENTREZID"]],
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

pdf(paste(out_dir,"dotplot_DEgenes_noLFCthr.pdf", sep='/'), height = 8, width = 10)
dotplot(ego, showCategory = 20)
dev.off()

pdf(paste(out_dir,"barplot_DEgenes_noLFCthr.pdf",sep='/'), height = 8, width = 10)
barplot(ego, showCategory = 20)
dev.off()

# pdf(paste(out_dir,"emapplot_DEgenes_noLFCthr.pdf",sep='/'), height = 8, width = 10)
# emapplot(ego, showCategory = 30)
# dev.off()
# 
# pdf(paste(out_dir,"cnetplot_DEgenes_noLFCthr.pdf",sep='/'), height = 8, width = 10)
# cnetplot(ego, showCategory = 8, categorySize="pvalue", foldChange = lfc)
# dev.off()
# 
# pdf(paste(out_dir,"cnetplot_DEgenes_noLFCthr_20.pdf",sep='/'), height = 8, width = 10)
# cnetplot(ego, showCategory = 20, categorySize="pvalue", foldChange = lfc)
# dev.off()

### COMPARE CLUSTER UP vs DOWN
target_genes_noLFCthr_up = res_no.na$symbol[res_no.na$padj < 0.05 & res_no.na$log2FoldChange > 0]
eg_up = bitr(target_genes_noLFCthr_up, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
lfc_up = as.vector(res_no.na$log2FoldChange[match(eg_up$SYMBOL, res_no.na$symbol)])
names(lfc_up) = eg_up$ENTREZID


target_genes_noLFCthr_down = res_no.na$symbol[res_no.na$padj < 0.05 & res_no.na$log2FoldChange < 0]
eg_down = bitr(target_genes_noLFCthr_down, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
lfc_down = as.vector(res_no.na$log2FoldChange[match(eg_down$SYMBOL, res_no.na$symbol)])
names(lfc_down) = eg_down$ENTREZID




list_genes = list(genes_up = eg_up$ENTREZID, genes_down = eg_down$ENTREZID)

ck <- compareCluster(geneCluster = list_genes, 
                     fun = "enrichGO", 
                     OrgDb=org.Mm.eg.db, 
                     ont="BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

head(as.data.frame(ck))
write.csv(as.data.frame(ck), paste(out_dir,"Up_vs_Down_Enrichment.tsv",sep='/'), quote=F, sep="\t")

pdf(paste(out_dir,"Up_vs_Down_Enrichment.pdf",sep="/"), height = 8, width = 10)
dotplot(ck, showCategory = 10)
dev.off()


### with LFC cutoff ###
# eg_withLFCthr = bitr(target_genes_withLFCthr, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
# lfc_withLFCthr = as.vector(res_lfc_cutoff_no.na$log2FoldChange[match(eg_withLFCthr$SYMBOL, res_lfc_cutoff_no.na$symbol)])
# names(lfc_withLFCthr) = eg_withLFCthr$ENTREZID
# 
# ego <- enrichGO(gene          = eg_withLFCthr[["ENTREZID"]],
#                 universe      = universe[["ENTREZID"]],
#                 OrgDb         = org.Mm.eg.db,
#                 ont           = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.05,
#                 qvalueCutoff  = 0.1,
#                 readable      = TRUE)
# 
# #pdf("output_DESeq2_analysis/dotplot_DEgenes_withLFCthr.pdf")
# dotplot(ego, showCategory = 20)
# #dev.off()
# #pdf("output_DESeq2_analysis/barplot_DEgenes_withLFCthr.pdf")
# barplot(ego, showCategory = 20)
# #dev.off()
# #pdf("output_DESeq2_analysis/emapplot_DEgenes_withLFCthr.pdf")
# emapplot(ego, showCategory = 50)
# #dev.off()
# #pdf("output_DESeq2_analysis/cnetplot_DEgenes_withLFCthr.pdf", width = 12, height = 10)
# cnetplot(ego, showCategory = 8, categorySize="pvalue", foldChange = lfc)
# #dev.off()

### only UPREGULATED genes ###
# target_genes_noLFCthr_up = res_no.na$symbol[res_no.na$padj < 0.01 & res_no.na$log2FoldChange > 0]
# eg_up = bitr(target_genes_noLFCthr_up, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
# 
# lfc_up = as.vector(res_no.na$log2FoldChange[match(eg_up$SYMBOL, res_no.na$symbol)])
# names(lfc_up) = eg_up$ENTREZID
# 
# ego_up <- enrichGO(gene          = eg_up[["ENTREZID"]],
#                 universe      = universe[["ENTREZID"]],
#                 OrgDb         = org.Mm.eg.db,
#                 ont           = "BP",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.05,
#                 qvalueCutoff  = 0.1,
#                 readable      = TRUE)
# 
# 
# pdf("output_DESeq2_analysis/enrichment_analysis_upregulatedGenes/dotplot_DEgenes_noLFCthr_up.pdf")
# dotplot(ego_up, showCategory = 20)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_upregulatedGenes/barplot_DEgenes_noLFCthr_up.pdf")
# barplot(ego_up, showCategory = 20)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_upregulatedGenes/emapplot_DEgenes_noLFCthr_up.pdf")
# emapplot(ego_up, showCategory = 30)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_upregulatedGenes/cnetplot_DEgenes_noLFCthr_up.pdf", width = 10, height = 8)
# cnetplot(ego_up, showCategory = 8, categorySize="pvalue", foldChange = lfc_up)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_upregulatedGenes/cnetplot_DEgenes_noLFCthr_20_up.pdf", width = 20, height = 20)
# cnetplot(ego_up, showCategory = 30, categorySize="pvalue", foldChange = lfc_up)
# dev.off()


### only DOWNREGULATED genes ###

# target_genes_noLFCthr_down = res_no.na$symbol[res_no.na$padj < 0.01 & res_no.na$log2FoldChange < 0]
# eg_down = bitr(target_genes_noLFCthr_down, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
# 
# lfc_down = as.vector(res_no.na$log2FoldChange[match(eg_down$SYMBOL, res_no.na$symbol)])
# names(lfc_down) = eg_down$ENTREZID
# 
# ego_down <- enrichGO(gene          = eg_down[["ENTREZID"]],
#                    universe      = universe[["ENTREZID"]],
#                    OrgDb         = org.Mm.eg.db,
#                    ont           = "BP",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.05,
#                    qvalueCutoff  = 0.1,
#                    readable      = TRUE)
# 
# 
# pdf("output_DESeq2_analysis/enrichment_analysis_downregulatedGenes/dotplot_DEgenes_noLFCthr_down.pdf")
# dotplot(ego_down, showCategory = 20)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_downregulatedGenes/barplot_DEgenes_noLFCthr_down.pdf")
# barplot(ego_down, showCategory = 20)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_downregulatedGenes/emapplot_DEgenes_noLFCthr_down.pdf", height = 10, width = 10)
# emapplot(ego_down, showCategory = 30)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_downregulatedGenes/cnetplot_DEgenes_noLFCthr_down.pdf", width = 10, height = 8)
# cnetplot(ego_down, showCategory = 8, categorySize="pvalue", foldChange = lfc_down)
# dev.off()
# pdf("output_DESeq2_analysis/enrichment_analysis_downregulatedGenes/cnetplot_DEgenes_noLFCthr_50_down.pdf", width = 20, height = 20)
# cnetplot(ego_down, showCategory = 30, categorySize="pvalue", foldChange = lfc_down)
# dev.off()
# # 
# # 
# # 
# # 
# # ### kegg enrichment analysis
# # kk <- enrichKEGG(gene         = eg[["ENTREZID"]],
# #                  organism     = 'mmu',
# #                  universe      = universe[["ENTREZID"]],
# #                  pvalueCutoff = 0.05,
# #                  qvalueCutoff = 0.1)
# # 
# # dotplot(kk)
# # 
# # library("pathview")
# # #my_path <- pathview(gene.data  = eg[["ENTREZID"]],
# # #                     pathway.id = 'mmu03040',
# # #                     species    = "mmu")
# # 
# # 
# # 
# # #########################
# # ### create RPKM table ###
# # #########################
# # 
# # 
# # library(edgeR)
# # 
# # lengths = read.csv("../featureCounts/Vogel_mDorsalTelencefalon-Ctrl-E14.5_RNA-Seq_Rep1.counts.txt", sep="\t", header = F, skip = 2)
# # rownames(lengths) = lengths$V1
# # lengths[,c(1,2,3,4,5,7)] = NULL
# # colnames(lengths) = c("length")
# # 
# # index_2=match(row.names(count),row.names(lengths))
# # lengths_def = lengths[index_2,]
# # 
# # 
# # 
# # y=DGEList(count=count, genes=data.frame(length=lengths_def))
# # y = calcNormFactors(y)
# # RPKM = as.data.frame(rpkm(y))
# # 
# # #RPKM$Average = apply(RPKM, 1, mean)
# # 
# # write.table(RPKM,
# #             'output_DESeq2_analysis/RPKM_expression.tsv',
# #             sep="\t",
# #             quote=F)
# # 
# # 
# # 
# # ### GAGE gene set enrichment analysis ###
# # #########################################
# # 
# # featureData_1 = featureData
# # RPKM_1 = RPKM
# # 
# # featureData_1 = featureData_1[match(unique(featureData_1$gene_name),featureData_1$gene_name),]
# # RPKM_1 = RPKM_1[row.names(featureData_1),]
# # 
# # RPKM_1 = RPKM_1[!row.names(RPKM_1) %in% dup_ID,]
# # 
# # index_4 = match(row.names(RPKM_1), row.names(featureData_1))
# # 
# # rownames(RPKM_1) = toupper(featureData_1$gene_name[index_4])
# # 
# # data("mm_GO")
# # 
# # gage_res = gage(RPKM_1, gsets = mm_GO, ref = c(1,2,3,4,5), same.dir = F)
# # 
# # greater=as.data.frame(head(as.data.frame(gage_res$greater)[as.data.frame(gage_res$greater)$q.val<0.1,],n=10))
# # greater$logtr_q.val=-log10(greater$q.val)
# # greater$gene_sets=rownames(greater)
# # greater=greater[order(greater$logtr_q.val),]
# # 
# # lower=as.data.frame(gage_res$less)[as.data.frame(gage_res$less)$q.val<0.1,]
# # lower$logtr_q.val=-log10(lower$q.val)
# # lower$gene_sets=rownames(lower)
# # lower=lower[order(lower$logtr_q.val),]
# # 
# # figure_title="GAGE gene sets enrichment analysis\nDb2 vs C1+C2"
# # 
# # gr=ggplot(data=greater, aes(x=reorder(gene_sets,logtr_q.val) , y=logtr_q.val, col="upregulated\ngene sets")) +
# #   geom_bar(stat="identity", col='red')+
# #  #scale_fill_manual(values=c("red"))+
# #   theme_minimal()+
# #   theme(plot.title = element_text(hjust = 0.5),  axis.title.x=element_blank(),axis.title.y =element_blank())
# # gr= gr + coord_flip(ylim = c(0, 4))
# # 
# # lw=ggplot(data=lower, aes(x=reorder(gene_sets,logtr_q.val), y=logtr_q.val)) +
# #   geom_bar(stat="identity",col='blue')+
# #   #scale_fill_manual(values=c("blue"))+
# #   theme_minimal()+
# #   theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.title.y =element_blank())
# # lw=lw + coord_flip(ylim = c(0, 4))
# # 
# # fig=ggarrange(gr,lw, #labels = c('UP-REGULATION','DOWN-REGULATION'),
# #           ncol = 1, nrow = 2, align = "hv",
# #           common.legend = TRUE)
# # fig=annotate_figure(fig,
# #                 top = text_grob(figure_title, color = "black", face = "bold", size = 10),
# #                 bottom = text_grob("-Log10(q value)", color = "black",
# #                                     hjust = 1, x = 1, face = "italic", size = 10),
# #                # left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
# #                 left  = "Significantly enriched canonical pathways(q.val < 0.1)"
# #   ) 
# # 
# # pdf(paste(path,'results_DESeq2/mRNA/Db2_vs_C1C2/GAGE_C2GO_PATHWASONLY_Db2_vs_C1C2_mRNA.pdf',sep=''), width = 13)
# # fig
# # dev.off()
# # 
# # 
# # 
# # ### PCA manually with prcomp ###
# # 
# # matrix_rlog=as.data.frame(t(assay(rld)))
# # variances=apply(matrix_rlog, 2, var)
# # variances=variances[order(variances, decreasing = T)]
# # matrix_rlog_1 = matrix_rlog[,names(matrix_rlog) %in% names(variances)[1:500]]
# # 
# # pca=prcomp(matrix_rlog)
# # pca_1=prcomp(matrix_rlog_1)
# # summary_pca=(summary(pca))
# # summary_pca_1=summary(pca_1)
# # plot(pca$rotation[,2]~pca$rotation[,1])
# # pdf('PCA_2vs3.pdf', height= 4, width = 8)
# # par(mfrow=c(1,2))
# # plot(pca_1$x[,2]~pca_1$x[,1], col=c(rep('blue',4),rep('red',4)), pch=20, xlab='PC1: ~ 87% variance', ylab='PC2: ~ 6% variance'); legend("topleft", legend=c("db/bd", "Controls"), col=c(2,4), pch=20)
# # plot(pca_1$x[,3]~pca_1$x[,2], col=c(rep('blue',4),rep('red',4)), pch=20, xlab='PC2: ~ 6% variance', ylab='PC3: ~ 3% variance'); legend("topright", legend=c("db/bd", "Controls"), col=c(2,4), pch=20)
# # dev.off()
# # 
# # loadings <- data.frame(gene_ens=rownames(pca$rotation), pca$rotation[,1:2]^2, stringsAsFactors=F)
# # loadings <- loadings[order(loadings$PC1, decreasing=T),]
# # loadings_def = merge(loadings, featureData, by='gene_ens')
# # loadings_def = loadings_def[order(loadings_def$PC1, decreasing = T),]
# # top_var_genes=as.character(head(unique(loadings_def$gene_name),n=1000))
# # write(top_var_genes, 'top_var_genes_PC1.txt', sep='\t')
# # #write(as.character(unique(loadings_def$gene_name)), 'background_genes.txt', sep='\t')
# # 
# # 
# # rownames(loadings) <- NULL
# # head(loadings)
# # png("screeplot.png", w=600, h=600)
# # screeplot(pca, type="lines", pch=16, main="Scree plot")
# # dev.off()
# # png("loadings.png", w=800, h=500)
# # par(mfrow=c(1,2))
# # plot(sort(loadings$PC1, decr=T), pch=20, main="PC1 loadings", ylab='loadings', xlab = 'gene index'); points(sort(loadings$PC1, decr=T)[1:1000], col=2, pch=20); legend("topright", legend=c("Top 1000 Genes", "Other Genes"), col=c(2,1), pch=20)
# # plot(sort(loadings$PC2, decr=T), pch=20, main="PC2 loadings", ylab='loadings', xlab = 'gene index'); points(sort(loadings$PC2, decr=T)[1:1000], col=2, pch=20); legend("topright", legend=c("Top 1000 Genes", "Other Genes"), col=c(2,1), pch=20)
# # dev.off()
# # par(mfrow=c(1,1))
# # 
# # biplot(pca)
# # summary(pca)
# # 
# # 
