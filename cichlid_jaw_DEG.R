source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")



###conduct array quality metrics to detect and remove outliers

#install.packages("BiocManager")
library("BiocManager")
BiocManager::install("genefilter")
BiocManager::install("affycoretools")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("Biobase")
BiocManager::install("rnaseqWrapper")
library(DESeq2) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)


#biocLite("BiocParallel")
library("BiocParallel", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library/")

sessionInfo()
packageVersion("DESeq2")

BPPARAM = MulticoreParam(workers=5)

#set your working directory
setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression") #you will need to change to your own directory

# reads files
countFiles = list.files("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression", pattern="counts1$", full.names=TRUE)
countFiles

directory <- "/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression"
sampleFiles <- grep(".counts1",list.files(directory),value=TRUE)


# import sample table
sampleTable <- read.table("paper1_speciesinfo2.txt", header=T)

#start analysis
#dds = DESeqDataSetFromHTSeqCount(countFiles, sampleData=sampleTable, design= ~ condition)

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 5, ]

#run analysis

dds <- DESeq(dds) 

#res <- results(dds)
#write.table(res, file = "DESeq_allsamples_fullresults.txt", sep = "\t", quote=F)

#save data
save(dds, file="./DESeq2_all.RData")

#plot dispersion plot
png("dispersion_plot.png")
plotDispEsts(dds)
dev.off()

#matrix of normalized counts

rld <- rlog(dds) 
vsd <- varianceStabilizingTransformation(dds)
matrix <- (assay(rld))
write.table(matrix, file = "DESeq_NormalizedCounts.txt", sep = "\t", quote=F)

#########################################################################

#plot SD of transformed data using shifted logarithm transformation

library("vsn") 
notAllZero <- (rowSums(counts(dds))>5)


png("DESeq_meanSDplot_normalized_shiftedlog.png") 
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
dev.off()

png("DESeq_meanSDplot_normalized_regularizedlog.png")
meanSdPlot(assay(rld[notAllZero,]))
dev.off()

png("DESeq_meanSDplot_normalized_varianceStabilizingTransformation.png")
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

##########################################################################

#heatmap of normalized counts

library("pheatmap") 
#select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
#nt <- normTransform(dds)
#log2.norm.counts <- assay(nt)[select,] 
#df <- as.data.frame(colData(dds)[,("condition")]) 

#png("heatmap_20_highest_expr_genes_log2normalized.png")
#pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#dev.off()

#png("heatmap_20_highest_expr_genes_rld.png")
#pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#dev.off()

#png("heatmap_20_highest_expr_genes_vsd.png")
#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#dev.off()

#########################################################################

#make contrasts

res_LdTd <- results(dds, contrast=c("condition", "Ld", "Td"), alpha=0.05)
write.table(res_LdTd, file = "LdTd_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)

########################################################################

#summarise results

res_LdTd1 <- results(dds, contrast=c("condition", "Ld", "Td"), alpha=0.01)
write.table(res_LdTd1, file = "LdTd_DESeq_fullresults_0.01.txt", sep = "\t", quote=F)

#summary of results

#sum <- summary(res_LdTd)
#write.table(sum, file = "LdTd_DESeq_summary_0.05.txt", sep = "\t", quote=F)

genesDE <- sum(res_LdTd$padj < 0.05, na.rm=TRUE)
write.table(genesDE, file = "LdTd_DESeq_genesDE_0.05.txt", sep = "\t", quote=F)

#sum1 <- summary(res_LdTd1)
#write.table(sum1, file = "LdTd_DESeq_summary_0.01.txt", sep = "\t", quote=F)

genesDE1 <- sum(res_LdTd1$padj < 0.01, na.rm=TRUE)
write.table(genesDE1, file = "LdTd_DESeq_genesDE_0.01.txt", sep = "\t", quote=F)

##########################################################################

#plot MA plot

png("DEseq_MAplot_0.05.png")
plotMA(res_LdTd, main="DESeq2_0.05", ylim=c(-2,2))
dev.off()

png("DEseq_MAplot_0.01.png")
plotMA(res_LdTd1, main="DESeq2_0.01", ylim=c(-2,2))
dev.off()

##########################################################################

#heatmaps of sample to samples distances: clustering

png("DEseq_clustering_heatmap_ofsamples_rld.png")
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer") 
sampleDistMatrix <- as.matrix(sampleDists) 
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
colnames(sampleDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=colors)
dev.off()


png("DEseq_clustering_heatmap_ofsamples_vsd.png")
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer") 
sampleDistMatrix <- as.matrix(sampleDists) 
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
colnames(sampleDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=colors)
dev.off()

#######################################################################


#principal component plot of samples

png("PCA_samples_rld.png")
plotPCA(rld, intgroup=("condition"))
dev.off()

png("PCA_samples_vsd.png")
plotPCA(vsd, intgroup=("condition"))
dev.off()

#netapp/home/tmp/RtmpAeSB6R/downloaded_packages



