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
library(ggplot2)


#biocLite("BiocParallel")
library("BiocParallel", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library/")

sessionInfo()
packageVersion("DESeq2")


#set your working directory
setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression") #you will need to change to your own directory

# import sample table
sampleTable <- read.table("paper1_speciesinfo2.txt", header=T)

#Make counts matrix
files <- list.files(path="/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression")
genes <- read.table(files[1], header=FALSE, sep="\t")[,1]     # gene names
countData    <- do.call(cbind,lapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))
countData    <- cbind(genes,countData)
colnames(countData) <- countData[1, ]
countData <- countData[-1, ]
data_num <- as.data.frame(apply(countData, 2, as.numeric))
df <- data_num[,-1]
rownames(df) <- countData[,1]
rm(countData)
rm(data_num)
countData<-df

v=setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw")


#Make colData
colData<-sampleTable
g<-colData

dds=DESeqDataSetFromMatrix(countData=countData,
                           colData = g,
                           design = ~condition)

vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("condition"),force=T)

###### Deseq #######

totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, ylab="raw counts")

dds<-DESeqDataSetFromMatrix(countData=countData, colData=g, design=~condition)

#One-step DESeq
dds<-DESeq(dds)

head(dds)
res<- results(dds)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot Snails")
g






rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$condition)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valstress)
head(rldpvals)
dim(rldpvals)
# [1] 19717    13
table(complete.cases(rldpvals))
#FALSE  TRUE 
#17202  2515




####### Morphometrics #######

library(geomorph)
library(tidyverse)
library(ggplot2)
library(factoextra)

morpho <- readland.tps("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/GBE_final_landmarks.tps", specID = "imageID")


morpho2D<-two.d.array(morpho)
class<-read.csv("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/morpho_classifiers.csv")
class$Species<-as.factor(class$Species)
rownames(morpho2D) <- class$Species



#Procrustes superimposition
morpho.gpa<-gpagen(morpho, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                   max.iter = 30, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)

coords<-morpho.gpa$coords # a 3D array of Procrustes coordinates
coords2D<-two.d.array(coords) # a 2D array of Procrustes coordinates
rownames(coords2D)<-class$Species
t<-as.data.frame(coords2D)

morpho.gpa$Csize # a vector of centroid sizes

pca<-prcomp(t)



#Extract PC scores for significant PCs
pc_scores<-pca$x[,1:2]
pc_scores<-as.data.frame(pc_scores)
pc_scores$Species<-class$Species


fviz_pca_biplot(pca, label = "ind",
                col.var = "black",
                geom = "point",
                pointsize=4,
                col.ind=class$Species,
                mean.point=FALSE,
                repel = TRUE,
                addEllipses = FALSE,
                labelsize=3,
                alpha.var=0.5,
                arrowsize=0.25,
                invisible = "var")
