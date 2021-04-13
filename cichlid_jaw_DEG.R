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


#set your working directory
setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression") #you will need to change to your own directory

# import sample table
sampleTable <- read.table("paper1_speciesinfo2.txt", header=T)

#Make counts matrix
files <- list.files(path="/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression")
genes <- read.table(files[1], header=FALSE, sep="\t")[,1]     # gene names
countData    <- do.call(cbind,lapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))
countData    <- cbind(genes,df)
colnames(countData) <- countData[1, ]
countData <- countData[-1, ]
countData<-as.data.frame(countData)

#Make colData
colData<-sampleTable



