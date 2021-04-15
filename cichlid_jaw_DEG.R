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
setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw") #you will need to change to your own directory

# import sample table
sampleTable <- read.table("paper1_speciesinfo2.txt", header=T)

#Make counts matrix
files <- list.files(path="/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression")
setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw/gene_expression")
genes <- read.table(files[1], header=FALSE, sep="\t")[,1]     # gene names
countData<-do.call(cbind,lapply(files,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))
countData<-cbind(genes,countData)
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
### Outliers ####
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

setwd("/Users/13018/Documents/Ecological_Genomics/cichlid_jaw")
write.csv(countData, file="counts_WGCNA.csv", quote=F, row.names=TRUE)





###### WGCNA ######

dds<-DESeqDataSetFromMatrix(countData=countData, colData=g, design=~ condition) 
dds<-DESeq(dds)

res<- results(dds)
#filter for contigs with average(baseMean) >3
res3<-res[res$baseMean>3, ]
dim(res) #26993
dim(res3) #18257

# get rlog data (better transformation when size factors vary across samples)
rld <- rlogTransformation(dds, blind=TRUE, fitType="local")
head(rld)
rld_wg=(assay(rld))
head(rld_wg)
nrow(rld_wg)
#26993
rldFiltered=(assay(rld))[(rownames((assay(rld))) %in% rownames(res3)),]
nrow(rldFiltered)
#18257
write.csv( rldFiltered,file="jaw_wgcna_allgenes.csv",quote=F,row.names=T)
#now we have our filtered data to take into WGCNA

####First part of tutorial:Data input and cleaning
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

dat=read.csv("jaw_wgcna_allgenes.csv")
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#18257
datExpr0 = as.data.frame(t(dat))

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #if TRUE, no outlier genes, if false run the script below

### Outlier detection incorporated into trait measures. 
traitData= read.csv("WGCNA_traits.csv", row.names=1)
dim(traitData)
head(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
dim(datExpr0)
rownames(datExpr0)
rownames(traitData)=rownames(datExpr0)
traitData$Sample= NULL 
# datTraits=allTraits
datTraits=traitData

table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="signed")
# this calculates the whole network connectivity we choose signed because we care about direction of gene expression
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

save(datExpr0, datTraits, file="jaw_Samples_Traits_ALL.RData")

#Figure out proper SFT
# Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5)); #may need to adjust these power values to hone in on proper sft value
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9 (but still under)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=15 #smallest value to plateau at ~0.85
adjacency=adjacency(datExpr0, power=softPower,type="signed") #must change method type here too!!
#heatmap(adjacency, labRow=FALSE, labCol=FALSE)

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

save(TOM,file="TOM.RDS")

library(flashClust)
geneTree= flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(10,6)
# pdf(file="dendrogram_thresh16.5_signed_1868.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
# dev.off()
#each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes

minModuleSize=90 #we only want large modules
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

dynamicColors= labels2colors(dynamicMods)
#plot dendrogram and colors underneath, pretty sweet
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 15)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_crep_nomerge.RData")

lnames = load(file = "Network_crep_nomerge.RData")
#plot
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

MEDissThres= 0.38
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_signed_0.6.RData")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "jaw_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network_signed_0.6.RData");
lnames = load(file = "Network_crep_nomerge.RData");
lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#represent module trait correlations as a heatmap
quartz()
sizeGrWindow(12,8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#Gene relationship to trait and important modules:
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraits$Tooth); #change Lipidrobust to your trait name
names(weight) = "Tooth"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

#Gene-trait significance correlation plots
# par(mfrow=c(2,3))
module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for Tooth Shape",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#names(dis)
sizeGrWindow(8,7);
which.module="black" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

#quartz()
# par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")
#this is a cool plot where you can see that genes in this module are upregulated in the pH7.5 treatment

##############################heatmap of module expression with bar plot of trait of interest by sample...
#here we just have binary traits, but if you have a continuous trait this code is cool
sizeGrWindow(8,7);
which.module="black" #pick module of interest
which.trait="Tooth" #change trait of interest here
datTraits=datTraits[order((datTraits$PC1),decreasing=T),]#change trait of interest here

trait=datTraits[, paste(which.trait)]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset
genes=genes[rownames(datTraits),]

#quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(trait, col=which.module, main="", cex.main=2,
        ylab="fvfm",xlab="sample")#change trait of interest here

#Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="black"]) #black  blue brown green  grey  pink   red 

c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size
table(moduleColors)
#moduleColors
#black         blue       coral2     darkgrey  floralwhite       grey60        ivory 
#1338         2928         1167         1088          895          206          564 
#mediumorchid   orangered4       purple    steelblue 
#478         1557          615         1749
head(c.vsd)
write.csv(c.vsd,"rlog_black.csv",quote=F)

#Gene relationship to trait and important modules: Gene Significance and Module membership
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_black.csv", row.names=1)
head(vsd)
gg=read.table("Crep454_iso2gene.tab", sep="\t")
head(gg)
library(pheatmap)

############################################
whichModule="ivory"
top=100

datME=MEs
vsd <- read.csv("Crep_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]
# attaching gene names
summary(hubs)

gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$V1) { 
    counts=counts+1
    gn=gg[gg$V1==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn) 
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames
length(hubs)

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#quartz()
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))



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
                invisible = "var") + labs(title="Shape Variation of the Jaw Apparatus",x="PC1 (48.4%)",y="PC2 (11.9%)")

