## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "GSE107741/data")
resultsDir <- file.path(workingDir, "GSE107741/results")


if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}

installifnot("pd.mogene.1.0.st.v1")
installifnot("mogene10sttranscriptcluster.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")


## ------------------------------------------------------------------------
targets <- read.delim("~/Documents/BDBI/2n/3r_trim/Omics_Techniques/Camps_Judit-OmicsTechniques/Exercise_3/GSE107741/data/targets.txt")
class(targets)
print(targets)


## ------------------------------------------------------------------------
expression <- as.matrix(read.delim("~/Documents/BDBI/2n/3r_trim/Omics_Techniques/Camps_Judit-OmicsTechniques/Exercise_3/GSE107741/data/expression.txt", row.names=1, comment.char="#"))
class(expression)
head(expression)


## ------------------------------------------------------------------------
if (!require(GEOquery)) {
  BiocManager::install("GEOquery")
}
require(GEOquery)
gse <- getGEO("GSE107741")
class(gse)
names(gse)

rawData <- gse[[1]]

names(pData(rawData))


## ------------------------------------------------------------------------
colnames(rawData) <-rownames(pData(rawData)) <- targets$SampleName
sampleNames <- as.character(targets$SampleName)
sampleColor <- as.character(targets$Colors)


## ------------------------------------------------------------------------
#BOXPLOT
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)



## ------------------------------------------------------------------------
#HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
#PRINCIPAL COMPONENT ANALYSIS
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=targets$sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)


## ------------------------------------------------------------------------
# eset_rma <- rma(rawData)


## ------------------------------------------------------------------------
BiocManager::install("affydata")
eset <- normalize(rawData)


## ------------------------------------------------------------------------
#BOXPLOT
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
#HIERARQUICAL CLUSTERING

clust.euclid.average <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
#PRINCIPAL COMPONENT ANALYSIS

plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

## ------------------------------------------------------------------------
plotPCA(exprs(eset), labels= targets$sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)


## ------------------------------------------------------------------------

require (limma)


## ------------------------------------------------------------------------
#CONTRAST MATRIX.lINEAR MODEL
treat <- targets$Type
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design)<-c("Fertile", "Subfertile")
rownames(design) <- targets$SampleName 
print(design)

#COMPARISON
cont.matrix1 <- makeContrasts (
  Fertile.vs.Subfertile = Fertile-Subfertile,
  levels=design)
comparison1 <- "Effect of Induction"


## ------------------------------------------------------------------------
#MODEL FIT
fit1 <- lmFit(eset, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1 <- eBayes(fit.main1)


## ------------------------------------------------------------------------
#FILTER BY FALSE DISCOVERY RATE AND FOLD CHANGE
topTab <-  topTable (fit.main1, number=nrow(fit.main1), coef="Fertile.vs.Subfertile", adjust="fdr",lfc=abs(3))


## ------------------------------------------------------------------------
volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep="\n"))


