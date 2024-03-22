#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

# loading functions

# from system get the value of the environment variabl of export AlsammanPopulationGeneticsPath
functionsPath<-Sys.getenv("AlsammanPopulationGeneticsPath")
functionsPath<-paste(functionsPath,"DiversityAnalysis/phylogeneticAnalysis/functions",sep="")
# read all r source scripts in function folder
sourcesList<-list.files(paste(functionsPath,sep=""),pattern="*.r",full.names=TRUE)
# if the list is empty then exit
if (length(sourcesList)==0) {
  print(paste("No functions found in you should use the configureMe.sh script to set the path to the functions folder",sep=""))
  quit()
}
# source all functions
for (i in 1:length(sourcesList)) {
  source(sourcesList[i])
}






#####################################################################################################
# What is this system linux or windows
system<-Sys.info()[1]
# if system is windows then change the path separator
if (system=="Windows") {
  pathSeparator<-"\\"
} else {
  pathSeparator<-"/"
}

#######################################################################################


args = commandArgs(trailingOnly=TRUE)


if (length(args) < 3) {
  print(paste("Usage: Rscript phylogenetic_analysis.r <vcffile> <ResultsFolder> <maf> <missing.rate>",sep=" "))
  quit()
}

#Loading libraries
library(SNPRelate)
library(ggplot2)
library(ape)
library(RColorBrewer)



vcffile<-args[1]
metaFile<-args[2]
ResultsFolder<-args[3]
outFolder<-args[4]

maf<-0.05
missing.rate<-0.2

if (length(args) > 4) {
  maf<-as.numeric(args[5])
  missing.rate<-as.numeric(args[6])
} else {
  print(paste("Using default values for maf  > ",maf," and missing.rate > ",missing.rate,sep=" "))
}


vcffile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/BreadWheatBefore_Filtered.vcf"
metaFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/passportData.order.representative.tsv"
ResultsFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Results"

ResultsFolder<-paste(ResultsFolder,"/phylogeneticAnalysis",sep="")
# create the folder if not exists
if (!file.exists(ResultsFolder)) {
  dir.create(ResultsFolder, recursive = TRUE)
}

# get the name of the file
outfilename<-basename(vcffile)
outputTitle<-"Phylogenetic tree"

gdsfilename <- paste(outfilename,".gds",sep="")
# remove file if exists
if (file.exists(gdsfilename)) {
    file.remove(gdsfilename)
}

snpgdsVCF2GDS(vcffile, gdsfilename)
snpgdsSummary(gdsfilename)
genofile <- openfn.gds(gdsfilename)

### Calculate the distance matrix and phylogenetic tree
#Calculate Distance matrix
dissMatrix  <-  snpgdsDiss(genofile, maf=maf, missing.rate= missing.rate, autosome.only = FALSE, num.thread=4, verbose=TRUE)

# Hierarchical clustering
snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)

# Phylogenetic tree
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, 
    col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, 
    verbose=TRUE)
# close genotype file
closefn.gds(genofile)


pdf(paste(ResultsFolder,"/",outfilename,".pdf",sep=""), width=30, height=8)
snpgdsDrawTree(
    cutTree, 
    main = "Phylogenetic tree",
    y.label.kinship=T,
    cex=1, 
    font=1, 
    cex.lab=0.5,
    lwd=10,
    cex.axis=2)
dev.off()

# as a png  
png(paste(ResultsFolder,"/",outfilename,".png",sep=""), width=30, height=8, units="in", res=300)
snpgdsDrawTree(
    cutTree, 
    main = "Phylogenetic tree",
    y.label.kinship=T,
    cex=1, 
    font=1, 
    cex.lab=0.5,
    lwd=10,
    cex.axis=2)
dev.off()


write.tree(as.phylo(as.hclust(cutTree$dendrogram)), file=paste(ResultsFolder,"/",outfilename,".newick",sep=""))

# samples 
# save samples as text file
write.table(cutTree$samples, file=paste(ResultsFolder,"/",outfilename,".samples.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
# save samp.order
write.table(cutTree$samp.order, file=paste(ResultsFolder,"/",outfilename,".samp.order.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)
# save samp.group
write.table(cutTree$samp.group, file=paste(ResultsFolder,"/",outfilename,".samp.group.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F)

# save the distance matrix
dissMatrixTable<-as.data.frame(dissMatrix$diss)
rownames(dissMatrixTable)<-dissMatrix$sample.id
colnames(dissMatrixTable)<-dissMatrix$sample.id
write.table(dissMatrixTable, file=paste(ResultsFolder,"/",outfilename,".dissMatrix.tsv",sep=""), sep="\t", quote=F, row.names=T, col.names=T)

# stats of the phylogenetic tree
TreeSummary<-summary(cutTree)
# save the stats
write.table(TreeSummary, file=paste(ResultsFolder,"/",outfilename,".TreeSummary.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T)

# read the metadata file
metaData<-read.table(metaFile, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F, row.names=1)
samplesNames<-rownames(metaData)

# get the the values of colors in hex
# 2 hex colors
colors<-brewer.pal(12, "Set3")

# all the columns to factor except the first one
metaData<-data.frame(lapply(metaData, as.factor))
# add row names
rownames(metaData)<-samplesNames


create_iTOL_annotation_file<-function(metaData, metCol, outfilename, ResultsFolder) {
  # create the iTOL annotation file
# DATASET_COLORSTRIP
# SEPARATOR SPACE
# DATASET_LABEL Species
# COLOR #ff0000
# COLOR_BRANCHES 0
# LEGEND_TITLE Species
# STRIP_WIDTH 30
# LEGEND_SHAPES 1 1 1 1 1 1
# LEGEND_COLORS #E41A1C #377EB8 #4DAF4A #984EA3 #FF7F00 #FFFF33
# LEGEND_LABELS Unknown Tv Ta_._aestivum Ta_._compactum Av Tt_._durum

  sink(paste(ResultsFolder,"/",outfilename,".colors.txt",sep=""))
  cat("DATASET_COLORSTRIP\n")
  cat("SEPARATOR TAB\n")
  cat("DATASET_LABEL\t",colnames(metaData)[metCol],"\n",sep="")
  cat("COLOR\t#ff0000\n")
  cat("COLOR_BRANCHES\t0\n")
  cat("LEGEND_TITLE\t",colnames(metaData)[metCol],"\n",sep="")
  cat("STRIP_WIDTH\t30\n")
  shapes<-rep(1,length(levels(metaData[,metCol])))
  shapes<-paste(shapes,collapse="\t")
  cat("LEGEND_SHAPES\t",shapes,"\n",sep="")
  cat("LEGEND_COLORS\t",paste(colors[1:length(levels(metaData[,metCol]))],collapse="\t"),"\n",sep="")
  cat("LEGEND_LABELS\t",paste(levels(metaData[,metCol]),collapse="\t"),"\n",sep="")


  cat("DATA\n")
  for (i in 1:nrow(metaData)) {
    cat(rownames(metaData)[i],"\t",colors[as.numeric(metaData[i,metCol])],"\t", as.character(metaData[i,metCol]),"\n",sep="")
  }
  sink()

  # create a range strip
  sink(paste(ResultsFolder,"/",outfilename,".range.txt",sep=""))
  cat("TREE_COLORS\n")
  cat("SEPARATOR TAB\n")
  cat("DATA\n")
  for (i in 1:nrow(metaData)) {
    cat(rownames(metaData)[i],"\t","range","\t",colors[as.numeric(metaData[i,metCol])],"\t", as.character(metaData[i,metCol]),"\n",sep="")
  }
  sink()


}

for (i in 1:ncol(metaData)) {
  
  create_iTOL_annotation_file(metaData, i, colnames(metaData)[i], ResultsFolder)
}


# calculate the sub distance matrix in each group
library(reshape2)
getDistanceData<-function(metCol,metaData,dissMatrixTable)
{
  categories<-levels(metaData[,metCol])
  # create a list of the categories
  distCatList<-list()
  for (i in 1:length(categories)) {
    cat("Calculating the distance matrix for category ",categories[i],"\n")
    # get the samples in the category
    samplesInCategory<-rownames(metaData[metaData[,metCol]==categories[i],])
    # get the sub distance matrix
    subDissMatrix<-dissMatrixTable[samplesInCategory,samplesInCategory]
    # remove the upper triangle
    subDissMatrix[upper.tri(subDissMatrix)]<-NA
    subsamplesNames<-rownames(subDissMatrix)
    # remove the diagonal
    diag(subDissMatrix)<-NA
    # convert to matrix
    subDissMatrix<-as.matrix(subDissMatrix)
    # add row names
    rownames(subDissMatrix)<-subsamplesNames
    # add column names
    colnames(subDissMatrix)<-subsamplesNames
    subDissMatrix<-melt(subDissMatrix)
    # remove the NA
    subDissMatrix<-subDissMatrix[!is.na(subDissMatrix$value),]
    # add to the list
    distCatList[[i]]<-subDissMatrix
    # save the sub distance matrix
    #write.table(subDissMatrix, file=paste(ResultsFolder,"/",outfilename,".",categories[i],".dissMatrix.tsv",sep=""), sep="\t", quote=F, row.names=T, col.names=T)
  }
  # add names to the list
  names(distCatList)<-categories
  # convert to a data frame with the category as a column
  catDataFrame<-data.frame(sample1=character(), sample2=character(), value=numeric(), category=character(), stringsAsFactors=F)
  
  for (i in 1:length(distCatList)) {
    catDataFrame<-rbind(catDataFrame, cbind(distCatList[[i]], category=names(distCatList)[i]))
  }
  return(catDataFrame)

}


DistCatFrame<-as.data.frame(matrix(NA, nrow = 0, ncol = 5))
colnames(DistCatFrame)<-c("sample1", "sample2", "value", "category", "col")
for (i in 1:ncol(metaData)) {
print(colnames(metaData)[i])
 ColCatDataFrame<-getDistanceData(colnames(metaData)[i],metaData,dissMatrixTable)
 ColCatDataFrame$col<-colnames(metaData)[i]
 DistCatFrame<-rbind(DistCatFrame, ColCatDataFrame)
}



# make it a function
getDistanceBoxplot<-function(col,metaData,DistCatFrame)
{
 dplot<-ggplot(DistCatFrame[DistCatFrame$col==colnames(metaData)[col],], aes(x=category, y=value, fill=category)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=paste("Distance boxplot for ",colnames(metaData)[col],sep=""), x="Category", y="Distance")+
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
    # remove _ from the category
    scale_x_discrete(labels = function(x) gsub("_"," ",x))
  return(dplot)
}


# convert to function
getTtestMatrixPlot<-function(col,metaData,DistCatFrame)
{
  # create a matrix of comparisons
  catnume<-unique(metaData[,col])
  compMatrix<-matrix(NA, nrow = length(catnume), ncol = length(catnume))
  rownames(compMatrix)<-catnume
  colnames(compMatrix)<-catnume
  # calculate the t-test for each comparison
  for (i in 1:nrow(compMatrix)) {
    for (j in 1:ncol(compMatrix)) {
      if (i!=j) {
        catDataFrame<-DistCatFrame[DistCatFrame$col==colnames(metaData)[col],]
        catDataFrame<-catDataFrame[catDataFrame$category %in% c(rownames(compMatrix)[i],colnames(compMatrix)[j]),]
        ttest<-t.test(catDataFrame$value~catDataFrame$category)
        cat("Comparing ",rownames(compMatrix)[i]," and ",colnames(compMatrix)[j],"\n")
        cat("p-value = ",ttest$p.value,"\n")
        compMatrix[i,j]<-ttest$p.value
      }
    }
  }
  # replace p-values with stars * < 0.05, ** < 0.01, *** < 0.001
  compMatrixStars<-compMatrix
  # repalce > 0.05 with 0
  compMatrixStars[compMatrixStars>0.05]<-NA
  compMatrixStars[compMatrixStars<0.001]<-3
  compMatrixStars[compMatrixStars<0.01]<-2
  compMatrixStars[compMatrixStars<0.05]<-1
  # replace NA with 0
  compMatrixStars[is.na(compMatrixStars)]<-0
  # replace _ with " "
  colnames(compMatrixStars)<-gsub("_"," ",colnames(compMatrixStars))
  rownames(compMatrixStars)<-gsub("_"," ",rownames(compMatrixStars))
  # convert colnames to character
  colnames(compMatrixStars)<-as.character(colnames(compMatrixStars))
  rownames(compMatrixStars)<-as.character(rownames(compMatrixStars))

  # plot the matrix using ggplot and color with 4 colors
  matPlot<-ggplot(melt(compMatrixStars), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="T-test p-values matrix", x="Category", y="Category")+
    # make it categorical
    scale_fill_gradientn(colours=c("white", "blue", "green", "red"), limits=c(0,3), breaks=c(0,1,2,3), labels=c("Not significant", "< 0.05", "< 0.01", "< 0.001"))+
    # remove x and y titles
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  return(matPlot)
}



# create for all columns
for (coln in 1:ncol(metaData)) {
  dplot<-getDistanceBoxplot(coln,metaData,DistCatFrame)
  tplot<-getTtestMatrixPlot(coln,metaData,DistCatFrame)
  mplot<-grid.arrange(dplot, tplot, ncol=2, widths=c(3, 1))
  ggsave(paste(ResultsFolder,"/",outfilename,".",colnames(metaData)[coln],".distance.boxplot.png",sep=""), mplot, width=20, height=5, units="in", dpi=300)
}
