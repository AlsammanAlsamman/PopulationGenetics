#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
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

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} #else if (length(args)==1) {
  # default output file
  #args[2] = "out.txt"
#}

# loading libraries
library(SNPRelate)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)




VCFFILE<-args[1]
resultsFolder<-args[2]
dataFolder<-args[3]
MetaDataFile<-args[4]

# if data folder is not provided then use the results folder
if (length(args)==2) {
  dataFolder<-resultsFolder
}

# VCFFILE<-"/home/samman/Documents/publishing/Drought_Ahmed_Master/NewAnalysis/Data/DR.vcf"
# resultsFolder<-"/home/samman/Documents/publishing/Drought_Ahmed_Master/NewAnalysis/PCA"
# dataFolder<-resultsFolder
# MetaDataFile="/home/samman/Documents/publishing/Drought_Ahmed_Master/NewAnalysis/Data/PassPort.txt"



# Check folders 
# if folder path does not end with path separator then add it
if (substr(resultsFolder,nchar(resultsFolder),nchar(resultsFolder))!=pathSeparator) {
  resultsFolder<-paste(resultsFolder,pathSeparator,sep="")
}
if (substr(dataFolder,nchar(dataFolder),nchar(dataFolder))!=pathSeparator) {
  dataFolder<-paste(dataFolder,pathSeparator,sep="")
}

snpgdsVCF2GDS(VCFFILE, paste("BreadGenotype",".gds",sep=""), method="biallelic.only", verbose=TRUE)
genofile <- openfn.gds(paste("BreadGenotype",".gds",sep=""))
PCA <- snpgdsPCA(genofile,num.thread=4, autosome.only=FALSE) # The autosome.only=FALSE option is used to include the X chromosome in the analysis
# close the gds file
closefn.gds(genofile)
PCA.eigenvect<-as.data.frame(PCA$eigenvect)
colnames(PCA.eigenvect)<-paste0("PC",1:ncol(PCA.eigenvect))
PCA.eigenval<-PCA$eigenval
# remove NAs
PCA.eigenval<-PCA.eigenval[!is.na(PCA.eigenval)]
# calculate the percentage of variance explained by each PC
PCA.eigenval.percentage<-PCA.eigenval/sum(PCA.eigenval)*100
# plot the eigenvalues
eigenvalues<-data.frame(PC=1:length(PCA.eigenval),Eigenvalue=PCA.eigenval.percentage)

#create folder for SNPrelate PCA results
dir.create(paste(resultsFolder,"/SNPrelate_PCA",sep=""))


p<-ggplot(eigenvalues,aes(x=PC,y=Eigenvalue))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="PC",y="Eigenvalue (%)")

ggsave(paste(resultsFolder,"/SNPrelate_PCA/PCA_eigenvalues.png",sep=""),p,width=10,height=5)

PC1.variance<-PCA.eigenval.percentage[1]
PC2.variance<-PCA.eigenval.percentage[2]

# plot the PCA with out metadata
p<-ggplot(PCA.eigenvect,aes(x=PC1,y=PC2))+geom_point()+theme_bw()+labs(x=paste0("PC1 (",round(PC1.variance,2),"%)"),y=paste0("PC2 (",round(PC2.variance,2),"%)"))

ggsave(paste(resultsFolder,"/SNPrelate_PCA/PCA.png",sep=""),p,width=10,height=5)

# if metadata file exists then plot the PCA with metadata

# check if metadata file exists
MetaData<-read.csv(MetaDataFile,header=TRUE,sep="\t",stringsAsFactors = FALSE, row.names = 1)
# add sample.id column to metadata
rownames(PCA.eigenvect)<-PCA$sample.id
# merge metadata with eigenvectors
PCA.eigenvect.meta<-merge(PCA.eigenvect,MetaData,by="row.names")

# foreach column in metadata, color the PCA plot
for (i in 1:ncol(MetaData)){
  # convert to factor
  metCol<-colnames(MetaData)[i]
  metCol<-as.factor(PCA.eigenvect.meta[,metCol])
  
  p<-ggplot(PCA.eigenvect.meta,aes(x=PC1,y=PC2))+
  geom_point(aes(color=metCol))+
  theme_bw()+labs(x=paste0("PC1 (",round(PC1.variance,2),"%)"),y=paste0("PC2 (",round(PC2.variance,2),"%)"))+
  # legend at bottom
  theme(legend.position="left")
  ggsave(paste(resultsFolder,"/SNPrelate_PCA/PCA_",colnames(MetaData)[i],".png",sep=""),p,width=5,height=5)
  # one plot with legend at all
  p<- p + theme(legend.position="none")
  ggsave(paste(resultsFolder,"/SNPrelate_PCA/PCA_",colnames(MetaData)[i],"_no_legend.png",sep=""),p,width=5,height=5)
}