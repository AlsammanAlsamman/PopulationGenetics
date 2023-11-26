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

# if data folder is not provided then use the results folder
if (length(args)==2) {
  dataFolder<-resultsFolder
}

# VCFFILE<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/Adult/Adult.vcf"
# resultsFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/RESULTS/Adults"
# dataFolder<-resultsFolder

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

PCA <- snpgdsPCA(genofile,num.thread=4)
PCA.eigenvect<-as.data.frame(PCA$eigenvect)
colnames(PCA.eigenvect)<-paste0("PC",1:ncol(PCA.eigenvect))

PCA.eigenval<-PCA$eigenval
# remove NAs
PCA.eigenval<-PCA.eigenval[!is.na(PCA.eigenval)]

# calculate the percentage of variance explained by each PC
PCA.eigenval.percentage<-PCA.eigenval/sum(PCA.eigenval)*100


# plot

# plot the eigenvalues
eigenvalues<-data.frame(PC=1:length(PCA.eigenval),Eigenvalue=PCA.eigenval.percentage)

#create folder for SNPrelate PCA results
dir.create(paste(resultsFolder,"/SNPrelate_PCA",sep=""))


p<-ggplot(eigenvalues,aes(x=PC,y=Eigenvalue))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="PC",y="Eigenvalue (%)")
p
ggsave(paste(resultsFolder,"/SNPrelate_PCA/PCA_eigenvalues.png",sep=""),p,width=10,height=5)

PC1.variance<-PCA.eigenval.percentage[1]
PC2.variance<-PCA.eigenval.percentage[2]

# plot the PCA with out metadata
p<-ggplot(PCA.eigenvect,aes(x=PC1,y=PC2))+geom_point()+theme_bw()+labs(x=paste0("PC1 (",round(PC1.variance,2),"%)"),y=paste0("PC2 (",round(PC2.variance,2),"%)"))
p
ggsave(paste(resultsFolder,"/SNPrelate_PCA/PCA.png",sep=""),p,width=10,height=5)
