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
library(pheatmap)




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


dataFolder


snpgdsVCF2GDS(VCFFILE, paste("BreadGenotype",".gds",sep=""), method="biallelic.only", verbose=TRUE)

genofile <- openfn.gds(paste("BreadGenotype",".gds",sep=""))

# calculate the kinship matrix by IBS
# kinshipMatrix<-snpgdsIBDKING(genofile, num.thread=4)
kinshipMatrix<-snpgdsIBS(genofile, num.thread=4, autosome.only = F)
# kinshipMatrix
# kinshipMatrix$kinship

kinshipMatrix<-kinshipMatrix$ibs
rownames(kinshipMatrix)<-genofile$sample.id
colnames(kinshipMatrix)<-genofile$sample.id
# save the kinship matrix as a text file
# write.table(kinshipMatrix,file=paste(resultsFolder,"/kinshipMatrix.txt",sep=""),sep="\t",quote=F,col.names=NA,row.names=NA)

# plot the kinship matrix
png(paste(resultsFolder,"/kinshipMatrix.png",sep=""),width=1000,height=1000)
pheatmap(kinshipMatrix,main="Kinship Matrix")
dev.off()

pdf(paste(resultsFolder,"/kinshipMatrix.pdf",sep=""),width=10,height=10)
pheatmap(kinshipMatrix,main="Kinship Matrix")
dev.off()



# save session info to a text file
sink(paste(resultsFolder,"/sessionInfo.txt",sep=""))
sessionInfo()
sink()
