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
library(openxlsx)





VCFFILE<-args[1]
resultsFolder<-args[2]

if (length(args)==3) {
  print("Metadata file is provided")
  metaFile<-args[3]
}



snpgdsVCF2GDS(VCFFILE, paste("BreadGenotype",".gds",sep=""), method="biallelic.only", verbose=TRUE)

genofile <- openfn.gds(paste("BreadGenotype",".gds",sep=""))

# calculate the kinship matrix by IBS
# kinshipMatrix<-snpgdsIBDKING(genofile, num.thread=4)
kinshipMatrixObject<-snpgdsIBS(genofile, num.thread=4, autosome.only = F)
# kinshipMatrix
# kinshipMatrix$kinship

kinshipMatrix<-kinshipMatrixObject$ibs
rownames(kinshipMatrix)<-kinshipMatrixObject$sample.id
colnames(kinshipMatrix)<-kinshipMatrixObject$sample.id
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

# add metadata to the kinship matrix

print(metaFile)
# if metadata file is not provided then stop
if (length(args)<3) {
  stop("Metadata file is not provided", call.=FALSE)
}

metadata<-read.table(metaFile,header=T,sep="\t",stringsAsFactors = F, row.names = 1)



pdf(paste(resultsFolder,"/kinshipMatrixWithMetadata.pdf",sep=""),width=15,height=15)
pheatmap(kinshipMatrix,main="Kinship Matrix",annotation_col=metadata, annotation_row = metadata, 
# remove row and column names
show_colnames = F, show_rownames = F)
dev.off()

png(paste(resultsFolder,"/kinshipMatrixWithMetadata.png",sep=""),width=1500,height=1500)
pheatmap(kinshipMatrix,main="Kinship Matrix",annotation_col=metadata, annotation_row = metadata,
# remove row and column names
show_colnames = F, show_rownames = F)
dev.off()

# close the gds file
closefn.gds(genofile)

# plot heatmap using  heatmap.2 function
png(paste(resultsFolder,"/kinshipMatrix.png",sep=""),width=1000,height=1000)
heatmap(kinshipMatrix,main="Kinship Matrix",trace="none",col=brewer.pal(9,"Blues"))
dev.off()

# save the kinship matrix as a text file
write.table(as.data.frame(kinshipMatrix),file=paste(resultsFolder,"/kinshipMatrix.txt",sep=""),sep="\t")

# save the kinship matrix as a xlsx file
write.xlsx(as.data.frame(kinshipMatrix),paste(resultsFolder,"/kinshipMatrix.xlsx",sep=""))

# 