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

if (length(args) < 3) {
  print(paste("Usage: Rscript plot_GWAS_GEMMA.r <gwasfile> <traitName> <ResultFolder> <fdrcutoff> <qthreshold> <chrTable>",sep=""))
  quit()
}

# snpmetaFile<- args[1]
# metFile<-args[2]

snpmetaFile<- "/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/moremoremoreResults/SNP_vs_Meta/snpMeta.txt"
metaFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/moremoremoreResults/SNP_vs_Meta/passportData.order.representative.tsv"

# read the snp meta file
snpMeta<-read.table(snpmetaFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)
Meta<-read.table(metaFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1)

# create a matrix with m x snp columns

snpmmat<-matrix(nrow = ncol(snpMeta), ncol = ncol(Meta))
# add the column names and row names
rownames(snpmmat)<-colnames(snpMeta)
colnames(snpmmat)<-colnames(Meta)



# calculate the chi square test for each snp with each meta data
for (i in 1:ncol(snpMeta)) {
  for (j in 1:ncol(Meta)) {
  mcol<-colnames(Meta)[j]
  snpcol<-colnames(snpMeta)[i]
  mvalues<-unique(Meta[,mcol])
  snpvalues<-unique(snpMeta[,snpcol])
  snpvalues
  # create a matrix with the number of rows equal to the number of unique values in the meta data
  smmat<-matrix(nrow = length(mvalues), ncol = length(snpvalues))
  for (m in 1:length(mvalues)) {
    for (s in 1:length(snpvalues)) {
      # get the number of rows where the meta data is equal to the current value and the snp data is equal to the current value
      smmat[m,s]<-nrow(snpMeta[which(Meta[,mcol]==mvalues[m] & snpMeta[,snpcol]==snpvalues[s]),])
      # add the row and column names
      rownames(smmat)<-mvalues
      colnames(smmat)<-snpvalues
    }
  }
  pvalue<-NA
  # perform the fisher test
  tryCatch({

  pvalue<-  fisher.test(smmat)$p.value
  }, error=function(e) {
    pvalue<-NA
  })
# add the pvalue to the matrix
  snpmmat[i,j]<-pvalue

  }
}

# remove all rows and columns with NA values
snpmmat<-snpmmat[!apply(snpmmat, 1, function(x) all(is.na(x))),]
snpmmat<-snpmmat[,!apply(snpmmat, 2, function(x) all(is.na(x)))]

# melt the matrix
library(reshape2)
snpmmat<-melt(snpmmat)
snpmmat
# select the rows with pvalue < 0.05
snpmmat<-snpmmat[which(snpmmat$value<0.1),]
snpmmat


# add the row names as a column "Taxa"
Meta$Taxa<-rownames(Meta)
snpMeta$Taxa<-rownames(snpMeta)

# select first meta col
mCol<-1
snpCol<-1

mColNames<-colnames(Meta)
snpColNames<-colnames(snpMeta)
# get the meta data
metaData<-Meta[,c("Taxa","K5")]
colnames(Meta)
snpData<-snpMeta[,c("Taxa","TaDArTAG007312")]
snpData






#merge the two data frames
mergedData<-merge(metaData, snpData, by="Taxa", all=TRUE)
# barplot with stacked bars showing the allele frequency in each group
library(ggplot2)
pdf("/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/moremoremoreResults/SNP_vs_Meta/snpMeta.pdf")
ggplot(mergedData, aes(x=mergedData[,2], fill=mergedData[,3])) + 
  geom_bar(position="stack") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title="SNP vs Meta", x=colnames(mergedData)[2], y="Count")+
  # legend title
  guides(fill=guide_legend(title="Allele State"))
dev.off()

# is there an association between the allele state and the meta data
# chi square test
chisq.test(mergedData[,2], mergedData[,3])
