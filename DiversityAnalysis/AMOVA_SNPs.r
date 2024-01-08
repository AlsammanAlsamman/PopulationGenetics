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
}

if (length(args) < 5) {
  print(paste("Usage: Rscript AMOVA.r vcffile metadata popName subPopName ResultFolder"))
  quit()
}

# vcffile<-args[1]
# metadata<-args[2]
# popName<-args[3]
# subPopName<-args[4]
# ResultFolder<-args[5]

#######################################################################################
filesep<-"\t"
# if the file name contains the word "csv" then the file separator is comma
if (grepl("csv",metadata)) {
  filesep<-","
}

#######################################################################################

# Load the required libraries
library(poppr)
library(vegan)
library(adegenet)
library(vcfR)
library(ape)
library(pegas)
# reading the input file
library(dplyr)
library(tidyr)
library(reshape2)

vcffile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_sample.vcf"
SNPAlleleListFile = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/TargetSNPs_test.txt"
ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/snp_amova_results"

# if the result folder does not exist then create it
if (!dir.exists(ResultFolder)) {
  dir.create(ResultFolder)
}


# get a matrix of the genotypes
convertGeno2Binary<-function(geno) {
  # convert the genotypes to binary
  geno[geno=="0/0"]<-0
  geno[geno=="0/1"]<-1
  geno[geno=="1/0"]<-1
  geno[geno=="1/1"]<-2
  geno[geno=="./."] <- NA
  return(geno)
}

# reading the SNP list
SNPList<-read.table(SNPAlleleListFile,header=FALSE,sep="\t")
# reading the vcf file
vcfR<-read.vcfR(vcffile)
# get the genotypes
geno<-vcfR@gt
# add the SNP names
rownames(geno)<-as.character(vcfR@fix[,3])
geno<-geno[rownames(geno) %in% SNPList[,1],]
# convert the genotypes to binary
geno<-convertGeno2Binary(geno)
# group samples by marker value
geno<-t(geno)
# change the values to characters , 0->A, 1->H, 2->B
geno[geno==0]<-"A"
geno[geno==1]<-"H"
geno[geno==2]<-"B"
geno[is.na(geno)]<-"N"
# convert to data frame
# save the snpMetadata to file in the run folder
snpMetaFile<-paste(ResultFolder,pathSeparator,"snpMeta.txt",sep="")
# remove the first row in the geno
geno<-geno[-1,]
# write the snp metadata to file
write.table(geno,snpMetaFile,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
# get the target SNPs
targetSNPs<-colnames(geno)
# Read The VCF file and convert it to genind object
vcfData <- read.vcf(vcffile)
# convert the vcf object to genind object
g<-loci2genind(vcfData)
# convert the genind object to genclone object
g<-as.genclone(g)
# set the ploidy to 2
ploidy(g)<-2
# filter the genind object for missing data and minor allele frequency (MAF)
g<-missingno(g, cutoff = 0.1)
# get samples from the metadata file
vcfSamples<-rownames(g@tab)
# order the samples in the metadata file as in the genind object
geno<-geno[vcfSamples,]

# convert geno to data frame
geno<-as.data.frame(geno)
strata(g)<-geno

amovacc<-poppr.amova(g,~TaDArTAG000007, method="ade4")

# Results
amovaResults<-data.frame(Df=0,SumSq=0,MeanSq=0,SNP=NA)
# components
amovaComponents<-data.frame(Sigma=0,Percentage=0,SNP=NA)
# statphi
amovaStatphi<-data.frame(Phi=0,SNP=NA)

# do the same for each SNP
for (i in 1:length(targetSNPs)) {
  # create a folder for the current SNP
  currentSNP<-targetSNPs[i]
  snpFormular<-paste("~",currentSNP,sep="")
  snpFormular<-as.formula(snpFormular)
  
  # run the amova
  snpAmovacc<-poppr.amova(g,snpFormular, method="ade4")

  # add the results to the data frame
  # results
  SNP.amovaResults<-snpAmovacc$results
  SNP.amovaResults$SNP<-currentSNP
  colnames(SNP.amovaResults)<-c("Df","SumSq","MeanSq","SNP")
  amovaResults<-rbind(amovaResults,SNP.amovaResults)
  # components
  SNP.amovaComponents<-snpAmovacc$componentsofcovariance
  SNP.amovaComponents$SNP<-currentSNP
  colnames(SNP.amovaComponents)<-c("Sigma","Percentage","SNP")
  amovaComponents<-rbind(amovaComponents,SNP.amovaComponents)
  # statphi
  SNP.amovaStatphi<-snpAmovacc$statphi
  SNP.amovaStatphi$SNP<-currentSNP
  colnames(SNP.amovaStatphi)<-c("Phi","SNP")
  amovaStatphi<-rbind(amovaStatphi,SNP.amovaStatphi)
}

# remove the first row
amovaResults<-amovaResults[-1,]
amovaComponents<-amovaComponents[-1,]
amovaStatphi<-amovaStatphi[-1,]

# write the results to file
amovaResultsFile<-paste(ResultFolder,pathSeparator,"amovaResults.csv",sep="")
write.table(amovaResults,amovaResultsFile,sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

amovaComponentsFile<-paste(ResultFolder,pathSeparator,"amovaComponents.csv",sep="")
write.table(amovaComponents,amovaComponentsFile,sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

amovaStatphiFile<-paste(ResultFolder,pathSeparator,"amovaStatphi.csv",sep="")
write.table(amovaStatphi,amovaStatphiFile,sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

# set.seed(1999)
#monpopsignif   <- randtest(amovacc, nrepet = 999)