#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################




#rm(list=ls())
# library
library(LEA)
library(ggplot2)
library(parallel)
library(vcfR)

# use args
args<-commandArgs(trailingOnly = TRUE)


DataFile<-args[1]
resultFolder<-args[2]
maxK<-args[3]
maxiter<-args[4]


# DataFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/Data/Field/Field_filtered.vcf"
# resultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/test"

# create folders results and outFolder if not exist
dir.create(file.path(resultFolder),showWarnings = FALSE)

# read vcf file
vcf<-read.vcfR(DataFile)
genonames<-colnames(vcf@gt)
# remove the first column
genonames<-genonames[-1] # nolint

TaxaNames<-data.frame(Taxa=genonames)
# save the TaxaNames
write.table(TaxaNames,file=paste(resultFolder,"/","TaxaNames.txt",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

# convert vcf to geno
GenoFile<-vcf2geno(DataFile,paste0(resultFolder,"/Genotypes.geno"))
# save the IDs order

# if maxK is not provided then use 12
if (maxK=="") {
    maxK<-12
}
# if maxiter is not provided then use 10000
if (maxiter=="") {
    maxiter<-10000
}

structure.iterations<-as.numeric(maxiter) # number of iterations
maxK<-as.numeric(maxK)          # maximum number of clusters
#detect number of cores
ncores<-detectCores()   # number of cores to use for parallel processing

# Do the analysis using LEA , The snfm folder will be created in Data folder not in Results folder
project = snmf(GenoFile,
              K = 1:maxK,
              entropy = TRUE,iterations = structure.iterations,
              repetitions = 10,CPU=ncores-3,
              project = "new")

# save the project
# plot cross-entropy criterion of all runs of the project
pdf(paste(resultFolder,"/","cross-entropy.pdf",sep=""))
plot(project, cex = 1.2, col = "lightblue", pch = 19)
dev.off()

# create folder for best ks
dir.create(file.path(resultFolder,"/bestKs"),showWarnings = FALSE)

# load project using LEA
for (kn in 2:maxK) {
    bestK = which.min(cross.entropy(project, K = kn))
    # move the best K file 
    snmfFolder<-paste(resultFolder,"/Genotypes.snmf",sep="")
    # # read the best K file
    StrutcMat<-read.table(paste(snmfFolder,"/K",kn,"/run",
                                bestK,"/Genotypes","_r",bestK,".",kn,".Q",sep=""))
    # # save the best K file to bestKs folder with run number as name
    write.table(StrutcMat,file=paste(resultFolder,"/bestKs/","K",kn,".Q",sep=""),
    quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}