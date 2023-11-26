#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################



#!/usr/bin/env Rscript
rm(list=ls())
# library
library(LEA)
library(ggplot2)
library(parallel)
library(vcfR)


# use args
args<-commandArgs(trailingOnly = TRUE)


DataFile<-args[1]
genoDataname<-args[2]

resultFolder<-args[3]
outFolder<-args[4]

# print the arguments
print("DataFile")
print(DataFile)
print("genoDataname")
print(genoDataname)
print("resultFolder")
print(resultFolder)
print("outFolder")
print(outFolder)

# # get geno  
# geno<-vcfR2genind(vcf)
# geno<-rownames(geno@tab)


# path<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_2_13_9_2023"
# DataFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_2_13_9_2023/Data/genotype/Adult.vcf" # in 0,1,2 format (012) generated from vcftools
# resultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_2_13_9_2023/results"
# outFolder<-"populationStructureAdult"
# genoDataname<-"AdultGenotype"



# create folders results and outFolder if not exist
dir.create(file.path(resultFolder,outFolder),showWarnings = FALSE)

# read vcf file
vcf<-read.vcfR(DataFile)
genonames<-colnames(vcf@gt)
# remove the first column
genonames<-genonames[-1]

TaxaNames<-data.frame(Taxa=genonames)
# save the TaxaNames
write.table(TaxaNames,file=paste(resultFolder,"/",outFolder,"/","TaxaNames.txt",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


# convert vcf to geno
GenoFile<-vcf2geno(DataFile,paste0(resultFolder,"/",outFolder,"/",genoDataname,".geno"))
# save the IDs order



structure.iterations<-10000 #number of repeats
maxK<-12   #maximum number of Ks
#detect number of cores
ncores<-detectCores()   # number of cores to use for parallel processing


# Do the analysis using LEA , The snfm folder will be created in Data folder not in Results folder
project = snmf(GenoFile,
              K = 1:maxK,
              entropy = TRUE,iterations = structure.iterations,
              repetitions = 10,CPU=ncores-3,
              project = "new")


#  project = load.snmfProject("results/populationStructureAdult/AdultGenotype.snmfProject")


# plot cross-entropy criterion of all runs of the project
pdf(paste(resultFolder,"/",outFolder,"/","cross-entropy.pdf",sep=""))
plot(project, cex = 1.2, col = "lightblue", pch = 19)
dev.off()

# create folder for best ks
dir.create(file.path(resultFolder,"/",outFolder,"bestKs"),showWarnings = FALSE)

# load project using LEA

for (kn in 2:maxK) {
    bestK = which.min(cross.entropy(project, K = kn))
    print(bestK)
    # move the best K file 
    snmfFolder<-paste(resultFolder,"/",outFolder,"/",genoDataname,".snmf",sep="")
    # # read the best K file
    StrutcMat<-read.table(paste(snmfFolder,"/K",kn,"/run",
                                bestK,"/",genoDataname,"_r",bestK,".",kn,".Q",sep=""))
    # # save the best K file to bestKs folder with run number as name
    write.table(StrutcMat,file=paste(resultFolder,"/",outFolder,"/bestKs/","K",kn,".Q",sep=""),
    quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}