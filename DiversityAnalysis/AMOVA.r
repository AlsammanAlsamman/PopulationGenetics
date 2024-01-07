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
library(ggplot2)
# brewer
library(RColorBrewer)
library(ape)
library(pegas)
# reading the input file

vcffile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_sample.vcf"
metadataFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/sample_meta.tsv"
popName<-"Pedigree"
subPopName<-"P1.IG"
ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/amova_results/"

# reading the vcf file
vcf<-read.vcf(vcffile)
vcfR<-read.vcfR(vcffile)
# reading the metadata file
metadata<-read.table(metadataFile,header = T,sep = filesep,stringsAsFactors = F, row.names = 1)
# select the required columns
metadata<-metadata[,c(popName,subPopName)]

# # convert the genclone object
# vcfR<-read.vcfR(vcffile)
# # select the required samples
# vcfR<-vcfR[rownames(metadata),]
# rownames(metadata)

colnames(vcf)<-paste(getCHROM(vcfR),getPOS(vcfR),sep="_")

# check if there is any duplicated column names
if (anyDuplicated(colnames(vcf))) {
  print("There is duplicated marker locations")
  print("Removing duplicated marker locations")
  # remove duplicated column names
  vcf<-vcf[,!duplicated(colnames(vcf))]
}


# select from the vcf file the required samples
vcf<-vcf[rownames(metadata),]
# select the required samples
metadata<-metadata[rownames(vcf),]
# convert to genlight object
genIndObj<-loci2genind(vcf,ploidy=1)
# select from the genIndObj the required samples
metadata<-metadata[rownames(genIndObj@tab),]
# check sevearl things
# is the metadata and the genIndObj have the same number of samples
if (nrow(metadata)!=nInd(genIndObj)) {
  print("The number of samples in the metadata and the genIndObj are not the same")
  print("Please check the metadata file")
  quit()
}


# is the metadata and the genIndObj have the same samples
if (!all(rownames(metadata)==rownames(genIndObj))) {
  print("The samples in the metadata and the genIndObj are not the same")
  print("Please check the metadata file")
  quit()
}

















# order the metadata according to the genIndObj
metadata<-metadata[rownames(genIndObj@tab),]
pop<-metadata[,popName]
subPop<-metadata[,subPopName]
# create strata
strata.df<-data.frame(pop=pop,subPop=subPop)

strata(genIndObj)<-strata.df
setPop(genIndObj)<-~pop/subPop
genCloneObj<-as.genclone(genIndObj)
# change polyploidy to 1
ploidy(genCloneObj)<-1

table(strata(genIndObj, ~pop/subPop, combine = FALSE))  # Subpopulations
genIndObj
pop(genIndObj)
pop
subPop
# usign Using quasieuclid correction method
genIndObjAmova <- poppr.amova(genIndObj, ~pop/subPop)

genIndObjAmova
set.seed(1999)
genIndObjSignif   <- randtest(genIndObjAmova, nrepet = 999)
genIndObjSignif
plot(genIndObjSignif)


# using the sample data


# reading the vcf file
vcf<-read.vcf(vcffile)
# tak only thf first 1000 markers
vcf<-vcf[1:1000,]
# convert the genclone object
genIndObj<-loci2genind(vcf,ploidy=1)
# create random population
pop<-sample(1:3,nInd(genIndObj),replace = T)
# create random subpopulation
subPop<-sample(1:3,nInd(genIndObj),replace = T)
# add the pop and subpop to the strata
strata.df<-data.frame(pop=pop,subPop=subPop)
strata(genIndObj)<-strata.df
setPop(genIndObj)<-~pop/subPop

# do the amova
genIndObjAmova <- poppr.amova(genIndObj, ~pop/subPop)
# do the randome test
genIndObjSignif   <- randtest(genIndObjAmova, nrepet = 999)
genIndObjSignif
# plot the results
plot(genIndObjSignif)





#######################################################################################
# genIndObjAmovaCC <- poppr.amova(genIndObjAmova, ~pop/subPop, clonecorrect = TRUE)
# # ploidy(genlight) 
# x.dist <- dist(genlight)
# tree <- aboot(genlight, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
# cols <- brewer.pal(n = nPop(genlight), name = "Dark2")
# plot.phylo(tree, cex = 0.8, font = 2, adj = 0) # , tip.color =  cols[pop(gl.rubi)]
# nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
# #legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
# legend('topleft', legend = c("CA","OR","WA"), fill = cols, border = FALSE, bty = "n", cex = 2)
# axis(side = 1)
# title(xlab = "Genetic distance (proportion of loci that are different)")
# # converting the vcf file to genlight object
# genlight<-vcfR2genlight(vcf)
# ## This is a genclone object