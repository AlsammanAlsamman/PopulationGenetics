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

  # VCFFILE<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Morocco_meeting_final/Data/BreadWheatBefore_filtered_selectedSamples.vcf.recode.vcf"
  # resultsFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Morocco_meeting_final/Results/PCA"
  # dataFolder<-resultsFolder
  # MetaDataFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Morocco_meeting_final/Data/passportData.order.representative.tsv"

# Check folders 
# if folder path does not end with path separator then add it
if (substr(resultsFolder,nchar(resultsFolder),nchar(resultsFolder))!=pathSeparator) {
  resultsFolder<-paste(resultsFolder,pathSeparator,sep="")
}
if (substr(dataFolder,nchar(dataFolder),nchar(dataFolder))!=pathSeparator) {
  dataFolder<-paste(dataFolder,pathSeparator,sep="")
}

library(SNPRelate)

VCFFILE<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/BreadWheatBefore_Filtered.vcf"
snpgdsVCF2GDS(VCFFILE, paste("BreadGenotype",".gds",sep=""), method="biallelic.only", verbose=TRUE)

dataFolder <- "/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data"
gdsPath <- paste(dataFolder, "BreadGenotype.gds", sep = "")
# # open the gds file
# seqVCF2GDS(VCFFILE, verbose = TRUE, out.fn = gdsPath)

gdsFile <- seqOpen(gdsPath, readonly = FALSE)

# close the gds file
closefn.gds(genofile)
closefn.gds(gdsFile)
# for %
library(SNPRelate)
library(dplyr)
library(SeqArray)
library(SNPRelate)
library(pander)
library(scales)
library(magrittr)
library(tidyverse)

snpgdsLDpruning(gdsFile, sample.id = NULL, autosome.only = TRUE, maf = 0.05, missing.rate = 0.25, method = "corr", ld.threshold = 0.4, verbose = TRUE)


seqResetFilter(gdsFile)
gdsFile




