#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  print(paste("Usage: Rscript markerDenisty.r vcfFile ResultFolder windowSize maxMarkerCount",sep=" "))
  quit()
}
#######################################################################################
# What is this system linux or windows
system<-Sys.info()[1]
# if system is windows then change the path separator
if (system=="Windows") {
  pathSeparator<-"\\"
} else {
  pathSeparator<-"/"
}
#######################################################################################
# from system get the value of the environment variabl of export AlsammanPopulationGeneticsPath
functionsPath<-Sys.getenv("AlsammanPopulationGeneticsPath")
functionsPath<-paste(functionsPath,"SNP_QC/markerDenisty/functions",sep="")
# read all r source scripts in function folder
sourcesList<-list.files(paste(functionsPath,sep=""),pattern="*.r",full.names=TRUE)
# if the list is empty then exit
if (length(sourcesList)==0) {
  print(paste("No functions found in you should use the configureMe.sh script to set the path to the functions folder",sep=""))
  quit()
}

# source all functions
for (i in 1:length(sourcesList)) {
  source(sourcesList[i])
}
#######################################################################################
# load libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
#######################################################################################

# plot marker denisty across genome
vcfFile<-args[1]
ResultFolder<-args[2]
windowSize<-as.numeric(args[3])
maxMarkerCount<-as.numeric(args[4])

# vcfFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData.vcf"
# ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS"

# if the window size is not provided then use 50Mb
if (is.null(windowSize)) {
  windowSize<-1000000
  print(paste("Window size is not provided, using default window size of 1Mb",sep=""))
}
if (is.null(maxMarkerCount)) {
  maxMarkerCount<-50
  print(paste("Max marker count is not provided, using default max marker count of 50",sep=""))
}

snps<-read.table(vcfFile,sep="\t",header=F,blank.lines.skip=TRUE,
                 comment.char = "#")
# get only the first 2 columns
snps<-snps[,1:2]
# rename columns
colnames(snps)<-c("chrom","pos")
# rename NA chromosomes to Unknown
snps$chrom[is.na(snps$chrom)]<-"Unknown"
# calculate marker density
MarkerDen<-MarkerDenisty(snps, AcrossGenome=TRUE, ForeachChromosome=TRUE, window.size=windowSize)

#
DenistyPlot<-MarkerDenistyByWindowPlot(MarkerDen$chromDensity, maxMarkerCount)
DenistyPlot<-DenistyPlot 

# save plot
ggsave(paste(ResultFolder,"/MarkerDensityByWindow.png",sep=""),DenistyPlot, width = 10, height = 10, units = "in", dpi = 300)
# pdf
ggsave(paste(ResultFolder,"/MarkerDensityByWindow.pdf",sep=""),DenistyPlot, width = 10, height = 10, units = "in", dpi = 300)



# Future work
# # barplot using ggplot2
# p2<-ggplot(MarkerDen$eachChromTotal, aes(x = chrom, y = count, fill = chrom)) +
#   geom_bar(stat = "identity")
# # flip
# p2<-p2 + coord_flip()
# p2<-p2 + # remove x-axis label and ticks
#   theme(axis.title.x=element_blank())
# ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "left")
# p<-MarkerDenistyPlot(MarkerDen, 15)
