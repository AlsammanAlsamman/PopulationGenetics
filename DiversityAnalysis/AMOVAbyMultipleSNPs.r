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

if (length(args) < 3) {
  print(paste("Usage: Rscript AMOVAbySNPs.r vcffile SNPAlleleList ResultFolder"))
  quit()
}

vcffile<-args[1]
SNPAlleleListFile<-args[2]
ResultFolder<-args[3]

#######################################################################################

# Load the required libraries
library(poppr)
library(vegan)
library(adegenet)
library(vcfR)
library(ape)
library(ggplot2)
library(pegas)
library(dplyr)
library(tidyr)
library(reshape2)

vcffile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/someSNPs.vcf"
SNPAlleleListFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/TargetSNPs.txt"
ResultFolder="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Results/snpAMOVA"

# create a subfolder for the current run using the current date and time
dateTimestamp<-Sys.time()

# remove the spaces and the colons from the date and time
dateTimestamp<-gsub(" ","_",dateTimestamp)
dateTimestamp<-gsub(":","_",dateTimestamp)
ResultFolder<-paste(ResultFolder,pathSeparator,dateTimestamp,sep="")

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
# check if there are data for the SNPs in the list
if (nrow(geno)==0) {
  print("No data for the SNPs in the list")
  quit()
}
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
g<-missingno(g, cutoff = 0.05)
# get samples from the metadata file
vcfSamples<-rownames(g@tab)
# order the samples in the metadata file as in the genind object
geno<-geno[vcfSamples,]

# convert geno to data frame
geno<-as.data.frame(geno)
strata(g)<-geno

# Results
amovaResults<-data.frame(Df=0,SumSq=0,MeanSq=0,SNP=NA,rowname=NA)
# components
amovaComponents<-data.frame(Sigma=0,Percentage=0,SNP=NA,rowname=NA,varianceCategory=NA)
# statphi
amovaStatphi<-data.frame(Phi=0,SNP=NA,rowname=NA)

print("Starting AMOVA")
# do the same for each SNP
for (i in 1:length(targetSNPs)) {
  print(paste("Running AMOVA for SNP ",targetSNPs[i],sep=""))
  # create a folder for the current SNP
  currentSNP<-targetSNPs[i]
  snpFormular<-paste("~",currentSNP,sep="")
  snpFormular<-as.formula(snpFormular)
  
  # run the amova
  snpAmovacc<-poppr.amova(g,snpFormular, method="ade4", threads = 5)
  # add the results to the data frame
  # results
  SNP.amovaResults <- snpAmovacc$results
  SNP.amovaResults$SNP<-currentSNP
  SNP.amovaResults$rowname<-rownames(SNP.amovaResults)
  colnames(SNP.amovaResults)<-c("Df","SumSq","MeanSq","SNP","rowname")
  amovaResults<-rbind(amovaResults,SNP.amovaResults)
  
  # components
  SNP.amovaComponents<-snpAmovacc$componentsofcovariance
  SNP.amovaComponents$SNP<-currentSNP
  SNP.amovaComponents$rowname<-rownames(SNP.amovaComponents)
  SNP.amovaComponents$varianceCategory<-c("Between Populations","Between samples within populations","Within samples","Total")
  colnames(SNP.amovaComponents)<-c("Sigma","Percentage","SNP","rowname","varianceCategory")
  amovaComponents<-rbind(amovaComponents,SNP.amovaComponents)
  
  # statphi
  SNP.amovaStatphi<-snpAmovacc$statphi
  SNP.amovaStatphi$SNP<-currentSNP
  SNP.amovaStatphi$rowname<-rownames(SNP.amovaStatphi)
  colnames(SNP.amovaStatphi)<-c("Phi","SNP","rowname")
  amovaStatphi<-rbind(amovaStatphi,SNP.amovaStatphi)

  # sink the results to file
  amovasnpFile<-paste(ResultFolder,pathSeparator,"amova_",currentSNP,".txt",sep="")
  sink(amovasnpFile)
  print(snpAmovacc)
  sink()
}

# remove the first row
amovaResults<-amovaResults[-1,]
amovaComponents<-amovaComponents[-1,]
amovaStatphi<-amovaStatphi[-1,]

print("AMOVA is done")
print("Writing the results to file")
# write the results to file
amovaResultsFile<-paste(ResultFolder,pathSeparator,"amovaResults.csv",sep="")
write.table(amovaResults,amovaResultsFile,sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

amovaComponentsFile<-paste(ResultFolder,pathSeparator,"amovaComponents.csv",sep="")
write.table(amovaComponents,amovaComponentsFile,sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

amovaStatphiFile<-paste(ResultFolder,pathSeparator,"amovaStatphi.csv",sep="")
write.table(amovaStatphi,amovaStatphiFile,sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

# plot the variance components
# exclude the total
amovaComponents2plot<-amovaComponents[amovaComponents$varianceCategory!="Total",]

# for each SNP
snpsNames<-unique(amovaComponents2plot$SNP)
for (i in 1:length(snpsNames)) {
  currentSNP<-snpsNames[i]
  print(paste("Plotting variance components for SNP ",currentSNP,sep=""))
  # get the data for the current SNP
  currentSNPData<-amovaComponents2plot[amovaComponents2plot$SNP==currentSNP,]
  # abbreviate the variance category names
  currentSNPData$AbbrevCats<-abbreviate(currentSNPData$varianceCategory)
  # plot the data
  currentSNPPlot<-ggplot(currentSNPData,aes(x=varianceCategory,y=Percentage,fill=varianceCategory)) + 
  geom_bar(stat="identity",position="dodge") + 
  # change the legend title
  labs(fill="Variance Category")+
  # change the x axis title
  xlab("Variance Category")+
  # change the y axis title
  ylab("Percentage")+
  # no image title
  ggtitle("")+
  # make the abbreviations as the x axis labels
  scale_x_discrete(labels=currentSNPData$AbbrevCats)
  
  # save the plot to file
  currentSNPPlotFile<-paste(ResultFolder,pathSeparator,"amovaComponents_",currentSNP,".png",sep="")
  ggsave(currentSNPPlotFile,currentSNPPlot,width=8,height=4)
}

# plot using facet_grid
amovaComponents2plot$AbbrevCats<-abbreviate(amovaComponents2plot$varianceCategory)

FacetPlot<-ggplot(amovaComponents2plot,aes(x=varianceCategory,y=Percentage,fill=varianceCategory)) + 
  geom_bar(stat="identity",position="dodge") + 
  # change the legend title
  labs(fill="Variance Category")+
  # change the x axis title
  xlab("Variance Category")+
  # change the y axis title
  ylab("Percentage")+
  # no image title
  ggtitle("")+
  # make the abbreviations as the x axis labels
  scale_x_discrete(labels=amovaComponents2plot$AbbrevCats)+
  facet_grid(.~SNP)

# save the plot to file
FacetPlotFile<-paste(ResultFolder,pathSeparator,"amovaComponents_Facet.png",sep="")
ggsave(FacetPlotFile,FacetPlot,width=8,height=4)

# plot as multiple lines
multiLine<-ggplot(amovaComponents2plot,aes(x=varianceCategory)) + 
  geom_line(aes(y=Percentage,group=SNP,colour=SNP)) + 
  geom_point(aes(y=Percentage,group=SNP,colour=SNP)) + 
  # change the legend title
  labs(colour="SNP")+
  # change the x axis title
  xlab("Variance Category")+
  # change the y axis title
  ylab("Percentage")+
  # no image title
  ggtitle("")+
  # make the abbreviations as the x axis labels
  scale_x_discrete(labels=amovaComponents2plot$AbbrevCats)

# save the plot to file
multiLineFile<-paste(ResultFolder,pathSeparator,"amovaComponents_multiLine.png",sep="")
ggsave(multiLineFile,multiLine,width=15,height=10)
# as pdf
multiLineFile<-paste(ResultFolder,pathSeparator,"amovaComponents_multiLine.pdf",sep="")
ggsave(multiLineFile,multiLine,width=15,height=10)

# rotate and plot as multiple 
acrossSNPsMultiLine<-ggplot(amovaComponents2plot,aes(x=SNP)) + 
  geom_line(aes(y=Percentage,group=varianceCategory,colour=varianceCategory)) + 
  geom_point(aes(y=Percentage,group=varianceCategory,colour=varianceCategory)) + 
  # change the legend title
  labs(colour="Variance Category")+
  # change the x axis title
  xlab("SNP")+
  # change the y axis title
  ylab("Percentage")+
  # no image title
  ggtitle("")+
  # make the abbreviations as the x axis labels
  scale_x_discrete(labels=amovaComponents2plot$SNP)+
  # rotate the x axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# save the plot to file
acrossSNPsMultiLineFile<-paste(ResultFolder,pathSeparator,"amovaComponents_acrossSNPsMultiLine.png",sep="")
ggsave(acrossSNPsMultiLineFile,acrossSNPsMultiLine,width=10,height=10)
# as pdf
acrossSNPsMultiLineFile<-paste(ResultFolder,pathSeparator,"amovaComponents_acrossSNPsMultiLine.pdf",sep="")
ggsave(acrossSNPsMultiLineFile,acrossSNPsMultiLine,width=10,height=10)

# set.seed(1999)
#monpopsignif   <- randtest(amovacc, nrepet = 999)