#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

# About checking the phenotypes of the KASP genotypes
# We have validated several SNPs using the KASP genotyping system and now we want to check the phenotypes of the genotypes

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
#######################################################################################
# Libraries ---------------------------------
# loop over the loci
library(reshape2)
library(rstatix)
library(ggplot2)
library(ggpubr)
# for summarise
library(dplyr)
# Functions ---------------------------------
allelefre<-function(snpcol)
{
  # convert ?  to "Uncallable"
  snpcol<-as.character(snpcol)
  snpcol[snpcol=="?"]<-"Uncallable"
  # remove "?"
  snpcol<-snpcol[snpcol!="?"]
  #getthe allele frequency of current locus
  temp<-table(snpcol)
  temp<-temp/sum(temp)
  # fomrat the output to 0.00
  temp<-format(temp, digits=1)
  temp<-as.numeric(temp)
  return(temp)
}
is_empty <- function(x) {
  if (length(x) == 0 & !is.null(x)) {
    TRUE
  } else {
    FALSE
  }
}
checkObersvations<-function(x)
{
  # check observations for alleles
  # if all the observations are the same then return false
  enoughtObs<-TRUE
  alleles <- unique(x$Genotype)
  # if the number of alleles is less than 2 then return false
  if(length(alleles)<2)
  {
    enoughtObs<-FALSE
  }
  for (i in 1:length(alleles))
  {
    if(length(x$Genotype[x$Genotype==alleles[i]])<3)
    {
      enoughtObs<-FALSE
    }
  }
  return(enoughtObs)
}
#######################################################################################

# Program ---------------------------------
setwd("/home/samman/Documents/ICARDA/MaloufKasp")
pathSeparator<-"/"
# File contains the genotypes profile using KASP system
KASPFile<-"/home/samman/Documents/ICARDA/MaloufKasp/Data/genotype.csv"
# phenotype file
phenotypeFile<-"/home/samman/Documents/ICARDA/MaloufKasp/Data/phenotype.tsv"
ResultFolder<-"/home/samman/Documents/ICARDA/MaloufKasp/Results"
# create the folder if it does not exist
if(!dir.exists(ResultFolder))
{
  dir.create(ResultFolder)
}

# read the file
KASP <- read.csv(KASPFile, header = TRUE, sep = ",", row.names = 1)  

# read the phenotype file
phenotype<-read.csv(phenotypeFile, header = TRUE, sep = "\t", row.names = 1)
# filter the phenotype file
phenotype.rownames<-rownames(phenotype)
phenotype.colnames<-colnames(phenotype)
# convert non numeric to NA
phenotype<-apply(phenotype,2,function(x) as.numeric(as.character(x)))
# fill the NA with the mean of the column
phenotype<-apply(phenotype,2,function(x) ifelse(is.na(x),mean(x,na.rm=TRUE),x))
# convert negative values to NA
# print warning if there is negative values
if(length(phenotype[phenotype<0])>0)
{
  print("Warning: Negative values were found in the phenotype file")
}
phenotype<-apply(phenotype,2,function(x) ifelse(x<0,NA,x))
# convert the phenotype to a data frame
phenotype<-data.frame(phenotype)
# add the row and column names
rownames(phenotype)<-phenotype.rownames
colnames(phenotype)<-phenotype.colnames

# remove the loci with only one allele
KASP.filtered<-KASP[,apply(KASP,2,function(x) length(unique(x)))>1]

############ Analyse the genotypes
# get the name of the phenotypes
phenotypesNames<-colnames(phenotype)
lociNames<-colnames(KASP.filtered)
# merge genotype and phenotype
PhenoGeno<-merge(KASP.filtered,phenotype,by="row.names",all.x=TRUE)
# add the row names
rownames(PhenoGeno)<-PhenoGeno[,1]
# remove the first column
PhenoGeno<-PhenoGeno[,-1]
# create folders for snps
# create the folder if it does not exist
# for(loci in lociNames)
# {
#   if(!dir.exists(paste0("Results",pathSeparator,loci)))
#   {
#     dir.create(paste0("Results",pathSeparator,loci))
#   }
# }
# loop over the phenotypes and loci
# get the significant loci
sigLoci<-ls()
for(i in 1:length(phenotypesNames))
{
  # get the name of the phenotype
  phenoName<-phenotypesNames[i]
  for (j in 1:length(lociNames))
  {
    # get the name of the locus
    lociName<-lociNames[j]
    # select the locus and phenotype
    phenotemp<-PhenoGeno[,c(lociName,phenoName)]
    colnames(phenotemp)<-c("Genotype","Phenotype")
    # remove rows with ? as genotype
    phenotemp<-phenotemp[phenotemp$Genotype!="?",]
    # check if there is enough observations
    if(!checkObersvations(phenotemp))
    {
      next
    }
    stat.test <- phenotemp %>% t_test(Phenotype ~ Genotype)
    stat.test <- stat.test %>% add_xy_position(x = "Genotype")
    
    # is there a significant difference between alleles
    if(length(stat.test$p.adj.signif[stat.test$p.adj.signif<1])>0)
    {
      print("Significance were found")
      # save the locus name
      sigLoci<-c(sigLoci,lociName)
      
      # plot the results
      plotTitle<-paste0("Phenotype: ",phenoName," Locus: ",lociName)
      bxp <- ggboxplot(phenotemp, x = "Genotype", y = "Phenotype", fill = "Genotype", 
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"))
      bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)
      # add the title
      # bxp<-bxp+labs(title=plotTitle)
      # add the number of observations to the plot
      # bxp<-bxp+
      # annotate("text",
      #      x = 1:length(table(phenotemp$Genotype)),
      #      y = aggregate(phenotemp$Phenotype, by = list(phenotemp$Genotype), FUN = mean)$x,
      #      label = table(phenotemp$Genotype),
      #      col = "red",
      #      vjust = - 1)  
      # save the plot
      ggsave(paste0(phenoName,"_",lociName,".png"), width = 10, height = 10)
    }
  }
}
# table the significant loci
sigLociTable<-table(sigLoci)
# data frame of the significant loci
sigLociTable<-data.frame(sigLociTable)
colnames(sigLociTable)<-c("Locus","Count")
# sort the table
sigLociTable<-sigLociTable[order(sigLociTable$Count,decreasing = TRUE),]
# save the table
write.csv(sigLociTable,file=paste0(ResultFolder,pathSeparator,"sigLociTable.csv"),row.names = FALSE)
