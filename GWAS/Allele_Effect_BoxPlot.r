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


if (length(args) < 2) {
  print(paste("Usage: Allele_Effect_BoxPlot.r <vcffile> <phenotyfile> <snplist> <ResultFolder>",sep=""))
  quit()
}

# loading libraries
library(SNPRelate)
library(reshape2)
library(ggplot2)

# # read arguments
# VCFFILE<-args[1]
# SNPListFile<-args[2]
# traitDataFile<-args[3]
# ResultFolder<-args[4]
# sepPlots<-args[5]
# if not specified then set to FALSE
if (is.na(sepPlots)) {
  sepPlots<-FALSE
}

VCFFILE<-"/home/samman/Documents/publishing/Khaled_Fausarium/NewAnalysis/Data/FW_filtered.vcf.recode.vcf"
SNPListFile<-"/home/samman/Documents/publishing/Khaled_Fausarium/NewAnalysis/1_vcf2gwasout/TargetSNPs.txt"
traitDataFile<-"/home/samman/Documents/publishing/Khaled_Fausarium/NewAnalysis/Data/Phenotype.tsv"
ResultFolder<-"/home/samman/Documents/publishing/Khaled_Fausarium/NewAnalysis/1_vcf2gwasout/boxplot"

#setwd("/home/samman/Documents/publishing/Drought_Ahmed_Master/NewAnalysis/someResults")

#  read SNP list
SNPList<-read.table(SNPListFile,header=FALSE,sep="\t")
# Read VCF file
snpgdsVCF2GDS(VCFFILE, paste("BreadGenotype",".gds",sep=""), method="biallelic.only", verbose=TRUE)
genofile <- openfn.gds(paste("BreadGenotype",".gds",sep=""))

# get sample.id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# snp.allele
snp.allele <- read.gdsn(index.gdsn(genofile, "snp.allele"))
# snp.rs.id
snp.rs.id  <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
# genotype
genotype <- read.gdsn(index.gdsn(genofile, "genotype"))
rownames(genotype) <- as.character(sample.id)
colnames(genotype) <- as.character(snp.rs.id)
closefn.gds(genofile)
colnames(genotype)

# read trait data
traitData<-read.table(traitDataFile,header=TRUE,sep="\t", row.names=1)
rownames(traitData)<-as.character(rownames(traitData))
# select phenotypes values for each allele
traitData$taxa<-rownames(traitData)
# melt the data
traitData<-melt(traitData,id.vars="taxa")

# select only the SNPs in the SNP list
genotype<-genotype[,colnames(genotype) %in% as.character(SNPList$V1)]
traitData

# select for first marker
markers<- SNPList$V1
phenoGenoData.total<-data.frame(matrix(ncol = 5, nrow = 0))
# phenotype<-as.character(traitData$variable[1])
for (m in markers)
{
  print(paste("Plotting for marker ",m,sep=""))
  # get genotype for the marker
  genotype1<-genotype[,m]
  # get marker index
  markerIndex<-which(colnames(genotype)==m)
  # get marker allele
  markerAllele<-snp.allele[markerIndex]
  # select alelles
  allele0Genotypes<-names(genotype1[genotype1==0])
  allele1Genotypes<-names(genotype1[genotype1==1])
  allele2Genotypes<-names(genotype1[genotype1==2])

  # does marker name contain ">"
  if (grepl(">",m)) {
    # alleles values AA, AB, BB
    # get alleles contains letter > letter
    allelestrlist <- unlist(strsplit(m,":"))
    allelestr <- allelestrlist[length(allelestrlist)]
    alleles<-unlist(strsplit(allelestr,">"))
    # alleles<-unlist(strsplit(markerAllele,">"))
    allele0<-paste(alleles[1],alleles[1],sep="")
    allele1<-paste(alleles[1],alleles[2],sep="")
    allele2<-paste(alleles[2],alleles[2],sep="")
  } else
  {
    # alleles values AA, AB, BB
    alleles<-unlist(strsplit(markerAllele,"/"))
    allele0<-paste(alleles[1],alleles[1],sep="")
    allele1<-paste(alleles[1],alleles[2],sep="")
    allele2<-paste(alleles[2],alleles[2],sep="")
  }

  # create a table for the three alleles
  allelesTable<-data.frame(Genotype = allele0Genotypes, allele=rep(allele0,length(allele0Genotypes)))
  # add allele1
  allelesTable<-rbind(allelesTable,data.frame(Genotype = allele1Genotypes, allele=rep(allele1,length(allele1Genotypes))))
  # add allele2
  allelesTable<-rbind(allelesTable,data.frame(Genotype = allele2Genotypes, allele=rep(allele2,length(allele2Genotypes))))

  # merge with trait data
  phenoGenoData<-merge(traitData,allelesTable,by.x="taxa",by.y="Genotype")
  # plot marker effect
  
  p<-ggplot(phenoGenoData, aes(x=allele, y=value, fill=allele)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  # use facet_wrap to plot each phenotype in a separate panel
  facet_wrap(~variable, ncol=10) +
  # use theme_bw() to remove the gray background
  theme_bw() + 
  # add a title
  ggtitle(paste("Marker :",m,sep=""))+ 
  # no x-axis title
  xlab("") 
  ggsave(paste(ResultFolder,pathSeparator,m,".png",sep=""),p, width = 10, height = 2)
  # add to total data
  phenoGenoData$marker<-m
  phenoGenoData.total<-rbind(phenoGenoData.total,phenoGenoData)
}


# # plot all markers in one plot
# p<-ggplot(phenoGenoData.total, aes(x=allele, y=value, fill=allele)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2) +
#   # use facet_wrap to plot each phenotype in a separate panel
#   facet_grid(~marker+variable) +
#   # use theme_bw() to remove the gray background
#   theme_bw() + 
#   # no x-axis title
#   xlab("") 
#   ggsave(paste(ResultFolder,pathSeparator,"AllMarkers.png",sep=""),p, width = 20, height = 10)


