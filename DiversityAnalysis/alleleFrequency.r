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

################################## Libraries ############################################

library(SNPRelate)
library(reshape2)
library(ggplot2)

#######################################################################################
# The allele SNP frequency across the population
# The input is the vcf file and the SNP allele list 

# vcfFile = args[1]
# SNPAlleleList = args[2]
# ResultFolder = args[3]

vcfFile = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/GenotypeData_filtered.vcf"
SNPAlleleListFile = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/Data/TargetSNPs_test.txt"
ResultFolder = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS/AlleleFreq"

#  read SNP list
SNPList<-read.table(SNPAlleleListFile,header=FALSE,sep="\t")
SNPList

# Read VCF file
snpgdsVCF2GDS(vcfFile, paste("GenoTypeData",".gds",sep=""), method="biallelic.only", verbose=TRUE)
genofile <- openfn.gds(paste("GenoTypeData",".gds",sep=""))

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
#delete the gds file
file.remove(paste("GenoTypeData",".gds",sep=""))
# select for first marker
markers<- SNPList$V1
# loop over markers
# Allele frequency total Table
alleleFreqTotal<-data.frame(matrix(ncol = 4, nrow = 0))
colnames(alleleFreqTotal)<-c("Marker","Allele","Frequency","AlleleValue")

for (m in 1:length(markers)) 
{
  # check if the marker is in the genotype
  if (!(markers[m] %in% colnames(genotype))) {
    print(paste("Marker ",markers[m]," is not in the genotype",sep=""))
    next
  }
  print(paste("Plotting for marker ",markers[m],sep=""))
  # get genotype for the marker
  markerkGenotype<-genotype[,markers[m]]
  # get marker index
  markerIndex<-which(colnames(genotype)==markers[m])
  # get marker allele
  markerAllele<-snp.allele[markerIndex]
  allele0Genotypes<-names(markerkGenotype[markerkGenotype==0])
  allele1Genotypes<-names(markerkGenotype[markerkGenotype==1])
  allele2Genotypes<-names(markerkGenotype[markerkGenotype==2])
  # save genotypes for each allele in text file
  write.table(allele0Genotypes,paste(ResultFolder,pathSeparator,markers[m],"_allele0.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(allele1Genotypes,paste(ResultFolder,pathSeparator,markers[m],"_allele1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(allele2Genotypes,paste(ResultFolder,pathSeparator,markers[m],"_allele2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)


  # get the allele frequency
  allele0Freq<-length(allele0Genotypes)/length(sample.id)
  allele1Freq<-length(allele1Genotypes)/length(sample.id)
  allele2Freq<-length(allele2Genotypes)/length(sample.id)
  # get the allele value
  # split the allele
  alleles<-strsplit(markerAllele,"/")
  allele0Value<-paste(alleles[[1]][1],alleles[[1]][1],sep="")
  allele1Value<-paste(alleles[[1]][1],alleles[[1]][2],sep="")
  allele2Value<-paste(alleles[[1]][2],alleles[[1]][2],sep="")
  
  # plot the allele frequency
  alleleFreq<-data.frame(Marker=rep(markers[m],3),Allele=c("0","1","2"),Frequency=c(allele0Freq,allele1Freq,allele2Freq),AlleleValue=c(allele0Value,allele1Value,allele2Value)) 
  # add to the total table
  alleleFreqTotal<-rbind(alleleFreqTotal,alleleFreq)
  
  # plot the allele frequency
  plotAllele<- ggplot(alleleFreq,aes(x=AlleleValue, y=Frequency, fill=AlleleValue)) +
    geom_bar(stat="identity") +
    theme_bw()+
    ggtitle(paste("Allele Frequency for Marker ",markers[m],sep=""))+
    xlab("Allele")+
    ylab("Frequency")+
    theme(plot.title = element_text(hjust = 0.5))
  # save the plot
  ggsave(paste(ResultFolder,pathSeparator,markers[m],".png",sep=""),plot=plotAllele,width=5,height=5)
}
# plot marker allele frequency
# sa
alleleFreqTotal
# facet_grid
plotAlleles<-ggplot(alleleFreqTotal, aes(x=AlleleValue, y=Frequency, fill=AlleleValue)) +
  geom_bar(stat="identity") +
  facet_grid(.~Marker) +
  theme_bw()+
  ggtitle("Allele Frequency for all Markers")+
  xlab("Allele")+
  ylab("Frequency")+
  theme(plot.title = element_text(hjust = 0.5))

# save the plot
ggsave(paste(ResultFolder,pathSeparator,"AlleleFrequency.png",sep=""),plot=plotAlleles,width=10,height=5)







