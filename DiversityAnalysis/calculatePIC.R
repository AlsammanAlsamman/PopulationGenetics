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

if (length(args) < 2) {
  print(paste("Usage: Rscript calculatePIC.r <vcfFile> <ResultFolder>",sep=""))
  quit()
}

vcfFile<-args[1]
ResultFolder<-args[2]


# load libraries

library(inline)
library(SNPRelate)
library(ggplot2)
library(xlsx)


############################### Functions ##########################################
# calculate PIC in a vector of frequencies in C
PICvec<-cfunction(c(freq="numeric"),"
    // get the length of the vector
    int n = LENGTH(freq);
    // create a 
    double *p = REAL(freq);
    // create a matrix of outer products
    double o[n][n];
    double PIC;

    // Calculate the outer product of p^2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            o[i][j] = p[i] * p[i] * p[j] * p[j];
        }
    }

    // Calculate PIC
    double sum_p_squared = 0.0;
    double sum_upper_tri_o = 0.0;

    for (int i = 0; i < n; i++) {
        sum_p_squared += p[i] * p[i];
        for (int j = i + 1; j < n; j++) {
            sum_upper_tri_o += 2 * o[i][j];
        }
    }
  PIC = 1 - sum_p_squared - sum_upper_tri_o;
  return ScalarReal(PIC);
")

#######################################################################################
# vcfFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Anna_Nassima/5_analysis/Data/BreadGenotype_filtered.vcf"
# ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Anna_Nassima/5_analysis/Results"


# read vcf file using SNPRelate
snpgdsVCF2GDS(vcfFile, paste("Genotype",".gds",sep=""), method="biallelic.only", verbose=TRUE)

# convert to genlight
genofile <- openfn.gds(paste("Genotype",".gds",sep=""))
# filter data to remove loci with missing data and MAF < 0.05


# get genotype
genotype<-read.gdsn(index.gdsn(genofile, "genotype"))
# get marker names
markerNames<-read.gdsn(index.gdsn(genofile, "snp.rs.id"))
# get sample names
sampleNames<-read.gdsn(index.gdsn(genofile, "sample.id"))

# close the gds file
closefn.gds(genofile)

rownames(genotype)<-as.character(sampleNames)
colnames(genotype)<-as.character(markerNames)

# calculate the frequency of each allele in each marker
alleleFreq<-apply(genotype,2,function(x) table(x)/length(x))
# remove sites with only one allele (monomorphic sites)
alleleFreq<-alleleFreq[lapply(alleleFreq,length)>1]
# calculate PIC
PIC<-lapply(alleleFreq,PICvec)
# convert to data frame
PIC<-data.frame(do.call(rbind,PIC))

colnames(PIC)<-c("PIC")
PIC$Marker<-rownames(PIC)
# plot PIC as histogram
p<-ggplot(PIC,aes(x=PIC))+geom_histogram(binwidth = 0.01)+theme_bw()+labs(x="PIC",y="Frequency")
ggsave(paste(ResultFolder,pathSeparator,"PIC_histogram.png",sep=""),p,width=10,height=5)
# as boxplot
p<-ggplot(PIC,aes(x="PIC",y=PIC))+geom_boxplot()+theme_bw()+labs(x="PIC",y="Frequency")
ggsave(paste(ResultFolder,pathSeparator,"PIC_boxplot.png",sep=""),p,width=10,height=5)
# as table
PIC<-PIC[,c("Marker","PIC")]
write.table(PIC,file=paste(ResultFolder,pathSeparator,"PIC.tsv",sep=""),sep="\t",quote=FALSE,row.names = FALSE)
# as xlsx
write.xlsx(PIC,file=paste(ResultFolder,pathSeparator,"PIC.xlsx",sep=""),row.names = FALSE)