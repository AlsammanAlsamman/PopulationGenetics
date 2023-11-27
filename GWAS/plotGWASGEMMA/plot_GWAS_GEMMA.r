#!/usr/bin/env Rscript

# plot LEA output of the ewas analysis
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpubr)
# from system get the value of the environment variabl of export AlsammanPopulationGeneticsPath
functionsPath<-Sys.getenv("AlsammanPopulationGeneticsPath")
functionsPath<-paste(functionsPath,"GWAS/plotGWASGEMMA/functions",sep="")
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

# use args to get the input and output files
args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 3) {
  print(paste("Usage: Rscript plot_GWAS_GEMMA.r <gwasfile> <traitName> <ResultFolder> <fdrcutoff> <thresholdType> <qthreshold> <chrTable>",sep=""))
  quit()
}

gwasfile<-args[1]
traitName<-args[2]
ResultFolder<-args[3]


# default fdr cutoff
fdrcutoff<-0.1
# default qthreshold
qthreshold<-0.001

# ChromsomesTable
chrTable<-""

if (length(args) == 4) {
  fdrcutoff<-as.numeric(args[4])
}

thresholdType<-"quantile"

# if args[5] is not empty then use it as the threshold type
if (args[5]!="") {
  thresholdType<-args[5]
}
if (thresholdType=="quantile" & length(args) == 6) {
   qthreshold<-as.numeric(args[6])
}
  
if (length(args) == 7) {
   chrTable<-args[7] 
}

# default fdr cutoff
if (length(args) != 4) {
  
  print(paste("The default FDR cutoff is ",fdrcutoff,sep=""))
}
# default qthreshold
if (length(args) != 5) {
  
  print(paste("The default qthreshold is ",qthreshold,sep=""))
}

# create the output folder if it does not exist
if (!file.exists(ResultFolder)) {
  dir.create(ResultFolder)
  print(paste("The folder ",ResultFolder," was created",sep=""))
}


# gwasfile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/Field/gwas_out/gy_BLUE_mod_sub_Field.converted.part5_Field_filtered.assoc.txt"
# traitName<-"Samman"
# ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Amounane_24_11_2023/RESULTS/test"


# The log file
logFile<-paste(ResultFolder,"/",traitName,".log",sep="")
# read the output of the gwas analysis
gwasData<-read.csv(gwasfile,header=TRUE,sep="\t")

# print to the log file
sink(logFile, append = FALSE, split = FALSE)
print(paste("The number of SNPs is ",nrow(gwasData),sep=""))
sink()

# calculate the FDR
gwasData$fdr<-p.adjust(gwasData$p_wald, method = "BH")
# calculate the threshold
threshold<-quantile(gwasData$p_wald, probs = qthreshold, na.rm = TRUE)

# select the log global value
log.global<- NA

loggvalues<-gwasData[gwasData$p_wald<threshold,]$p_wald
if (length(loggvalues)!=0) {
  loggvalues<--log10(loggvalues)
  log.global<-min(loggvalues)
} 

if(thresholdType!="quantile" & length(args) > 6) {
   log.global<-as.numeric(args[6]) 
   print(paste("The log global value is ",log.global,sep=""))
}

# select the logmin value according to the FDR
logSignvalues<-gwasData[gwasData$fdr<fdrcutoff,]$p_wald
log.Sign<-NA
if (length(logSignvalues)!=0) {
  logSignvalues<--log10(logSignvalues)
  log.Sign<-min(logSignvalues)
} 

sink(logFile, append = TRUE, split = FALSE)
# check the logmin value
if (is.na(log.Sign)) {
  print(paste("No significant SNPs with FRD < ",fdrcutoff," were found",sep=""))
}

# select the significant SNPs according to the logmin
SNPsFDR<-gwasData[gwasData$fdr<fdrcutoff,]

# print to the log file
print(paste("The number of significant SNPs with FDR < ",fdrcutoff," is ",nrow(SNPsFDR),sep=""))
sink()

# sort the data by p-value
gwasData<-gwasData[order(gwasData$p_wald),]
# write the data to a file
write.table(gwasData, file = paste(ResultFolder,"/",traitName,".tsv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)

# write the FDR SNPs to a file if there are any
if (nrow(SNPsFDR)>0) {
  # select the significant SNPs according to the fdr
  write.table(SNPsFDR, file = paste(ResultFolder,"/",traitName,"_FDR.tsv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
}

# write the SNPs with p-value < threshold to a file if there are any
if (nrow(gwasData[gwasData$p_wald<threshold,])>0) {
  # select the significant SNPs according to the fdr
  write.table(gwasData[gwasData$p_wald<threshold,], file = paste(ResultFolder,"/",traitName,"_threshold.tsv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
}

# convert NA chr to Un
gwasData$chr[is.na(gwasData$chr)]<-"Un"
# if there is a chromosome table then read it
if (chrTable!="") {
  chrTable<-read.table(chrTable,sep="\t")
  # set the column names
  colnames(chrTable)<-c("chr","chrName")
  # replace the chromosome names according to the table
  for (i in 1:nrow(chrTable)) {
    gwasData[gwasData$chr==chrTable$chr[i],]$chr<-chrTable$chrName[i]
  }
  # sort the data by chromosome
  gwasData$chr <- factor(gwasData$chr, levels = chrTable$chrName)
  print(paste("The chromosome table was used to replace the chromosome names",sep=""))
}

#5826041|F|0-30:G>A-30:G>A 5825086 8776047|F|0-63:C>T-63:C>T 11063413|F|0-57:T>C-57:T>C  5825195 5825384|F|0-47:T>C-47:T>C 10268227|F|0-54:C>T-54:C>T  10264793|F|0-44:C>G-44:C>G  10269791|F|0-26:C>T-26:C>T  10264807|F|0-44:G>C-44:G>C  5825194|F|0-9:A>T-9:A>T
#dashSNPs<-c("5826041|F|0-30:G>A-30:G>A","8776047|F|0-63:C>T-63:C>T","11063413|F|0-57:T>C-57:T>C","5825384|F|0-47:T>C-47:T>C","10268227|F|0-54:C>T-54:C>T","10264793|F|0-44:C>G-44:C>G","10269791|F|0-26:C>T-26:C>T","10264807|F|0-44:G>C-44:G>C","5825194|F|0-9:A>T-9:A>T","5825086")

dashSNPs <- c("5826041|F|0-30:G>A-30:G>A"
,"5825086"
,"8776047|F|0-63:C>T-63:C>T"
,"11063413|F|0-57:T>C-57:T>C"
,"5825195"
,"5825384|F|0-47:T>C-47:T>C"
,"10268227|F|0-54:C>T-54:C>T"
,"10264793|F|0-44:C>G-44:C>G"
,"10269791|F|0-26:C>T-26:C>T"
,"10264807|F|0-44:G>C-44:G>C"
,"5825194|F|0-9:A>T-9:A>T")


p<-plotmanahhten(gwasData,
                logmin = log.Sign,
                fdrcutoff = fdrcutoff,
                log.global = log.global,
                annotateSignf = T,
                TraitName=traitName,
                fdrSignSNPs = SNPsFDR,
                matchChrTable = T,
                dash.this = dashSNPs
                )


# plot qqplot with random colors
cols<-sample(c("red","blue","green","black","orange","purple","brown","cyan","magenta","yellow"))[1]

qq<-plotqq(gwasData,col=cols)
# merge the two plots
qqpp<-ggarrange(p, qq, ncol = 2, nrow = 1, widths = c(3, 1))

# create folder for the qqplot if it does not exist
if (!file.exists(paste(ResultFolder,"/qqplot",sep=""))) {
  dir.create(paste(ResultFolder,"/qqplot",sep=""))
  print(paste("The folder ",paste(ResultFolder,"/qqplot",sep="")," was created",sep=""))
}
# create folder for the ppplot if it does not exist
if (!file.exists(paste(ResultFolder,"/ppplot",sep=""))) {
  dir.create(paste(ResultFolder,"/ppplot",sep=""))
  print(paste("The folder ",paste(ResultFolder,"/ppplot",sep="")," was created",sep=""))
}

# create folder for pdf if it does not exist
if (!file.exists(paste(ResultFolder,"/pdf",sep=""))) {
  dir.create(paste(ResultFolder,"/pdf",sep=""))
  print(paste("The folder ",paste(ResultFolder,"/pdf",sep="")," was created",sep=""))
}

# # pdf too
ggsave(p, filename = paste(ResultFolder,"/pdf/",traitName,".pdf",sep=""), width = 10, height = 5, dpi = 300)
ggsave(qqpp, filename = paste(ResultFolder,"/pdf/",traitName,"_qqpp.pdf",sep=""), width = 15, height = 4, dpi = 300)
ggsave(qq, filename = paste(ResultFolder,"/pdf/",traitName,"_qq.pdf",sep=""), width = 15, height = 5, dpi = 300)


# # png
# merged plot
ggsave(qqpp, filename = paste(ResultFolder,"/",traitName,"_qqpp.png",sep=""), width = 10, height = 3, dpi = 300, device = "png")
# # separate plots
ggsave(p, filename = paste(ResultFolder,"/ppplot/",traitName,".png",sep=""), width = 10, height = 3, dpi = 300)
ggsave(qq, filename = paste(ResultFolder,"/qqplot/",traitName,"_qq.png",sep=""), width = 15, height = 4, dpi = 300, device = "png")


