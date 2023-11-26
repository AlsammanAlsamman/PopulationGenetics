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

if (length(args) < 3) {
  print(paste("Usage: Rscript BasicStats.r <vcfFile> <metadataFile> <ResultFolder>",sep=""))
  quit()
}

vcfFile<-args[1]
metadataFile<-args[2]
ResultFolder<-args[3]

# load libraries
#install.packages("dartR")
#install.packages("xlsx")
library(dartR)
library(adegenet)
library('vcfR')
library('poppr')
library(reshape2)
library(xlsx)
library(ggplot2)
# library(devtools)
# install_github("nikostourvas/PopGenUtils")
library("PopGenUtils")

# vcfFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Anna_Nassima/5_analysis/Data/BreadGenotype.vcf"
# ResultFolder<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Anna_Nassima/5_analysis/Results"

# log file name with date and time
logFileName<-paste0(ResultFolder,pathSeparator,"BasicStats_",Sys.Date(),"_",Sys.time(),".log")

# read metadata
metaData <- read.csv( metadataFile , sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names = 1)
# read using vcfR and convert to genlight
vcfData <- read.vcfR(vcfFile)
# convert to genlight
gdata.gen <- vcfR2genlight(vcfData)

# filter data to remove loci with missing data
gdata.gen<-gl.filter.allna(gdata.gen)

# approved Cols with more than 1 value
approvedCols<-apply(metaData,2,function(x) length(unique(x))>1)
# print the name of non approved columns
sink(logFileName, append=TRUE)
print(paste("The following columns have only one value and will be ignored:",paste(names(approvedCols[approvedCols==FALSE]),collapse=", ")))
sink()

# select only approved columns
approvedCols<-names(approvedCols[approvedCols==TRUE])
# for each column calculate basic stats

# col Stats
colStats<-list()

samplenames<-gdata.gen@ind.names
for(metaCol in approvedCols)
{
  print(paste("Calculating basic stats for column:",metaCol))
  pop(gdata.gen) <- as.factor(metaData[samplenames, metaCol])
  # calculate basic stats
  genBasicStats<-gl.basic.stats(gdata.gen)
  # save to list
  colStats[[metaCol]]<-genBasicStats
}

# save the object to file for later use
save(colStats,file=paste0(ResultFolder,pathSeparator,"BasicStats_",Sys.Date(),"_",Sys.time(),".RData"))

# Statistics table across all columns
statsTable<-data.frame(matrix(NA,nrow = length(colStats),ncol = 10))
rownames(statsTable)<-names(colStats)
colnames(statsTable)<-c("Ho","Hs","Ht","Dst","Htp","Dstp","Fst","Fstp","Fis","Dest")
for (metaCol in names(colStats))
{
  print(paste("Plotting basic stats for column:",metaCol))
  # statsTable[metaCol,]<-rbind(colStats[[metaCol]]$overall)
  statsTable[metaCol,]<-colStats[[metaCol]]$overall  
}
# save to tsv file  
write.table(statsTable,file=paste0(ResultFolder,pathSeparator,"BasicStats_",Sys.Date(),"_",Sys.time(),".tsv"),sep="\t",quote=FALSE)
# save as xlsx file
write.xlsx(statsTable,file=paste0(ResultFolder,pathSeparator,"BasicStats.xlsx"),row.names = TRUE)

# create a folder for plots basicStatsPlot if not exists
if (!dir.exists(paste0(ResultFolder,pathSeparator,"basicStatsPlot"))) {
  dir.create(paste0(ResultFolder,pathSeparator,"basicStatsPlot"))
  # create a pdf folder
  dir.create(paste0(ResultFolder,pathSeparator,"basicStatsPlot",pathSeparator,"pdf"))
}
# plot the basic stats for each column
for (metaCol in names(colStats))
{
  print(paste("Plotting basic stats for column:",metaCol))
  # get the stats
  genBasicStats<-colStats[[metaCol]]
  # get the perloc stats
  perLocStats<-genBasicStats$perloc
  perLocStats$loc<-rownames(perLocStats)
  # melt
  perLocStats<-melt(perLocStats, id.vars = "loc")
  # remove inf values
  perLocStats<-perLocStats[!is.infinite(perLocStats$value),]  
  # remove NA values
  perLocStats<-perLocStats[!is.na(perLocStats$value),]
  # plot
  p<-ggplot(perLocStats, aes(x=variable, y=value, fill=variable)) + 
    geom_boxplot()+ 
    # geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(metaCol)+
    # no legend
    guides(fill=FALSE) +
    # x title
    ylab("Values per loci")+
    # y title
    xlab("Genetic Index")+
    # no plot title
    theme(plot.title = element_blank())+
    # bigger font size
    theme(text = element_text(size=15))
 
  # save to pdf folder
  ggsave(paste0(ResultFolder,pathSeparator,"basicStatsPlot",pathSeparator,"pdf",pathSeparator,"BasicStats_",metaCol,".pdf"),p,width = 10, height = 10)
  # png
  ggsave(paste0(ResultFolder,pathSeparator,"basicStatsPlot",pathSeparator,"BasicStats_",metaCol,".png"),p,width = 10, height = 10)
}


# create a folder for plots acrossChromosome if not exists
if (!dir.exists(paste0(ResultFolder,pathSeparator,"acrossChromosome"))) {
  dir.create(paste0(ResultFolder,pathSeparator,"acrossChromosome"))
  # create a pdf folder
  dir.create(paste0(ResultFolder,pathSeparator,"acrossChromosome",pathSeparator,"pdf"))
}

# plot the basic stats for each column according to chromosome
for (metaCol in names(colStats))
{
  print(paste("Plotting basic stats for column:",metaCol))
  # get the stats
  genBasicStats<-colStats[[metaCol]]
  # get the perloc stats
  perLocStats<-genBasicStats$perloc
  perLocStats$loc<-rownames(perLocStats)
  # melt
  perLocStats<-melt(perLocStats, id.vars = "loc")
  # remove inf values
  perLocStats<-perLocStats[!is.infinite(perLocStats$value),]  
  # remove NA values
  perLocStats<-perLocStats[!is.na(perLocStats$value),]
  names(genBasicStats)
  # position info
  positionInfo<-data.frame(loc=gdata.gen@loc.names,chromosome=gdata.gen@chromosome,position=gdata.gen@position)
  # merge
  perLocStats<-merge(perLocStats,positionInfo,by="loc")
  # plot
  p<-ggplot(perLocStats, aes(x=variable, y=value, fill=variable)) + 
    geom_boxplot()+ 
    # geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(metaCol)+
    # no legend
    guides(fill=FALSE) +
    # x title
    ylab("Values per loci")+
    # y title
    xlab("Genetic Index")+
    # no plot title
    theme(plot.title = element_blank())+
    # bigger font size
    theme(text = element_text(size=15))+
    # facet according to chromosome
    facet_wrap(~chromosome, ncol = 2)
  # save to pdf folder
  ggsave(paste0(ResultFolder,pathSeparator,"acrossChromosome",pathSeparator,"pdf",pathSeparator,"BasicStats_",metaCol,".pdf"),p,width = 20, height = 30)
  # png
  ggsave(paste0(ResultFolder,pathSeparator,"acrossChromosome",pathSeparator,"BasicStats_",metaCol,".png"),p,width = 20, height = 30)
}


# calculate the MAF for each locus
MAFvalues<-maf(vcfData)
MAFvalues<-as.data.frame(MAFvalues)
# plot the MAF
p<-ggplot(MAFvalues, aes(x=Frequency)) + 
  geom_histogram(binwidth = 0.01)+
  ggtitle("MAF")+
  # no legend
  guides(fill=FALSE) +
  # x title
  ylab("Frequency")+
  # y title
  xlab("MAF")+
  # no plot title
  theme(plot.title = element_blank())+
  # bigger font size
  theme(text = element_text(size=15))
# save to pdf folder
ggsave(paste0(ResultFolder,pathSeparator,"MAF.pdf"),p,width = 10, height = 10)
# png
ggsave(paste0(ResultFolder,pathSeparator,"MAF.png"),p,width = 10, height = 10)
