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


# load libraries
library(ggplot2)
# install.packages("UpSetR")
library(UpSetR)
library(pheatmap)

#TargetFolder<-args[1]
#fdrcutoff<-args[2]
#minFreq<-args[3]
#outFolder<-args[4]

TargetFolder<-"/home/samman/Documents/publishing/ewas/Paper/Results/ewas_tables"
fdrcutoff<-0.06
minFreq<-3
outFolder<-"/home/samman/Documents/publishing/ewas/Paper/Results/ewas_tables/venn"
pathSeparator<-"/"

# create the output folder if it does not exist
if (!dir.exists(outFolder)) {
  dir.create(outFolder)
}


numberOfFiles<-length(list.files(TargetFolder, pattern = "*.tsv", full.names = TRUE))

filterCol<-"fdr"
# create a list to store the data
ListInput<-list()
# loop over the files in the folder
for (file in list.files(TargetFolder, pattern = "*.tsv", full.names = TRUE)) {
  # file basename
  traitName<-basename(file)
  traitName<-gsub(".tsv","",traitName)
  # read the file
  data<-read.csv(file,header=T,sep="\t",stringsAsFactors = F)
  #add a column for the chromosome and position
  data$rs<-paste(data$chr,data$pos,sep=":")
  # calculate fdr using pvalues column
  data<-data[data[,filterCol]<=fdrcutoff,]
  # add the data to the list
  ListInput[[traitName]]<-data$rs
}


# plot the venn diagram
png(paste(outFolder,pathSeparator,"venn.png",sep=""),width = 10,height = 10, units = "in", res = 300)
upset(fromList(ListInput),nsets = numberOfFiles)
dev.off()
 
# plot the venn diagram
pdf(paste(outFolder,pathSeparator,"venn.pdf",sep=""),width = 10,height = 10)
upset(fromList(ListInput),nsets = numberOfFiles)
dev.off()

# get intersection from list
# convert the list to a data frame
oneList<-unlist(ListInput)
# count the number of occurrences
oneList<-table(oneList)
oneList<-as.data.frame(oneList)
# remove rows with count less than 2
oneList<-oneList[oneList$Freq>=minFreq,]

# create a matrix of the rs ids
rsMatrix<-matrix(nrow = nrow(oneList),ncol = numberOfFiles)
# loop over the list and fill the matrix
for (i in 1:nrow(oneList)) {
  # get the rs id
  rs<-oneList[i,1]
  # loop over the list
  for (j in 1:numberOfFiles) {
    # get the trait name
    traitName<-names(ListInput)[j]
    # get the rs ids
    rsIds<-ListInput[[traitName]]
    # check if the rs id is in the list
    if (rs %in% rsIds) {
      # add the fdr
      rsMatrix[i,j]<-1
    } else {
      rsMatrix[i,j]<-0
    }
  }
}

# add row names
rownames(rsMatrix)<-oneList[,1]
# add column names
colnames(rsMatrix)<-names(ListInput)


# create Table for sum of rows
rsMatrixFreq<-cbind(rowSums(rsMatrix))
colnames(rsMatrixFreq)<-"Freq"
rsMatrixFreq<-as.data.frame(rsMatrixFreq)
rsMatrixFreq$Freq<-as.factor(rsMatrixFreq$Freq)


# get chr information from marker name using chr
rsMatrixChr<-gsub(":.*","",rownames(rsMatrix))
# create chrannotation
chranno <- data.frame( chr = rsMatrixChr)
rownames(chranno)<-rownames(rsMatrix)
# merge with rsMatrixFreq
rsMatrixFreq<-merge(rsMatrixFreq,chranno,by="row.names")
rownames(rsMatrixFreq)<-rsMatrixFreq[,1]
rsMatrixFreq<-rsMatrixFreq[,-1]
rsMatrixFreq$chr<-paste("chr",rsMatrixFreq$chr,sep="")

png(paste(outFolder,pathSeparator,"vennMatrix.png",sep=""),width = 10,height = 10, units = "in", res = 300)
pheatmap(rsMatrix, 
        annotation_row = rsMatrixFreq,
        legend_breaks = c(0, 1),
        legend_labels = c("absent", "present"), 
        color = colorRampPalette(c("white", "black"))(2),
        # smaller row font size
        fontsize_row = 4)
dev.off()

# without row names
png(paste(outFolder,pathSeparator,"vennMatrixNoRowNames.png",sep=""),width = 10,height = 10, units = "in", res = 300)
pheatmap(rsMatrix, 
        annotation_row = rsMatrixFreq,
        legend_breaks = c(0, 1),
        legend_labels = c("absent", "present"), 
        color = colorRampPalette(c("white", "black"))(2),
        show_rownames = FALSE)
dev.off()

# pdf 
pdf(paste(outFolder,pathSeparator,"vennMatrix.pdf",sep=""),width = 10,height = 20)
pheatmap(rsMatrix, 
        annotation_row = rsMatrixFreq,
        legend_breaks = c(0, 1),
        legend_labels = c("absent", "present"), 
        color = colorRampPalette(c("white", "black"))(2))
dev.off()

# without row names
pdf(paste(outFolder,pathSeparator,"vennMatrixNoRowNames.pdf",sep=""),width = 10,height = 10)
pheatmap(rsMatrix, 
        annotation_row = rsMatrixFreq,
        legend_breaks = c(0, 1),
        legend_labels = c("absent", "present"), 
        color = colorRampPalette(c("white", "black"))(2),
        show_rownames = FALSE)
dev.off()
# save the matrix as a table
write.table(rsMatrix,paste(outFolder,pathSeparator,"vennMatrix.txt",sep=""),sep="\t",quote = F,row.names = T,col.names = T)

# save freq table
write.table(rsMatrixFreq,paste(outFolder,pathSeparator,"vennMatrixFreq.txt",sep=""),sep="\t",quote = F,row.names = T,col.names = T)
