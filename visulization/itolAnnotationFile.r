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

# from system get the value of the environment variabl of export AlsammanPopulationGeneticsPath
functionsPath<-Sys.getenv("AlsammanPopulationGeneticsPath")


colorsFile <- "visulization/colorsList.txt"
dataFile <- "/home/samman/Documents/publishing/Lentil/newTree/Data.txt"
outAnnotationFile <- "/home/samman/Documents/publishing/Lentil/newTree/annotation.txt"

# read the colors file
colors <- read.csv(paste0(functionsPath,colorsFile), sep="\t", header = FALSE)
# read the data file
# The data should like this
# sample country/group

data <- read.csv(dataFile, sep="\t")
# add colnames
colnames(data) <- c("sample","country")
# add color as factor of country
data$color <- as.factor(data$country)
data$color<-colors$V1[as.numeric(data$color)]

colors$V1[as.numeric(data$color)]
as.numeric(data$color)
# add range
data$range <- "range"
# reorder the columns
data <- data[,c("sample","range","color","country")]
# replace - or space with .
data$sample <- gsub(" ",".",data$sample)
data$sample <- gsub("-",".",data$sample)

head(data)
#TREE_COLORS			
# SEPARATOR TAB		
# DATA	

# print the annotation file
sink(outAnnotationFile)
cat("TREE_COLORS\n")
cat("SEPARATOR TAB\n")
cat("DATA\n")
write.table(data, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
sink()

