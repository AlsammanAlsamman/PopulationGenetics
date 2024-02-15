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

# The idea is to make only representative values where values with frequency more than 20 are considered

dataFile  <- "/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Morocco_meeting_final/Data/passportData.order.tsv"
data <- read.table(dataFile, header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)

# loop through the columns and make a representative value

for (i in 1:ncol(data)) {
  # table the column
  table <- table(data[,i])
  # select the values with frequency more than 20
  tableLess <- table[table < 20]
  # replace the values with "Others"
  data[,i][data[,i] %in% names(tableLess)] <- "Others"
}
data

# save the data
write.table(data, file="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_Morocco_meeting_final/Data/passportData.order.representative.tsv", sep="\t", quote=FALSE, row.names=TRUE)
