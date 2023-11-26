#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


scriptDescription = "This script converts categorical values to numerical values. 
The script takes 4 arguments: input file, column name, replace dictionary, and output file.
The replace dictionary is a comma separated list of key:value pairs.
The script will replace the key with the value in the column. 
The script will exit with an error if there are still categorical values in the column after the replacement."
scriptUsage = "Rscript ConvertCategorical2Numerical.r input_file column_name replace_dictionary output_file [separator]"


args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  # print the script description
  cat(scriptDescription, "\n")
  # print the script usage
  cat(scriptUsage, "\n")
  # exit R
  q()
}


infile = args[1]
colname = args[2]
replaceDict = args[3]
outfile = args[4]
sep = ""

if (length(args) > 4) {
  sep = args[5]
} else {
  sep = ","
}


# infile = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/AdultsTraits_Ordered.csv"
# colname = "STB.Reaction"
# replaceDict = "HR:7,R:6,MR:5,I:4,MS:3,S:2,HS:1"
# outfile = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_3_25_10_2023/Data/AdultsTraits_Ordered2.csv"
# sep=","

# read the input file

df = read.csv(infile, header = TRUE, sep = sep, stringsAsFactors = FALSE)

# convert the replaceDict Table to a dictionary
replaceDict = strsplit(replaceDict, ",")[[1]]
replaceDict = strsplit(replaceDict, ":")

# replace the values
for (i in 1:length(replaceDict)) {
  df[df[,colname] == replaceDict[[i]][1], colname] = replaceDict[[i]][2]
}

# is there any categorical value left?
if (length(unique(df[,colname])) > 1) {
  stop("There are still categorical values in the column ", colname, ".n", call.=FALSE)
  # exit R
  q()
}

# save the output file
write.table(df, file = outfile, sep = sep, row.names = FALSE, quote = FALSE)