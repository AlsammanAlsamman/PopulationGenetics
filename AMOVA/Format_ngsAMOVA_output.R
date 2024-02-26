#!/usr/bin/env Rscript
library(reshape2)


# ngsAMOVAoutFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/ngsAMOVA/amovaput.amova.csv"

# # use args
args<-commandArgs(trailingOnly=TRUE)
ngsAMOVAoutFile<-args[1]
# out folder is the same as the input file
resultsFolder<-dirname(ngsAMOVAoutFile)
# read the file
ngsAMOVAout<-read.csv(ngsAMOVAoutFile, sep=",", stringsAsFactors = FALSE)
colnames(ngsAMOVAout)<-c("calculation","level","value")
ngsAMOVAout


# select "SSD", "MSD", "df"
ngsAMOVAout_table<-ngsAMOVAout[ngsAMOVAout$calculation %in% c("df","MSD","SSD" ),]
# convert to table
ngsAMOVAout_table<-dcast(ngsAMOVAout_table, level~calculation, value.var="value")
#  make the Total df equal -
ngsAMOVAout_table$df[ngsAMOVAout_table$level=="Total"]<-"-"

# select Phi Variance_coefficient Variance_component Percentage_variance
ngsAMOVAout_variance<-ngsAMOVAout[ngsAMOVAout$calculation %in% c("Percentage_variance","Variance_component"),]
# convert to table
ngsAMOVAout_variance<-dcast(ngsAMOVAout_variance, level~calculation, value.var="value")
# add Total
ngsAMOVAout_variance<-rbind(ngsAMOVAout_variance, c("Total",sum(ngsAMOVAout_variance$Percentage_variance),"-"))
ngsAMOVAout_variance
# add percentage to the table
ngsAMOVAout_variance<-cbind(ngsAMOVAout_table, ngsAMOVAout_variance$Percentage_variance, ngsAMOVAout_variance$Variance_component)

outfilename<-ngsAMOVAout_variance$level[1]
# remove Among, within, Total, and -
outfilename<-gsub("Among","",outfilename)
outfilename<-gsub("within","",outfilename)
outfilename<-gsub("Total","",outfilename)
outfilename<-gsub("_","",outfilename)
outfilename
# rename the columns
colnames(ngsAMOVAout_variance)<-c("level","df","MSD","SSD","Percentage_variance","Variance_component")
# rename the rows to "Among Populations", "Within Population", "Total"
ngsAMOVAout_variance$level<-c("Among Populations", "Within Population", "Total")
ngsAMOVAout_variance
# write the file
write.csv(ngsAMOVAout_variance,paste0(resultsFolder,"/",outfilename,"_ngsAMOVA.csv"), row.names=FALSE)