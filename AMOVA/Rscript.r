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


# AMOVA_OUT_FILE = args[1]
AMOVA_OUT_FILE = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/important_SNPs/ALL_SNP_AMOVA.csv"
#read the file
ngsAMOVAout<-read.csv(AMOVA_OUT_FILE, sep=",", stringsAsFactors = FALSE)
colnames(ngsAMOVAout)
# plot the the Percentage_variance across SNP
library(ggplot2)
library(reshape2)

ngsAMOVAout$SNP<-as.character(ngsAMOVAout$SNP)
ngsAMOVAout$Percentage_variance<-as.numeric(ngsAMOVAout$Percentage_variance)

# select "Among Populations", "Within Population" for the Percentage_variance

percentage_variance<-ngsAMOVAout[ngsAMOVAout$level %in% c("Among Populations", "Within Population"),]
percentage_variance
# reshape the data
percentage_variance<-dcast(percentage_variance, SNP~level, value.var="Percentage_variance")
# reshape the data
percentage_variance<-melt(percentage_variance, id.vars="SNP")

# boxplot
percentageBoxPlot
percentageBoxPlot<-ggplot(percentage_variance, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot() +
  labs(title="AMOVA", x="level", y="Percentage Variance") +
  # theme(legend.position="bottom")+
  # add the mean on the boxplot
  stat_summary(fun.y=mean, geom="text", aes(label=round(..y..,2)), vjust=-1) 

# out folder is the same as the input file
resultsFolder<-dirname(AMOVA_OUT_FILE)
# save as pdf
pdf(paste0(resultsFolder,"/","AMOVA_percentage_variance.pdf"))
print(percentageBoxPlot)
dev.off()




# filter the data to select snps with high percentage_variance with the level "Among"
SNPs_AMOVA_Bars<-ggplot(ngsAMOVAout, aes(x=SNP, y=Percentage_variance, fill=level)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="AMOVA", x="SNP", y="Percentage variance") #+
  # theme(legend.position="bottom") +
  #  make the legend title "source of variance"
  #  guides(fill=guide_legend(title="source of variance"))

# save as pdf
pdf(paste0(resultsFolder,"/","AMOVA_percentage_variance_bars.pdf"), width=10, height=5)
print(SNPs_AMOVA_Bars)
dev.off()





