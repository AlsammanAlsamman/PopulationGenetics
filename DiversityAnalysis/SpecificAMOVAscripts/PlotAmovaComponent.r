#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

TableCompFile<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS/AMOVA_SNPs/AMOVA_RESULTS/amovaComponents.csv"

CompVar<-read.table(TableCompFile, header=TRUE, sep=",")
colnames(CompVar)<-c("Sigma","Percentage","SNP","rowname","varianceCategory")
# plot
library(ggplot2)
library(ggpubr)

# remove Total
CompVar<-CompVar[CompVar$varianceCategory!="Total",]
# convert all negative values to 0
CompVar[CompVar<0]<-0

# Select the "Between samples within populations"
BSP<-CompVar[CompVar$varianceCategory=="Between samples within populations",]
# select SNP and Percentage
BSP<-BSP[,c("SNP","Percentage")]
# sort samples by percentage
BSP<-BSP[order(BSP$Percentage,decreasing = FALSE),]

# plot the variance components
p<-ggplot(CompVar,aes(x=SNP,y=Percentage,fill=varianceCategory))+
  geom_bar(stat="identity",position = "fill")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x="SNP",y="Percentage",fill="Variance Category")+
  ggtitle("Variance Components")+
  # make the x axis labels as the SNP names smaller
  theme(axis.text.x = element_text(size=6))+
  # make the x axis order as the order of the samples in the BSP
  scale_x_discrete(limits=BSP$SNP)
p
path<-"/home/samman/Documents/ICARDA_Zakaria/Tasks/Nassima_5_23_12_2023/RESULTS/AMOVA_SNPs/"
ggsave(paste0(path,"VarianceComponents.png"),p,width=20,height=5,dpi=600)


