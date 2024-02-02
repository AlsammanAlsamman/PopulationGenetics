#!/usr/bin/env Rscript

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

##install the GWLD packages

#devtools::install_github("Rong-Zh/GWLD/GWLD-R")

##using packages

library(GWLD)

#0) Loading example data with base-type genotype

data(duck)


data <- duck$SNP
names(data)
SNP <- data$genotype
head(SNP)
Info <- data$info
head(Info)
## Read data from vcf or plink format file

# Recode the genotypes in vcf file as 0, 1, 2, NA, or recode to another type

#vcf_data <- read.vcf("vcf format file", genotype="int")

#plink_data <- read.plink("plink format file’s prefix")

#1) recode the genotype with 0,1,2,NA(NA for missing values)

SNP <- codegeno(SNP, sep="/")
SNP[1:10,1:10]
#2) calculate with different methods (D, D’, r2, RMI, MI)

result <- GWLD(SNP, method = "r^2", cores = 1)
head(result)
##or use the following code `

r2 <- LD(SNP, method = "r^2", cores = 1)
head(r2)
rmi <- RMI(SNP, cores=1)


mi <- MI(SNP, cores=1)

#2.1) LD heatmap plot
SNP[1:10,1:10]
Info$POS
Info$ID
class(SNP)
p <- HeatMap(SNP, method = "RMI", SnpPosition = Info$POS, SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)
p
#3) decay from result(step 2)

rmi_decay <- decay(rmi, Info)

#or calculate circos from data

rmi_decay <- calc_decay(SNP, Info, method="RMI")
# convert Dist and Value to numeric
rmi_decay$Dist <- as.numeric(rmi_decay$Dist)
rmi_decay$Value <- as.numeric(rmi_decay$Value)
#4) circos from result(step 2)

rmi_circos <- circos(rmi, Info, threshold=0.2)

#or calculate circos from data

#rmi_circos <- calc_circos(SNP, Info, method="RMI", threshold=0.2)

#4.1) circos plot

circosdata <- duck$Circos
circosdata
circos.ideogram(circosdata$chr)

circos.linksnp(circosdata$linkdata, bg.col=rainbow(29))


############### using vcf

library(GWLD)

# Recode the genotypes in vcf file as 0, 1, 2, NA, or recode to another type
vcf_data <- read.vcf("/home/samman/Documents/publishing/Drought_Ahmed_Master/ReDoingGWAS/Data/GenotypicData.vcf", genotype="int")
head(vcf_data)

# remove multi-allelic SNPs that have , in the ALT column
vcf_data <- vcf_data[vcf_data[, 5] != ",", ]
# remove SNPs with . in the ALT column
vcf_data <- vcf_data[vcf_data[, 5] != ".", ]


# select only SNPs located on CA1 
vcf_data <- vcf_data[vcf_data[, 1] == "CA1", ]
# convert to data frame
vcf_data <- as.data.frame(vcf_data)
vcf_data$POS <- as.numeric(vcf_data$POS)
# sort the SNPs based on their position increasing
vcf_data <- vcf_data[order(vcf_data[, 2], decreasing = FALSE), ] 
vcf_data[1:10,1:10]
# SNPs are after the 9th column
SNPs <- vcf_data[, 10:ncol(vcf_data)]
Info <- vcf_data[, 1:9]

colnames(Info) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
Info<-as.data.frame(Info)
Info$ID <- paste0(Info$CHROM,":",Info$POS,"_",Info$REF,",",Info$ALT)

# extract only   CHROM    POS            ID
Info<-Info[,c("CHROM","POS","ID")]

# convert to numeric by applying as.numeric to each column
SNPs <- apply(SNPs, 2, as.numeric)
SNPs <-t(SNPs)
colnames(SNPs) <- Info$ID

# # Select Random 50 SNPs
# snprandom <- sample(1:nrow(SNPs), 50)
# SNPs<-SNPs[,snprandom]
# Info<-Info[snprandom,]

Info
#2) calculate with different methods (D, D’, r2, RMI, MI)
result <- GWLD(SNPs, method = "r^2", cores = 1)

##or use the following code `

r2 <- LD(SNPs, method = "r^2", cores = 1)

rmi <- RMI(SNPs, cores=1)
rmi

mi <- MI(SNPs, cores=1)
# rotate the INFO data frame upside down
Info <- Info[nrow(Info):1, ]
#2.1) LD heatmap plot
p <- HeatMap(SNPs, method = "RMI", SnpPosition = as.numeric(Info$POS), SnpName = Info$ID, cores = 1, color = "YellowTored", showLDvalues = F)
p
















#3) decay from result(step 2)

#rmi_decay <- decay(rmi, Info)

#or calculate circos from data

#rmi_decay <- calc_decay(SNP, Info, method="RMI")

#4) circos from result(step 2)

#rmi_circos <- circos(rmi, Info, threshold=0.2)

#or calculate circos from data

#rmi_circos <- calc_circos(SNPs, Info, method="RMI", threshold=0.2)

#4.1) circos plot

#circosdata <- duck$Circos
#duck$Circos
#circos.ideogram(circosdata$chr)
#circos.linksnp(circosdata$linkdata, bg.col=rainbow(29))









####################

