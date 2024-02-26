#!/bin/bash

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


################# To print error
#TOOL INPUT > OUT 2>error
#################

# vcftools select samples
inputVCF="GOATs.vcf"
keep="ListWanted.txt"
output="GOATs"
vcftools --vcf $inputVCF --keep $keep --recode --recode-INFO-all --out $output

inputVCF="GOATs.vcf"
output="GOATs"
# filter using vcftools
vcftools --gzvcf $inputVCF --maf 0.05 --recode --recode-INFO-all --out $output.MAF