#!/bin/bash

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


scriptUsage() {
    echo "This script is used to run ****"
    echo "Usage: ./Command.sh *prameters*"
}



# About the script
# Match the vcf file with the phenotypic data

phenotypicData=$1
vcfFile=$2
outFolder=$3
sep=$4

# get the sample names from the phenotypic data
if [[ -z $sep ]];
then
    sep="\t"
fi

# get the sample names from the phenotypic data
cut -f1 $phenotypicData | tail -n +2 > $outFolder/${phenotypicData##*/}.samples

# use vcftools to select the samples from the vcf file
vcftools --vcf $vcfFile --keep $outFolder/${phenotypicData##*/}.samples --recode --out $outFolder/${vcfFile##*/}.samples







