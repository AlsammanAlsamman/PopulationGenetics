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

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi

# vcfFile=$1

vcfFile="/home/samman/Documents/ICARDA/DrSawsan/Eman/Chickpea_MAGIC.vcf.gz"
fileList="/home/samman/Documents/ICARDA/DrSawsan/Eman/GenotypesList.txt"

# extract samples from vcf file using vcftools
vcftools --gzvcf $vcfFile --keep $fileList --recode --recode-INFO-all --out Chickpea_MAGIC



