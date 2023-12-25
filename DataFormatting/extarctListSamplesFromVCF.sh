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
    echo "Usage: ./extractListSamplesFromVCF.sh <vcfFile> <fileList> <outFile>"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi

vcfFile=$1
fileList=$2
outFile=$3

# extract samples from vcf file using vcftools
# if the is compressed, use --gzvcf
if [[ $vcfFile == *.gz ]]; then
    vcftools --gzvcf $vcfFile --keep $fileList --recode --recode-INFO-all --out $outFile
else
    vcftools --vcf $vcfFile --keep $fileList --recode --recode-INFO-all --out $outFile
fi

# get the new vcf file sample names
newVcfFile=$(ls $outFile".recode.vcf")
newVcfFileSampleNames=$(grep "#CHROM" $newVcfFile | cut -f10-)
echo $newVcfFileSampleNames | tr ' ' '\n' > $outFile".samples"


