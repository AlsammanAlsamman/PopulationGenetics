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

# inputVCF=$1
inputVCF="190117_ICARDA_GBS_771522_SNPs_2299_samples_fixed.vcf"
bgzip $inputVCF
tabix -p vcf $inputVCF".gz"


# # compress the vcf file using zip
# inputVCF="Barley_Selected_Genotypes.vcf"
# zip $inputVCF".zip" $inputVCF
