#!/bin/bash

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

 # calculate the LD between SNPs using PLINK


vcffile=Barley.vcf.recode.vcf
outfolder=plink

# convert vcf to plink
plink --vcf $vcffile --make-bed --out $vcffile --allow-extra-chr

plink --bfile $vcffile --r2 --ld-window-r2 0.2 --ld-window 99999 --ld-window-kb 99999 --out $outfolder/$vcffile --allow-extra-chr




insFolder="institutes"

vcffile="Barley.vcf.recode.vcf"

# loop over all files in the folder
for file in $insFolder/*.txt
do
    # from the vcffile name extract samples using vcftools
#    vcftools --vcf $vcffile --keep $file --recode --out $file
    # convert vcf to plink
     plink --vcf $file.recode.vcf --make-bed --out $file --allow-extra-chr
    # calculate LD
     plink --bfile $file --r2 --ld-window-r2 0.2 --ld-window 99999 --ld-window-kb 99999 --out $file\_\d --allow-extra-chr
done



#####################################################################################################


insFolder="instituteswheat"

vcffile="Wheat.vcf.recode.vcf"

# loop over all files in the folder
for file in $insFolder/*.txt
do
    # from the vcffile name extract samples using vcftools
    #vcftools --vcf $vcffile --keep $file --recode --out $file
    # convert vcf to plink
    #plink --vcf $file.recode.vcf --make-bed --out $file --allow-extra-chr
    # calculate LD
    #plink --bfile $file --r2 --ld-window-r2 0.2 --ld-window 99999 --ld-window-kb 99999 --out $file\_\d --allow-extra-chr
done


# compress folder using zip
zip -r LDs.zip LDs

################    
# selecting specific SNPs from the vcf file using vcftools

vcftools --vcf ISR003.txt.recode.vcf --snps IDsSNPs.txt --recode --out ISR003