#!/bin/bash

# converting vcf to plink

scriptUsage() {
    echo "This script is used to convert vcf to plink"
    echo "Usage: ./2Plink.sh vcfFile outputPrefix"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi

vcfFile=$1
outputPrefix=$2

plink --vcf $vcfFile -out $outputPrefix --allow-extra-chr
# You can use vcftools to convert VCF to plink format but it will not gave you the fam file
# vcftools --vcf $vcfFile --out $outputPrefix --plink