#!/bin/bash

scriptUsage() {
    echo "This script is used to run vcftools to filter the vcf file"
    echo "Usage: ./vcf_filter.sh <inputData> <outData>"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi

# get the input parameters
inputData=$1
outData=$2
vcftools --vcf $inputData --maf 0.05 --max-missing 0.5 --recode --recode-INFO-all --out $outData  --min-alleles 2 

# because vcftools add .recode.vcf to the output file, we need to remove it
mv $outData.recode.vcf $outData