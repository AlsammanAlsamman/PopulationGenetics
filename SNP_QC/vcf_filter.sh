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

inputData=gbs_diversity_AGENT_wheat_ChineseSpringV2_1_biosampleIDs_MAF0o0001_ENA_woDup.vcf
outData=filtered/gbs_diversity_AGENT_wheat_ChineseSpringV2_1_biosampleIDs_MAF0o0001_ENA_woDup.vcf

vcftools --vcf $inputData --maf 0.01 --max-missing 0.2 --recode --recode-INFO-all --out $outData  --min-alleles 2 

# because vcftools add .recode.vcf to the output file, we need to remove it
mv $outData.recode.vcf $outData