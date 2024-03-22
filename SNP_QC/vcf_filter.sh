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
maxMissing=0.8
minAlleles=2
maf=0.01
compreesed=0

# maxMissing=0.8
# minAlleles=2
# maf=0.01

# inputData="ICARDA_GBS_771522_SNPs_2299.gz"
# outData="ICARDA_GBS_771522_SNPs_2299_filtered"
# maxMissing=0.9
# minAlleles=2
# maf=0.05
# compreesed=1

vcftools $(if [ $compreesed -eq 1 ]; then echo "--gzvcf"; else echo "--vcf"; fi)\
    $inputData --maf $maf\
    --max-missing $maxMissing \
    --min-alleles $minAlleles \
    --recode --recode-INFO-all --out $outData

# use bcftools to ge stastistics
bcftools stats $outData.recode.vcf > $outData.stats


# because vcftools add .recode.vcf to the output file, we need to remove it
# mv $outData.recode.vcf $outData

#########################Info########################
# --max-missing <float>
# Exclude sites on the basis of the proportion of missing 
#data (defined to be between 0 and 1, where 0 allows sites 
#that are completely missing and 1 indicates no missing data allowed).