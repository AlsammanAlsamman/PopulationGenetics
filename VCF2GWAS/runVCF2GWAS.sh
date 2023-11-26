#!/bin/bash

scriptUsage() {
    echo "This script is used to run vcf2gwas tool"
    echo "Usage: ./runVCF2GWAS.sh vcfFile phenData outFolder [usePCA]"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi


vcfFile=$1
phenData=$2
outFolder=$3

usePCA=$4

echo "####  Running vcf2gwas ####"
# if usePCA is empty
if [ -z "$usePCA" ]; then
    echo "Running vcf2gwas without PCA"
    vcf2gwas  -v $vcfFile -pf $phenData -ap -lmm -o $outFolder
else
    echo "Running vcf2gwas with PCA"
    vcf2gwas  -v $vcfFile -pf $phenData -ap -lmm -cf "PCA" -ac -o $outFolder
fi

echo "####  vcf2gwas finished ####"
echo "####  Copying results to gwas_out folder ####"
# copy all associated files to gwas_out folder
mkdir -p $outFolder/gwas_out
cp $outFolder/Output/*/*/*/*.assoc.txt $outFolder/gwas_out