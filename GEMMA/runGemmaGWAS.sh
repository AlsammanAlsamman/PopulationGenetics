#!/bin/bash

scriptUsage() {
    echo "This script is used to run GWAS using GEMMA"
    echo "Usage: ./runGemmaGWAS.sh genoypPrefixPlink phenData kinship output"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi

genoypPrefix=$1
phenData=$2
kinship=$3
output=$4

gemma -bfile $genoypPrefix \
    -p $phenData \
    -k $kinship -lmm -o $output