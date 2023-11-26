#!/bin/bash

scriptUsage() {
    echo "This script is used to run ****"
    echo "Usage: ./Command.sh *prameters*"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi


# run GEMMA kinship

echo "Running GEMMA Kinship"

genoypPrefix=$1
phenData=$2
ouPut=$3
outFolder=$4

# print a file with 1 repeated rows for each sample
# if the outputfolder ends with "/", remove it
if [[ $outFolder == */ ]]; then
    outFolder=${outFolder%?}
fi

mkdir -p $outFolder"/kinship"

# gemma -bfile $genoypPrefix -gk -o $ouPut\.gk -p $phenData -outdir $outFolder"/kinship"

echo "Kinship matrix is saved in $ouPut\.gk.cXX.txt"

# mkdir for LMM
mkdir -p $outFolder"/LMM"

echo "Running GEMMA LMM"

gemma -bfile $genoypPrefix \
    -p $phenData \
    -k $outFolder"/kinship/"$ouPut\.gk.cXX.txt \
    -lmm -o $ouPut\.lmm \
    -outdir $outFolder"/LMM"