#!/bin/bash

scriptUsage() {
    echo "This script is used to calculate the kinship matrix using GEMMA"
    echo "Usage: ./runGemmaKinship.sh genoypPrefixPlink phenData ouPutFolder"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi


genoypPrefix=$1
phenData=$2
ouPut=$3

gemma -bfile $genoypPrefix -gk -o $ouPut -p $phenData

echo "Done!"