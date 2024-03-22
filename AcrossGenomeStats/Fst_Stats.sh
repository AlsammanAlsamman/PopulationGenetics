#!/bin/bash

# Calculate Fst for each population pair using vcftools

# Define the input and output directories

inputMeta="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/passportData.order.representative.tsv"
metaDir="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/metaDir"
vcfFile="/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/BreadWheatBefore_Filtered.vcf"

# loop through all folders in the metaDir
for folder in $metaDir/*; do
    # get the population name
    col=$(basename $folder)
    metafile=$folder/meta.txt
    # loop through meta.txt file in the metaDir
    while IFS=$'\t' read -r -a line; do
        # get the accession name
        ga=${line[0]}
        # get the accession name
        gb=${line[1]}
        echo "Calculating Fst for $ga and $gb in $col"
        gaSampleFile=$folder/"$ga.txt"
        gbSampleFile=$folder/"$gb.txt"
        cat $gaSampleFile
    done < $metafile
done