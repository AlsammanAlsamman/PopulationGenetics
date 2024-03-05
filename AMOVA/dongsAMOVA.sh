#!/bin/bash

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


popsFile="pops.txt"
metadataFile="snpMeta.tsv"
vcfFile="BreadWheatBefore_Filtered.vcf"

lineNumber=2
# read all lines in the file
while IFS= read -r line 
do
  # display $line or do somthing with $line
  echo "$line"
  formula="Individual~$line"

  # split columns by tab in metadata file and take the first column and the column of the line
  cut -f 1,$lineNumber $metadataFile > temp_metadata.tsv
  echo $lineNumber
  # increment the column number
  lineNumber=$((lineNumber+1))
  # change the name of metadata file
  metadataFile="temp_metadata.tsv"
  # do ngsAMOVA
  ./ngsAMOVA --in-vcf $vcfFile -doEM 1 -doAMOVA 1 --printDistanceMatrix 1 --minInd 2 -doDist 2 --maxEmIter 100 --em-tole "1e-10" -m ${metadataFile} -f ${formula}
  # format the output
  Rscript Format_ngsAMOVA_output.R amovaput.amova.csv
  metadataFile="snpMeta.tsv"

done <"$popsFile"












