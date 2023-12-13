#!/usr/bin/python

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

# Description: extract samples order from VCF file

import sys

# if no arguments provided print the help message
if len(sys.argv) == 1:
    print("Usage: getSamplesOrderFromVCF.py <inputFile> <outputFile>")
    sys.exit(1)

# check that we have 1 arguments
if len(sys.argv) != 3:
    print("Error: The script takes 1 arguments")
    sys.exit(1)



inputFile = sys.argv[1]
outputFile = sys.argv[2]

# loop through the file until we reach the line that starts with #CHROM
samplesOrder = []
with open(inputFile, 'r') as f:
    for line in f:
        if line.startswith("#CHROM"):
            line = line.strip()
            line = line.split("\t")
            # get the samples order
            samplesOrder = line[9:]
            break
# write the samples order to the output file
with open(outputFile, 'w') as f:
    f.write("\n".join(samplesOrder))