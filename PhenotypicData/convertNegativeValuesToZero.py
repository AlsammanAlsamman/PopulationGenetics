#!/usr/bin/python

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

import sys

# About
# The script changes all negative values to zero

# if no arguments provided print the help message
if len(sys.argv) == 1:
    print("Usage: python convertNegativeValuesToZero.py <input_file> <output_file> <splitter>")
    print("Example: python convertNegativeValuesToZero.py input.txt output.txt \",\"")
    sys.exit(1)
    
# check that we have 1 arguments
if len(sys.argv) < 2:
    print("Error: The script takes at least one argument")
    sys.exit(1)

# read the file name from the command line
filename = sys.argv[1]
# infile extension
extension = filename.split(".")[-1]
ofile = filename + ".converted" + "." + extension
splitter = "\t"
if len(sys.argv) > 4:
    ofile = sys.argv[2]
    splitter = sys.argv[3]

# open the file
linenumber = -1
# open output file
sys.stdout = open(ofile, 'w')
with open(filename) as f:
    for line in f:
        linenumber += 1
        if linenumber == 0:
            print(line.strip())
            continue
        linedata = line.strip().split(splitter)
        for i in range(1, len(linedata)):
            if float(linedata[i]) < 0:
                linedata[i] = "0"
        print(splitter.join(linedata))
sys.stdout.close()
