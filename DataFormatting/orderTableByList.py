#!/usr/bin/python

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################


# ordering table by list of values in a column

########## import modules
import sys

########## Checks

# if no arguments provided print the help message
if len(sys.argv) == 1:
    print("Usage: orderTableByList.py <inputFileTarget> <inputFileList> <outputFile>")
    print("Example: orderTableByList.py Target.txt List.txt Output.txt")
    sys.exit(1)
    
# check that we have 3 arguments
if len(sys.argv) != 4:
    print("Error: The script takes 3 arguments")
    sys.exit(1)

########## Main
inputFileTarget = sys.argv[1]
inputFileList = sys.argv[2]
outputFile = sys.argv[3]

# check if the input file is tab delimited
with open(inputFileTarget, 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split("\t")
        if len(line) > 1:
            break
        else:
            print("Error: The input file is not tab delimited")
            sys.exit(1)
#

# read the list
list = []
with open(inputFileList, 'r') as f:
    for line in f:
        list.append(line.strip())

# read the target file
header = ""
TargetFile = {}
lineCounter = 0

with open(inputFileTarget, 'r') as f:
    for line in f:
        lineCounter += 1
        if lineCounter == 1:
            header = line.strip()
            continue
        line = line.strip()
        line = line.split("\t")
        TargetFile[line[0]] = line[1:]

# write the output file
with open(outputFile, 'w') as f:
    f.write(header + "\n")
    for item in list:
        f.write(item + "\t" + "\t".join(TargetFile[item]) + "\n")