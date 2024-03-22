#!/usr/bin/python

############################################## About Author #########################################
# Created by: Alsamman M. Alsamman                                                                  #
# Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
# License: MIT License - https://opensource.org/licenses/MIT                                        #
# Disclaimer: The script comes with no warranty, use at your own risk                               #
# This script is not intended for commercial use                                                    #
#####################################################################################################

import sys
import os
import pandas as pd

# # if no arguments provided print the help message
# if len(sys.argv) == 1:
#     print("Usage:")
#     print("Example: ")
#     sys.exit(1)

# working_dir = sys.argv[1]
# metfile = sys.argv[2]

working_dir = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/metaDir"
metfile = "/home/samman/Documents/ICARDA_Zakaria/Tasks/Fatima_final_analysis_Feb_24_2023/Data/passportData.order.representative.tsv"

# read the input file
with open(metfile) as f:
    lines = f.readlines()

# split the input file into two files
colnames = lines[0].strip().split("\t")
# create a folder for each column
for col in colnames[1:]:
    # check if the folder exists
    if not os.path.exists(working_dir + "/" + col):
       os.makedirs(working_dir + "/" + col)

# read the input file
df = pd.read_csv(metfile, sep="\t", header=0, index_col=0)
# for each column select the category and save it in a separate file
for col in colnames[1:]:
    # get values of the column
    values = df[col].unique()
    # for each value save the corresponding rows
    for val in values:
        # select the rows
        df1 = df[df[col] == val]
        # save the rows names as a text file
        df1.index.to_series().to_csv(working_dir + "/" + col + "/" + val + ".txt", header=False, index=False)

# create a file contains teh combination of the columns

# get the combinations of the columns

for col in colnames[1:]:

    # create a meta file in the column folder
    metacol = working_dir + "/" + col + "/meta.txt"
    metafileobj = open(metacol, "w")
    # get values of the column
    values = df[col].unique()
    # for each value save the corresponding rows
    for i in range(len(values)):
        for j in range(i+1, len(values)):
            g_a = values[i]
            g_b = values[j]
            # print(g_a, g_b, sep="\t")
            # write in the column folder a meta file
            metafileobj.write(g_a + "\t" + g_b + "\n")
    metafileobj.close()