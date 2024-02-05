//############################################# About Author #########################################
// Created by: Alsamman M. Alsamman                                                                  #
// Emails: smahmoud [at] ageri.sci.eg or A.Alsamman [at] cgiar.org or SammanMohammed [at] gmail.com  #
// License: MIT License - https://opensource.org/licenses/MIT                                        #
// Disclaimer: The script comes with no warranty, use at your own risk                               #
// This script is not intended for commercial use                                                    #
//####################################################################################################

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *ModifyLine(char *line, int target_column)
{
    int col = 0;
    // split line by tab
    char *token = strtok(line, "\t");
    // loop through the line and print all the tokens
    while (token != NULL)
    {
        if (col == target_column)
        {
            printf(".");
        }
        else
        {
            printf("%s", token);
        }
        token = strtok(NULL, "\t");
        if (token != NULL)
        {
            printf("\t");
        }
        col++;
    }
    // // check if the number of columns is as expected
    // if (col != expected_columns)
    // {
    //     fprintf(stderr, "The number of columns is not as expected\n");
    //     fprintf(stderr, "Expected: %d and found: %d\n", expected_columns, col);
    //     exit(1);
    // }
}
int isContigLine(char *line)
{
    if (line[1] == '#' && line[2] == 'c' && line[3] == 'o' && line[4] == 'n' && line[5] == 't' && line[6] == 'i' && line[7] == 'g')
    {
        return 1;
    }
    return 0;
}
// only keep the first information in the INFO column in the vcf file
int main(int argc, char **argv)
{
    char *vcf_file = argv[1];
    int target_column = atoi(argv[2]);
    int remove_contig_info = 0;
    if (argc == 4)
    {
        remove_contig_info = atoi(argv[3]);
    }
    
    // int expected_columns = atoi(argv[3]);
    // for instance the INFO column is the 8th column
    FILE *file = fopen(vcf_file, "r");
    char line[1000000];
    while (fgets(line, 1024, file))
    {
        // if this line is a comment line, skip it
        if (line[0] == '#')
        {
            // if the line is a contig line, skip it
            if (remove_contig_info && isContigLine(line))
            {
                continue;
            }
            printf("%s", line);
            continue;
        }
        ModifyLine(line, target_column);
    }
    return 0;
}