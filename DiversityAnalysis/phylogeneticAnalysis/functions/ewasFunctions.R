#!/usr/bin/env Rscript


# c function to read the geno file and split it into a matrix, where each line contains single line of a locus information for all individuals
readGeno<-cfunction(signature(filename="character", indvn="integer",locn="integer"), '
    FILE *file;
    char *fname;
    int indv = INTEGER(indvn)[0];
    int loc = INTEGER(locn)[0];
    fname = (char*)CHAR(STRING_ELT(filename,0));
    file = fopen(fname, "r");
    // create a matrix to store the geno file and fill it with zeros
    SEXP result;
    PROTECT(result = allocMatrix(INTSXP, indv, loc));
    int *r = INTEGER(result);
    memset(r, 0, indv*loc*sizeof(int));
    // read the file line by line
    char line[1000000];
    int i = 0;
    while(fgets(line, sizeof(line), file)) {
        // alocate a string to store the line
        char *str;
        str = (char*)malloc(sizeof(char)*strlen(line));
        memcpy(str, line, strlen(line));
        // loop over the line and split it into characters
        for(int j = 0; j < strlen(str); j++) {
            // if the character is 0, 1 or 2, store it in the matrix
            if(str[j] == \'0\' || str[j] == \'1\' || str[j] == \'2\' || str[j] == \'9\') {
                // calculate the position of the character in the matrix to fill the rows first
                int pos = (i % loc) * indv + (i / loc);
                r[pos] = str[j] - \'0\';
                i++;
            }
        }
        // free the memory
        free(str);
    }
    // free the memory
    fclose(file);
    UNPROTECT(1);
    return result; 
    ')


# fill missing with major allele 
# fillMissing a cfunction to fill the missing data with the major allele that recieves a matrix and returns a matrix
fillMissing<-cfunction(signature(x="matrix"), '
    int *r = INTEGER(x);
    int loc = nrows(x);
    int indv = ncols(x);
    // loop over the matrix and print the values
    for(int i = 0; i < loc; i++) {
        // loop over the individuals and calculate the number of 0, 1 and 2
        int n0 = 0;
        int n1 = 0;
        int n2 = 0;
        for(int j = 0; j < indv; j++) {
           int pos = i  + j* loc;
            if(r[pos] == 0) {
                n0++;
            } else if(r[pos] == 1) {
                n1++;
            } else if(r[pos] == 2) {
                n2++;
            } 
        }
        // calculate the major allele
        int major = 0;
        if(n1 > n0) {
            major = 1;
        }
        if(n2 > n1) {
            major = 2;
        }
        // loop over the individuals and replace the missing data with the major allele
        for(int j = 0; j < indv; j++) {
         int pos = i  + j* loc;
         if(r[pos] == 9) {
            r[pos] = major;
         }
        }
    }
    return x;
    ')


