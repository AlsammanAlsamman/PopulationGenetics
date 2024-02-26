#!/usr/bin/env Rscript


############################### Functions ##########################################
# calculate PIC in a vector of frequencies in C
PICvec<-cfunction(c(freq="numeric"),"
    // get the length of the vector
    int n = LENGTH(freq);
    // create a 
    double *p = REAL(freq);
    // create a matrix of outer products
    double o[n][n];
    double PIC;

    // Calculate the outer product of p^2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            o[i][j] = p[i] * p[i] * p[j] * p[j];
        }
    }

    // Calculate PIC
    double sum_p_squared = 0.0;
    double sum_upper_tri_o = 0.0;

    for (int i = 0; i < n; i++) {
        sum_p_squared += p[i] * p[i];
        for (int j = i + 1; j < n; j++) {
            sum_upper_tri_o += 2 * o[i][j];
        }
    }
  PIC = 1 - sum_p_squared - sum_upper_tri_o;
  return ScalarReal(PIC);
")

