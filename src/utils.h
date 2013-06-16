#ifndef UTILS_H
#define UTILS_H


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>

#include "mkl.h"
#include "mkl_rci.h"
#include <omp.h>
#include "typemat.h"


void freeCOO (coo_t* coo);
void freeCSR (csr_t* csr);
coo_t* allocCOO (int nnz);
csr_t* allocCSR (int nnz,int n);
int copyCOO (coo_t* dest,coo_t* src);
int coo2csr (coo_t* coo,csr_t* csr);
void csr_zero2one_index(csr_t* csr);
void csr_one2zero_index(csr_t* csr);
void printVector (double * vec, int n);
void printMatrix (csr_t* mat);
void sortCSR (csr_t* mat);

#endif