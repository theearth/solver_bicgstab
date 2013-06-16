#ifndef BILU0_H
#define BILU0_H

#include "utils.h"

csr_t** partition (const csr_t* mat,const int nb, const int dimb);
int bilu0 (csr_t** bilu0, const int nb);
int ilu0 (csr_t* ilu);
csr_t** partitionMatrix (csr_t* mat, csr_t** bilu0, const int nb);

#endif