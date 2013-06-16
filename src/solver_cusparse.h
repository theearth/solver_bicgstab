#ifndef SOLVER_CUSPARSE_H
#define SOLVER_CUSPARSE_H

#include "utils.h"
#include "dev_utils.h"
#include <cusparse_v2.h>
#include <cublas_v2.h>


int bicgstab_cusparse(csr_t* mat, csr_t* ilu, double* b, double* x);
int bicgstab_block_cusparse(csr_t* mat, csr_t** ilu, int nb, double* b, double* x);
int bicgstab_fullblock_cusparse(csr_t** mat, csr_t** ilu, int nb, double* b, double* x);

#endif