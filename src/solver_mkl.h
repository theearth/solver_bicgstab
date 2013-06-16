#ifndef SOLVER_MKL_H
#define SOLVER_MKL_H

#include "utils.h"

int solve_bicgstab(csr_t* mat, csr_t* ilu, double* b, double* x);

int solve_bicgstab_block(csr_t* mat, csr_t** ilu, int nb, double* b, double* x);

#endif