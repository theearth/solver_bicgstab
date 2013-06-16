#ifndef DEV_UTILS_H
#define DEV_UTILS_H

#include "utils.h"
#include "cuda_runtime.h"

#define CUDA_SAFE_CALL_NO_SYNC( call) {                               \
    cudaError err = call;                                             \
    if( cudaSuccess != err) {                                         \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                __FILE__, __LINE__, cudaGetErrorString( err) );       \
        exit(EXIT_FAILURE);                                           \
} }
#define CUDA_SAFE_CALL( call)    CUDA_SAFE_CALL_NO_SYNC(call);

void init();
void copymat2dev (csr_t* dev, csr_t* host);
void freemat (csr_t* dev);
csr_t* allocmat2dev (int nnz, int n);

#endif