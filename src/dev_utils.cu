#include "dev_utils.h"

void init()
{
    cudaSetDevice(0);
}

void copymat2dev (csr_t* dev, csr_t* host)
{
    dev->n = host->n;
    dev->nnz = host->nnz;

    CUDA_SAFE_CALL (cudaMemcpy (dev->ia, host->ia, (dev->n+1)* sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL (cudaMemcpy (dev->ja, host->ja, dev->nnz * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL (cudaMemcpy (dev->val, host->val, dev->nnz * sizeof(double), cudaMemcpyHostToDevice));
}

void freemat (csr_t* dev)
{
    CUDA_SAFE_CALL (cudaFree (dev->ia));
    CUDA_SAFE_CALL (cudaFree (dev->ja));
    CUDA_SAFE_CALL (cudaFree (dev->val));
}

csr_t* allocmat2dev (int nnz, int n)
{
    csr_t* dev = (csr_t*) malloc (sizeof(csr_t));
    void *ia = NULL, *ja = NULL, *val = NULL;
    CUDA_SAFE_CALL (cudaMalloc ((void**)&ia, (n+1)*sizeof(int)));
    CUDA_SAFE_CALL (cudaMalloc ((void**)&ja, nnz * sizeof(int)));
    CUDA_SAFE_CALL (cudaMalloc ((void**)&val,nnz * sizeof(double)));
    dev->ia = (int*)ia;
    dev->ja = (int*)ja;
    dev->val = (double*)val;
    return dev;
}