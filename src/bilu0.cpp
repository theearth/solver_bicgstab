#include "bilu0.h"


csr_t** partition (const csr_t* mat,const int nb, const int dimb)
{
	int* ia = mat->ia;
	int* ja = mat->ja;
	double* val = mat->val;
	int n = mat->n;
	int nnz = mat->nnz;
	
	csr_t **bilu = (csr_t**) malloc(sizeof(csr_t)*nb);

if (nb>1)
{
	int *nnz_block = (int*) calloc(nb,sizeof(int));
	int* dim_block = (int*) calloc (nb,sizeof(int));
	int* offset_block = (int*) calloc (nb+1, sizeof(int));
	int dim_b;
	if (n % dimb == 0)
	    dim_b = dimb;
	else 
	    dim_b = 1;
	
	int cells = n / dim_b;
	int mod = cells % nb;
	
	for (int i = 0; i < nb; i++)
	{
		dim_block[i] = cells / nb;
		if (mod > 0)
		{
			dim_block[i]++;
			mod--;
		}
		dim_block[i] *= dim_b;
		offset_block[i+1] = offset_block[i] + dim_block[i];
	}
	
#pragma omp parallel for shared(nnz_block, dim_block, ia, ja) 
	for (int num_block = 0; num_block < nb; num_block++ )
		for (int i = offset_block[num_block]; i < offset_block[num_block+1]; i++)
			for (int j = ia[i]; j < ia[i+1]; j++)
			{
				if ( (ja[j]>=offset_block[num_block]) && (ja[j]<offset_block[num_block+1]) )
					nnz_block[num_block]++;
			}

	
	for (int num_block = 0; num_block < nb; num_block++ )
		bilu[num_block] = allocCSR(nnz_block[num_block], dim_block[num_block]);


#pragma omp parallel for shared(nnz_block, dim_block, ia, ja) 
	for (int num_block = 0; num_block < nb; num_block++ )
	{
		int num_ja = 0;
		bilu[num_block]->ia[0] = num_ja;
		for (int i = offset_block[num_block]; i < offset_block[num_block+1]; i++)
		{
			for (int j = ia[i]; j < ia[i+1]; j++)
			{
				if ( (ja[j]>=offset_block[num_block]) && (ja[j]<offset_block[num_block+1]) )
				{
					bilu[num_block]->ja[num_ja] = ja[j]-offset_block[num_block];
					bilu[num_block]->val[num_ja] = val[j];
					num_ja++;
				}
			}
			bilu[num_block]->ia[i+1-offset_block[num_block]] = num_ja;
		}
	}
	for (int i=0; i<nb;i++)
	    printf("ILU nnz[%d] = %d\n",i,nnz_block[i]);
	free(nnz_block);
	free(dim_block);
	free(offset_block);
}
else 
{
	bilu[0] = allocCSR(nnz, n);
	memcpy(bilu[0]->ia,mat->ia, sizeof(int) * (n+1));
	memcpy(bilu[0]->ja,mat->ja, sizeof(int) * (nnz));
	memcpy(bilu[0]->val,mat->val, sizeof(double) * (nnz));
}
	
	return bilu;
} 


csr_t** partitionMatrix (csr_t* mat, csr_t** bilu0, const int nb)
{

    csr_t **bmat = (csr_t**) malloc(sizeof(csr_t)*nb);
    int* offset_block = (int*) calloc (nb+1, sizeof(int));
    int* nnz_block = (int*) calloc (nb, sizeof(int));

    for (int i=0; i<nb; i++)
    {
        offset_block[i+1] = offset_block[i] + bilu0[i]->n;
        nnz_block[i] = mat->ia[offset_block[i+1]] - mat->ia[offset_block[i]];
	bmat[i] = allocCSR(nnz_block[i], bilu0[i]->n);
	bmat[i]->nnz = nnz_block[i];
	bmat[i]->n = bilu0[i]->n;
	
	printf ("Matrix A n[%d] = %d, nnz[%d] = %d\n",i,bmat[i]->n,i,bmat[i]->nnz);
	/*
	memcpy(bmat[i]->ia,&(mat->ia[offset_block[i]]), sizeof(int) * (bmat[i]->n+1));
	memcpy(bmat[i]->ja,&(mat->ja[mat->ia[offset_block[i]]]), sizeof(int) * (bmat[i]->nnz));
	memcpy(bmat[i]->val,&(mat->val[mat->ia[offset_block[i]]]), sizeof(double) * (bmat[i]->nnz));
	for (int j = 0; j <= bmat[i]->n; j++)
	    bmat[i]->ia[j] -= bmat[i]->ia[0];
	*/
	for (int j = 0; j <= bmat[i]->n; j++)
	{
	    bmat[i]->ia[j] = mat->ia[offset_block[i]+j] - mat->ia[offset_block[i]];
	}
	for (int j = 0; j < bmat[i]->nnz; j++)
	{
	    bmat[i]->ja[j] = mat->ja[mat->ia[offset_block[i]] + j];
	    bmat[i]->val[j] = mat->val[mat->ia[offset_block[i]] + j];
	}
    }

    free(nnz_block);
    free(offset_block);
    return bmat;
}

int bilu0 (csr_t** bilu0, const int nb)
{
	double dpar[128];
	int ipar[128],ierr, err=0;
	ipar[30] = 1;
	dpar[30] = 1.e-20;
	dpar[31] = 1.e-16;

#pragma omp parallel for private(ierr)
	for (int i=0;i<nb;i++)
	{
		csr_zero2one_index(bilu0[i]);
		dcsrilu0(&bilu0[i]->n,bilu0[i]->val,bilu0[i]->ia,bilu0[i]->ja,bilu0[i]->val,ipar,dpar,&ierr);
		if (ierr!=0)
		{
			err = ierr;
		}
	}
	return err;

}

int ilu0 (csr_t* ilu)
{
	double dpar[128];
	int ipar[128],ierr;
	ipar[30] = 1;
	dpar[30] = 1.e-20;
	dpar[31] = 1.e-16;
	double* val = (double*) malloc (sizeof(double) * ilu->nnz );
	csr_zero2one_index(ilu);
	dcsrilu0(&ilu->n,ilu->val,ilu->ia,ilu->ja,val,ipar,dpar,&ierr);
	free (ilu->val);

	ilu->val = val;
	

	return ierr;

}