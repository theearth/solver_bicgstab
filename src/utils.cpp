#include "utils.h"

void freeCOO (coo_t* mat)
{
	if (mat->ia) 
		free(mat->ia);
	if (mat->ja) 
		free(mat->ja);
	if (mat->val) 
		free(mat->val);
}

void freeCSR (csr_t* mat)
{
	if (mat->ia) 
		free(mat->ia);
	if (mat->ja) 
		free(mat->ja);
	if (mat->val) 
		free(mat->val);
}


coo_t* allocCOO (int nnz)
{
	coo_t* mat = (coo_t*) malloc (sizeof(coo_t));
	mat->ia = (int*) malloc(nnz * sizeof(int));
	mat->ja = (int*) malloc(nnz * sizeof(int));
	mat->val = (double*) malloc(nnz * sizeof(double));
	mat->nnz = nnz;
	return mat;
}

csr_t* allocCSR (int nnz, int n)
{
	csr_t* mat = (csr_t*) malloc (sizeof(csr_t));
	mat->ia = (int*) calloc(n+1, sizeof(int));
	mat->ja = (int*) malloc(nnz * sizeof(int));
	mat->val = (double*) malloc(nnz * sizeof(double));
	mat->nnz = nnz;
	mat->n = n;
	return mat;
}

int copyCOO (coo_t* mat_dest, coo_t* mat_src)
{
	if ((mat_src->ia)&&(mat_src->ja)&&(mat_src->val))
	{
		mat_dest->nnz = mat_src->nnz;
		mat_dest->n = mat_src->n;
		memcpy(mat_dest->ia,mat_src->ia,sizeof(int)*(mat_dest->nnz));
		memcpy(mat_dest->ja,mat_src->ja,sizeof(int)*(mat_dest->nnz));
		memcpy(mat_dest->val,mat_src->val,sizeof(double)*(mat_dest->nnz));
	}
	else 
		return 1;
	return 0;
}

int coo2csr (coo_t* coo,csr_t* csr)
{
	int i,j = 0,k = 0;
	
	memcpy(csr->ja,coo->ja,sizeof(int)*(csr->nnz));
	memcpy(csr->val,coo->val,sizeof(double)*(csr->nnz));
	
	for (i = 1; i<csr->n+1; i++)
	{
		csr->ia[i]+=csr->ia[i-1];
		for (j; j<coo->nnz;j++)
			if ( coo->ia[j] == i-1 )
				csr->ia[i]++;
			else break;
	}
	
	return 0;
}

void csr_zero2one_index(csr_t* csr)
{
	for(int i = 0; i < csr->n+1; i++)
		csr->ia[i]++;
	for (int i = 0; i < csr->nnz; i++)
		csr->ja[i]++;
}
void csr_one2zero_index(csr_t* csr)
{
	for(int i = 0; i < csr->n+1; i++)
		csr->ia[i]--;
	for (int i = 0; i < csr->nnz; i++)
		csr->ja[i]--;
}

void printVector (double * vec, int n)
{
	FILE * out = fopen("out_vector.txt", "w");
	for (int i = 0; i < n; i++)
		fprintf(out,"%.8le\n",vec[i]);
	fclose(out);
}

void printMatrix (csr_t* mat)
{
	FILE * out = fopen("out_mat.txt", "w");
	//csr_one2zero_index(mat);
	for (int j = 0; j < mat->n; j++)
		for (int k = mat->ia[j]; k < mat->ia[j+1]; k++)
			fprintf(out,"%d\t%d\t%.12le\n",j,mat->ja[k],mat->val[k]);
	fclose(out);
}

void sortCSR (csr_t* mat)
{
    for (int i = 0; i < mat->n; i++)
	for (int j = 0; j < (mat->ia[i+1] - mat->ia[i]); j++)
	    for (int k = j+1; k < (mat->ia[i+1] - mat->ia[i]); k++)
		if (mat->ja[mat->ia[i]+j] > mat->ja[mat->ia[i]+k])
	        {
	    	    int buf = mat->ja[mat->ia[i]+j];
	    	    mat->ja[mat->ia[i]+j] = mat->ja[mat->ia[i]+k];
	            mat->ja[mat->ia[i]+k] = buf;
	            double buf_val = mat->val[mat->ia[i]+j];
	    	    mat->val[mat->ia[i]+j] = mat->val[mat->ia[i]+k];
	            mat->val[mat->ia[i]+k] = buf_val;
	        }
}