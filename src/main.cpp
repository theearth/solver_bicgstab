#include "utils.h"
#include "mmio.h"
#include "bilu0.h"
#include "solver_mkl.h"
#include "solver_cusparse.h"
#include "dev_utils.h"

int main(int argc, char** argv)
{
	FILE *out;
	coo_t *ilu = NULL, *A = NULL;

	double *b = NULL;
	ilu = (coo_t*) malloc (sizeof(coo_t));
	A = (coo_t*) malloc (sizeof(coo_t));
	
	int device = 0; //0 - CPU, 1 - GPU
	int key_ilu = 0, key_A = 0, key_rhs = 0, dimb = 1, part = 1;
	
	if (argc == 1)
	{
		
		fprintf(stderr,"Usage: %s -A [martix-market-filename] -ilu [martix-market-filename] -rhs [martix-market-filename]\n",argv[0]);
		fprintf(stderr,"Optional: -d compute on GPU, -h compute on CPU(default), -b dim (where dim - number of phases on hydrodynamic model, default dim = 1)\n");
			exit(1);
	}
	else 
	    for (int i = 1; i < argc ; i++)
	    {
			if (!strcmp(argv[i],"-ilu"))
			{
			    key_ilu = 1;
			    if( read_coo_MM(argv[i+1], ilu) != 0)
					exit(1);
			}
			if (!strcmp(argv[i],"-A"))
			{
				key_A = 1;
				if( read_coo_MM(argv[i+1], A) != 0)
					exit(1);
			}
			if (!strcmp(argv[i],"-rhs"))
			{
				key_rhs = 1;
				if( read_coo_vector(argv[i+1], &b) != 0)
				exit (1);
			}
			if (!strcmp(argv[i],"-d"))
				device = 1;

			if (!strcmp(argv[i],"-b"))
				dimb = atoi(argv[i+1]);
		
			if (!strcmp(argv[i],"-p"))
			    part = atoi(argv[i+1]);
	    }
	
	if ( (key_A) && (key_rhs) )
	{
	    if (key_ilu == 0)
	    { 
			ilu = allocCOO (A->nnz);
			if( copyCOO (ilu,A) != 0) exit (1);
	    }
	    printf ("Load Matrix complete\n"); /*0-index based*/
	}
	else 
	{
	    printf ("Load Matrix Error\n");
	    exit(1);
	}

	
	csr_t* A_mat = allocCSR(A->nnz,A->n);
	csr_t* ilu_mat = allocCSR(ilu->nnz,ilu->n);
	
	coo2csr (A,A_mat);
	//sortCSR(A_mat);
	//printMatrix(A_mat);
	coo2csr (ilu,ilu_mat);
	//sortCSR(ilu_mat);
	freeCOO(A);
	freeCOO(ilu);
	printf ("Convert 2csr complete\n");
	//ilu0 and solve
	double *x = (double*)calloc(A_mat->n, sizeof(double));

	int key; //1 - block, 2 - usually
	int nb;

	if (device == 0)
	{
	    if (part > 1)
	    {
	        key = 1;
	        nb = part;
	    }
	    else 
		key = 2;
	}
	else if (part == 1)
	{
	    key = 3;
	}
	else 
	    {
			key = 4;
			nb = part;
	    }
	
	if ( key == 1 )
	{
		printf ("Usage block ilu0(%d)\n",nb);	
		/* Partition Matrix for ILU0 */
		double time_partition = omp_get_wtime();
		csr_t** bilu = partition(ilu_mat,nb,dimb);
		printf ("time partition %lf\n",omp_get_wtime()-time_partition);
		double time_bilu = omp_get_wtime();
		int err = bilu0(bilu,nb);
		printf ("time bilu0 %lf\n",omp_get_wtime()-time_bilu);
		csr_zero2one_index(A_mat);
		double time_solve = omp_get_wtime();
		err = solve_bicgstab_block(A_mat,bilu,nb, b, x);	
		printf ("Time Solve %lf\n",omp_get_wtime()-time_solve);
		for (int i = 0; i < nb; i++)
			freeCSR(bilu[i]);
	}
	else if ( key == 2 )
	{
		printf ("Usage ilu0\n");	
		double time_ilu = omp_get_wtime();
		int err = ilu0(ilu_mat);
		printf ("Time ILU0 %lf\n",omp_get_wtime()-time_ilu);
		csr_zero2one_index(A_mat);
		double time_solve = omp_get_wtime();
		err = solve_bicgstab(A_mat,ilu_mat, b, x);	
		printf ("Time Solve %lf\n",omp_get_wtime()-time_solve);
	}
	else if ( key == 3 )
	{
		printf ("Using block ilu0\n");
		init();
		int err = bicgstab_cusparse(A_mat,ilu_mat, b, x);
	}
	else if ( key == 4 )
	{
	    printf ("Using block ilu0(%d)\n",nb);	
	    /* Partition Matrix for ILU0 */
	    double time_partition = omp_get_wtime();
	    csr_t** bilu = partition(ilu_mat,nb,dimb);
	    csr_t** bA = partitionMatrix (A_mat,bilu,nb);
	    printf ("time partition %lf\n",omp_get_wtime()-time_partition);
	    //csr_zero2one_index(A_mat);
	    int err = bicgstab_fullblock_cusparse(bA, bilu, nb, b, x);
	    for (int i = 0; i < nb; i++)
	    {
			freeCSR(bilu[i]);
			freeCSR(bA[i]);
	    }
	}
	
	printVector(x,A_mat->n);
	
		
	freeCSR(A_mat);
	freeCSR(ilu_mat);
	free(A);
	free(ilu);
	free(A_mat);
	free(ilu_mat);
	free (b);
	free (x);
	return 0;
}