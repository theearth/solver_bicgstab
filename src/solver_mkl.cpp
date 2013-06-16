#include "solver_mkl.h"
#define TIMER

int solve_bicgstab_block(csr_t* mat, csr_t** ilu, int nb, double* b, double* x)
{
	int n = mat->n;
	int nnz = mat->nnz;
	int *offset_ilu = (int*) calloc(nb, sizeof(int));
	for ( int i = 1; i < nb; i++ )
		offset_ilu[i] = offset_ilu[i-1] + ilu[i-1]->n;
	
	double tol = 1e-6, floatone = 1.0;
	const int max_iter = 200;
    
	double *r, *p, *y, *zm1, *zm2, *rm2, *rm1, *rm3, nrm0, nrm;
	r = (double*) malloc (sizeof(double) * n);
	p = (double*) malloc (sizeof(double) * n);
	y = (double*) malloc (sizeof(double) * n);
	rm1 = (double*) malloc (sizeof(double) * n);
	rm2 = (double*) malloc (sizeof(double) * n);
	rm3 = (double*) malloc (sizeof(double) * n);
	zm1 = (double*) malloc (sizeof(double) * n);
	zm2 = (double*) malloc (sizeof(double) * n);

	double rho = 1.0, rho1, beta = 0.0, alpha = 0.0, omega, temp, temp1;

	char lower1 = 'L', lower = 'N', lower2 = 'U';
	char upper1 = 'U', upper = 'N', upper2 = 'N';
	
	#ifdef TIMER
	double timerLUSol = 0, timerLUSol1, timerSpMV = 0, timerSpMV1;
	double timerTotal = omp_get_wtime();
	#endif

	cblas_dcopy (n, b, 1, r, 1);
	cblas_dcopy (n, r, 1, p, 1);
	cblas_dcopy (n, r, 1, zm1, 1);
	nrm0 = cblas_dnrm2 (n, r, 1);

	for (int k = 0; k < max_iter; k++)
	{
		rho1 = rho;
		rho = cblas_ddot(n, zm1, 1, r, 1);
		if ( k > 0 )
		{
			beta = (rho / rho1) * (alpha / omega);
			cblas_daxpy (n, -omega, zm2, 1, p, 1);
			cblas_dscal (n, beta, p, 1);
			cblas_daxpy (n, floatone, r, 1, p, 1);
		}
		
		#ifdef TIMER
		timerLUSol1 = omp_get_wtime();
		#endif
		
		#pragma omp parallel
		{
		#pragma omp for 
		for (int i = 0; i < nb; i++)
			mkl_dcsrtrsv (&lower1, &lower, &lower2, &ilu[i]->n, ilu[i]->val, ilu[i]->ia, ilu[i]->ja, &p[offset_ilu[i]], &y[offset_ilu[i]]);

		#pragma omp for 
		for (int i = 0; i < nb; i++)
			mkl_dcsrtrsv (&upper1, &upper, &upper2, &ilu[i]->n, ilu[i]->val, ilu[i]->ia, ilu[i]->ja, &y[offset_ilu[i]], &rm2[offset_ilu[i]]);
		}
		#ifdef TIMER
		timerLUSol += omp_get_wtime() - timerLUSol1;
		timerSpMV1 = omp_get_wtime();
		#endif
		
		mkl_dcsrgemv (&lower, &n, mat->val, mat->ia, mat->ja, rm2, zm2);
		
		#ifdef TIMER
		timerSpMV += omp_get_wtime() - timerSpMV1;
		#endif
		
		temp = cblas_ddot(n, zm1, 1, zm2, 1);
		
		alpha = rho / temp;
		cblas_daxpy (n, -alpha, zm2, 1, r, 1);
		cblas_daxpy (n, alpha, rm2, 1, x, 1);
		nrm = cblas_dnrm2 (n, x, 1);

		if  ((nrm < tol) && ( nrm / nrm0 < tol )) {printf("  iteration = %3d, residual = %le \n", k+1, nrm / nrm0); break; }
		
		#ifdef TIMER
		timerLUSol1 = omp_get_wtime();
		#endif
		
		#pragma omp parallel
		{
		#pragma omp for
		for (int i = 0; i < nb; i++)
			mkl_dcsrtrsv (&lower1, &lower, &lower2, &ilu[i]->n, ilu[i]->val, ilu[i]->ia, ilu[i]->ja, &r[offset_ilu[i]], &y[offset_ilu[i]]);
		#pragma omp for 
		for (int i = 0; i < nb; i++)
			mkl_dcsrtrsv (&upper1, &upper, &upper2, &ilu[i]->n, ilu[i]->val, ilu[i]->ia, ilu[i]->ja, &y[offset_ilu[i]], &rm3[offset_ilu[i]]);
		}
		#ifdef TIMER
		timerLUSol += omp_get_wtime() - timerLUSol1;
		timerSpMV1 = omp_get_wtime();
		#endif
		mkl_dcsrgemv (&lower, &n, mat->val, mat->ia, mat->ja, rm3, y);
		
		#ifdef TIMER
		timerSpMV += omp_get_wtime() - timerSpMV1;
		#endif
		
		temp = cblas_ddot(n, y, 1, r, 1);
		temp1 = cblas_ddot(n, y, 1, y, 1);
		omega = temp / temp1;

		cblas_daxpy (n, omega, rm3, 1, x, 1);
		cblas_daxpy (n, -omega, y, 1, r, 1);
		nrm = cblas_dnrm2 (n, r, 1);
		if ((nrm < tol) && ( nrm / nrm0 < tol )) {printf("  iteration = %3d, residual = %le \n", k+1, nrm / nrm0); break; }
		printf("  iteration = %3d, residual = %le \n", k+1, nrm / nrm0);
	}
	#ifdef TIMER
	printf("time LUSol\t%lf\ntime SpMV\t%lf\n",timerLUSol,timerSpMV);
	printf("time total\t%lf\n",omp_get_wtime()-timerTotal);
	#endif
	
	free (r);
	free (rm1);
	free (rm2);
	free (rm3);
	free (zm1);
	free (zm2);
	free (p);
	free (y);
	free (offset_ilu);

	return 0;
}

int solve_bicgstab(csr_t* mat, csr_t* ilu, double* b, double* x)
{
	double tol = 1e-6, floatone = 1.0;
	const int max_iter = 200;
	int n = mat->n;
	int nnz = mat->nnz;
    
	double *r, *p, *y, *zm1, *zm2, *rm2, *rm1, *rm3, nrm0, nrm;
	r = (double*) calloc (n, sizeof(double));
	p = (double*) calloc (n, sizeof(double));
	y = (double*) calloc (n, sizeof(double));
	rm1 = (double*) calloc (n, sizeof(double));
	rm2 = (double*) calloc (n, sizeof(double));
	rm3 = (double*) calloc (n, sizeof(double));
	zm1 = (double*) calloc (n, sizeof(double));
	zm2 = (double*) calloc (n, sizeof(double));
	

	double rho = 1.0, rho1, beta = 0.0, alpha = 0.0, omega, temp, temp1;

	char lower1 = 'L', lower = 'N', lower2 = 'U';
	char upper1 = 'U', upper = 'N', upper2 = 'N';
	
	#ifdef TIMER
	double timerLUSol = 0, timerLUSol1, timerSpMV = 0, timerSpMV1, timerVector = 0, timerVector1;
	double timerTotal = omp_get_wtime();
	#endif
	
	cblas_dcopy (n, b, 1, r, 1);
	cblas_dcopy (n, r, 1, p, 1);
	cblas_dcopy (n, r, 1, zm1, 1);

	nrm0 = cblas_dnrm2 (n, r, 1);
	for (int k = 0; k < max_iter; k++)
	{
		rho1 = rho;
		#ifdef TIMER
		timerVector1 = omp_get_wtime();
		#endif
		rho = cblas_ddot(n, zm1, 1, r, 1);
		if ( k > 0 )
		{
			
			beta = (rho / rho1) * (alpha / omega);
			cblas_daxpy (n, -omega, zm2, 1, p, 1);
			cblas_dscal (n, beta, p, 1);
			cblas_daxpy (n, floatone, r, 1, p, 1);
		}
		
		#ifdef TIMER
		timerVector += omp_get_wtime() - timerVector1;
		timerLUSol1 = omp_get_wtime();
		#endif
		
		mkl_dcsrtrsv (&lower1, &lower, &lower2, &n, ilu->val, ilu->ia, ilu->ja, p, y);
		mkl_dcsrtrsv (&upper1, &upper, &upper2, &n, ilu->val, ilu->ia, ilu->ja, y, rm2);
		
		#ifdef TIMER
		timerLUSol += omp_get_wtime() - timerLUSol1;
		timerSpMV1 = omp_get_wtime();
		#endif
		
		mkl_dcsrgemv (&lower, &n, mat->val, mat->ia, mat->ja, rm2, zm2);
		
		#ifdef TIMER
		timerSpMV += omp_get_wtime() - timerSpMV1;
		timerVector1 = omp_get_wtime();
		#endif
		
		temp = cblas_ddot(n, zm1, 1, zm2, 1);
		alpha = rho / temp;
		cblas_daxpy (n, -alpha, zm2, 1, r, 1);
		cblas_daxpy (n, alpha, rm2, 1, x, 1);

		nrm = cblas_dnrm2 (n, r, 1);
		#ifdef TIMER
		timerVector += omp_get_wtime() - timerVector1;
		#endif
		if ((nrm < tol) && ( nrm / nrm0 < tol )) {printf("  iteration = %3d, residual = %le \n", k+1, nrm / nrm0); break; }
		
		#ifdef TIMER
		timerLUSol1 = omp_get_wtime();
		#endif
		
		mkl_dcsrtrsv (&lower1, &lower, &lower2, &n, ilu->val, ilu->ia, ilu->ja, r, y);
		mkl_dcsrtrsv (&upper1, &upper, &upper2, &n, ilu->val, ilu->ia, ilu->ja, y, rm3);
		
		#ifdef TIMER
		timerLUSol += omp_get_wtime() - timerLUSol1;
		timerSpMV1 = omp_get_wtime();
		#endif
		
		mkl_dcsrgemv (&lower, &n, mat->val, mat->ia, mat->ja, rm3, y);
		
		#ifdef TIMER
		timerSpMV += omp_get_wtime() - timerSpMV1;
		timerVector1 = omp_get_wtime();
		#endif
		
		temp = cblas_ddot(n, y, 1, r, 1);
		temp1 = cblas_ddot(n, y, 1, y, 1);
		omega = temp / temp1;
		
		cblas_daxpy (n, omega, rm3, 1, x, 1);
		cblas_daxpy (n, -omega, y, 1, r, 1);
		nrm = cblas_dnrm2 (n, r, 1);
		#ifdef TIMER
		timerVector += omp_get_wtime() - timerVector1;
		#endif
		if ((nrm < tol) && ( nrm / nrm0 < tol )) {printf("  iteration = %3d, residual = %le \n", k+1, nrm / nrm0); break; }
		printf("  iteration = %3d, residual = %le \n", k+1, nrm / nrm0);
	}

	#ifdef TIMER
	printf("time LUSol\t%lf\ntime SpMV\t%lf\ntime v-v oper\t%lf\n",timerLUSol,timerSpMV,timerVector);
	printf("time total\t%lf\n",omp_get_wtime()-timerTotal);
	#endif
	
	free (r);
	free (rm1);
	free (rm2);
	free (rm3);
	free (zm1);
	free (zm2);
	free (p);
	free (y);
	
	return 0;
}