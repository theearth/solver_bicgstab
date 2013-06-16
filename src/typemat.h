#ifndef TYPEMAT_H
#define TYPEMAT_H
/*   Matrix Format   */

/*------- COO format */
typedef struct coo_t_ 
{
  int n;
  int nnz;
  int *ia;
  int *ja;
  double *val;
} coo_t;

/*------- CSR format */
typedef struct csr_t_
{
  int n;
  int nnz;
  int *ia;
  int *ja;
  double *val;
} csr_t;

 


#endif