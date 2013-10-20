/*
 * authors :Kranta Marieta  &&  Vouronikou Mina  &&  Konstantinidou Theodosia
 *
 ******         This file contains the headers of mathematical functions   **********
 **/



//assign3

void dc_sweep(double **pinakasG,double *x1,double *y1,double *b,int size,int flag);
int *LU(double **A,int size);
double* reverse_b(int *p,double *b,int size);
void Cholesky(double **A,int size);
double* forw_sub(double **L, double *b, int size,double* y,int flag);
double* back_sub(double **U, double *y, int size,double* x,int flag);

//assign4
double *CG(double *x,double **A,double *b,double itol,int size,double *d);
void Bi_CG(double *x,double **A,double *b,double itol,int size,double *d,double *d1,double **S);
double norm(double *vec,int size);
double * precond_sol(double *d,double *r,int size);
double* product(double **A,double *p,int size);
double* product2(double *p,double **A,int size);