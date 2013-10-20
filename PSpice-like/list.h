/*
 * authors :Kranta Marieta  &&  Vouronikou Mina  &&  Konstantinidou Theodosia
 *
 ******         This file contains headers of functions that handle the list  **********
 **/

 #include "./CSparse/Include/cs.h"
 
int change_nodes(int m2);
void change_name(int m2);
void add_list(struct elem *element);
void print_list();
void clear_list();
void check_list();
void MNAG(double **pinakasG,int size,double *b,char **vector,double max,int m2,struct dc_nodes *elem1);
int mapa();
int sort(double *map,int k);
void dc_add_list(struct dc_nodes *element);
void dc_print_list();
void clear_list2();
int check_list2();
void s_p_a_r_s_e(int max,int m2,int flag,int flag1);

//assign 5
double *sparse_CG(double *x,cs *A,double *b,double itol,int size,double *d);
double *sparse_Bi_CG(double *x,cs *A,double *b,double itol,int size,double *d);
double *lu_sparse(cs *G,int plithos_nonz,int n,double *b);
double *chol_sparse(cs *G,int plithos_nonz,int n,double *b);
double* mul(cs *A,double *l,int size);
double* mul2(double *l,cs *A,int size);
double* diag(cs *A,int size);
void dc_sweep_sparse(cs *pinakasG,int non_z,int size,double *b,int flag);
//assign 6
double value_for_time(struct dc_nodes *element , double t, double finalTime);
double linear1(double t1, double t2, double val1, double val2,double t);
double linear2(double t1, double t2, double val1, double val2,double t);
void MNAC(double **pinakasC,int size,double *b,char **vector,double max,int m2,struct dc_nodes* elem1);
double *BE(cs *G,cs *C,double *x,double *b,double itol,int size,int plithos_nonz,double h);

double *TR(cs *G,cs *C,double *x,double *b,double *oldb,double itol,int size,int plithos_nonz,double h);
double * back_euler(double **G,double **C,double *x1,double *b,int size,double h,double **A,double *tmp,double *da,double *y);
double * trapezoidal(double **G,double **C,double *x1,double *b,double *b1,int size,double h,double **A,double **B,double *tmp,double *da,double *y);
void transient_analysis(double **pinakasG,double **pinakasC,int flag,double *x1,double *b,int size,double itol);
 void transient_sparse(cs *pinakasG,cs *pinakasC,int non_z,double *k,int size,double *b,int flag);