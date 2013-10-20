/*
 * authors :Kranta Marieta  &&  Vouronikou Mina  &&  Konstantinidou Theodosia
 *
 ******         This file contains the mathematical functions **********
 **/



#include "types.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "functions.h"
#include <math.h>
#include <float.h>


struct node1 *dc_root;



void dc_sweep(double **pinakasG,double *x1,double *y1,double *b,int size,int flag){
  double j;
  int i,k;
  int counter=0;
  double itol = 1e-3;
  double tempb,tempc,tempa;
  double *d,*d1;
  double **S;
  struct node1 *tmp;
  FILE *fp2 =fopen("myfile.txt","w+");
  
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){  
      
    counter++;
    printf("\n\n %d) DC_sweep for %s\n\n",counter,tmp->element->type);
    

    
    if(flag == 1) { fprintf(fp2,"KANW DC_SWEEP ME CHOLESKY PARAGONTOPOIHSH\n\n"); }
    else if(flag==2){ fprintf(fp2,"KANW DC_SWEEP ME LU PARAGONTOPOIHSH\n\n"); }
    
    else{

    	d = malloc(sizeof(double)*size);
	for(i=0;i<size;i++){
	  if(pinakasG[i][i]!= 0){
		  d[i] = 1/pinakasG[i][i];
	  }else{
		  d[i] = 1.0;
	  }
	}
	
	if(flag == 4){
	    fprintf(fp2,"KANW DC_SWEEP ME Bi-CG PARAGONTOPOIHSH\n\n");
 	    d1 = malloc(sizeof(double)*size);
 	    S = calloc(sizeof(double)*size,sizeof(double));
 	    for(i=0;i<size;i++){
 	      S[i]= calloc(sizeof(double)*size,sizeof(double));  
 	    }
	    
	    for(i=0;i<size;i++){
	      for(k=0;k<size;k++){
		  S[k][i] = pinakasG[i][k];
	      }
 	    }
     
 	    for(i=0;i<size;i++){
		d1[i] = 1/S[i][i];
 	    }
 	    
	}
	else{ fprintf(fp2,"KANW DC_SWEEP ME CG PARAGONTOPOIHSH\n\n"); }
    }
    
    
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){tempb=b[tmp->element->pos];}
      if(tmp->element->pos2 != -1){tempc=b[tmp->element->pos2];}
    }
    if(tmp->element->type[0] == 'V'){tempa=b[tmp->element->pos];}
    
    
    for(j=tmp->element->start;j<tmp->element->end;j=j+tmp->element->step){
      
      
      if(tmp->element->type[0] == 'I'){
	if(tmp->element->pos != -1){
	  b[tmp->element->pos]= j;
	  if(tmp->element->pos2 != -1){
	    b[tmp->element->pos2]=-j;
	  }
	}
      }
      
      if(tmp->element->type[0] == 'V'){ b[tmp->element->pos]=j;}
      
      if(flag == 1 || flag == 2){
	y1 = forw_sub(pinakasG,b,size,y1,flag);	   
	x1 = back_sub(pinakasG,y1,size,x1,flag);
      }
      else if(flag == 3 || flag == 4){
	  for(i=0;i<size;i++){
    
	    x1[i] = 0.0;
      
	  }
	if(flag == 3){
	  CG(x1,pinakasG,b,itol,size,d);
	  printf("\n");
	}
	else{  
	  Bi_CG(x1,pinakasG,b,itol,size,d,d1,S);
	  
	}
      }
      
      //fprintf(fp2,"to neo b[%d] einai : %lf\n\n",tmp->element->pos,b[tmp->element->pos]);
      fprintf(fp2,"Gia vhma %lf exoume %s[%d] = %lf\n",j,tmp->element->tplot,tmp->element->value_plot,x1[tmp->element->value_plot-1]);
      //fprintf(fp2,"Gia vhma %lf exoume %s[%d] = %lf\n",j,tmp->element->tplot,tmp->element->value_plot,x1[9]);
     
      
    }
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){ b[tmp->element->pos]= tempb;}
      if(tmp->element->pos2 != -1){b[tmp->element->pos2]= tempc;}
    }
    if(tmp->element->type[0] == 'V'){b[tmp->element->pos]=tempa;}
    
    
    
    for(i=0;i<size;i++){
      printf("x[%d] :%lf\n",i,x1[i]);
    }
    
    printf("\n");
    for(i=0;i<size;i++){
      printf("y[%d] : %lf\n",i,y1[i]);
    }
}
fclose(fp2);
}
/************************************************************************************************/

int *LU(double **A,int size){
  int k,i,m,j;
  double x,temp1;
  double *temp2;
  int *p,*ret;
  m=0;
  
  p = malloc(sizeof(int)*size);
  ret = p;
  for(i=0;i<size;i++){
    p[i]=i;    			//pinakasG metathesis p
  }
  
  for(k=0;k<size;k++){
    x=fabs(A[k][k]);
    
    for(i=k;i<size;i++){
      
      if(fabs(A[i][k])>=x){
	
	x = A[i][k];
	m = i;
	
	temp2 = A[k];
	A[k] = A[m];
	A[m] = temp2;
	
	
	temp1 = p[k];
	p[k] = p[m];
	p[m] = temp1;
      }
      
    }
    
    
    
    
    for(i=k+1;i<size;i++){
      A[i][k] = A[i][k]/A[k][k];    //Lik
    }
    for(i=k+1;i<size;i++){
      for(j=k+1;j<size;j++){
	A[i][j] =  A[i][j] -  A[i][k] * A[k][j];
      }
      
    }
    
    
    
  }
  printf("\nO LU pinakasG einai : \n");
  printf("\n");
  
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      printf("%.3lf\t",A[i][j]);
    }
    printf("\n");
  }   
  free(p);
  return(ret);
  
}
/*******************************************************************************************************************/
//reverse the b vector
double* reverse_b(int *p,double *b,int size){
  int i;
  double *temp;
  
  temp = malloc(sizeof(double)*size);
  
  for(i=0;i<size;i++){
    temp[i]=b[p[i]];
  }
  
  return(temp); 
}

/**********************************************************************************/
//Cholesky
void Cholesky(double **A,int size){
  int k,i,j;

  double m;
  
  
  for(k=0;k<size;k++){
    
   
    
    A[k][k] =(double)sqrt(A[k][k]);
   
    m = A[k][k];

    if(m<0){
      printf("Den einai thetika orismenos pinakasG\n"); 
      exit(0);
    }
  
    for(i=k+1;i<size;i++){ 
       A[i][k] = (A[i][k])/(A[k][k]);
    }
    
      for(j=k+1;j<size;j++){ 
	    for(i=j;i<size;i++){ 
	       A[i][j] =  A[i][j] -  (A[i][k] * A[j][k]);   
	    }
      }
  
  
 
  }
  
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      printf(" %lf\t",A[i][j]); 
    }
    printf("\n");
  }
}

/**********************************************************************************/
// Forward solve Ly = b


double* forw_sub(double **L, double *b, int size,double *y,int flag){
  
  int k,j;
  for(k=0;k<size;k++){
    y[k] = b[k];
    for(j=0;j<k;j++){
      
      y[k] = y[k] - (L[k][j]*y[j]);
    }
    if(flag == 1){
      y[k] =y[k]/L[k][k];
    }
  }
  return(y);
}


/**********************************************************************************/
// Backward solve Ux = y

double* back_sub(double **U, double *y, int size,double *x,int flag){  
  int k,j;
  
  for(k=size-1;k>=0;k--){
    x[k] = y[k];
    for(j=k+1;j<size;j++){
      
     
      if(flag == 1){
	x[k] = x[k] - U[j][k] * x[j];
	
      }
      else{
	x[k] = x[k] - U[k][j] * x[j]; 
      }
      
    }
    x[k] = x[k]/U[k][k];  
    
  }
  
  return(x);
}



/*****************************************************************************/

double *CG(double *x,double **A,double *b,double itol,int size,double *d){
	
	double alpha,beta,omega,lala,alpha2;
	int i,k;
	double *z,*p,*Ap,*r,*Ax,*pA;
	
	
	r = malloc(sizeof(double)*size);
	p = malloc(sizeof(double)*size);
	
	Ax = product(A,x,size);
	printf("\n");
	
	for(i=0;i<size;i++){
		r[i] = b[i] - Ax[i];
	}
		
	z = precond_sol(d,r,size);
	
	for(i=0;i<size;i++) p[i] = z[i]; //copy z to p
	
	k=0;
	
	while (1){
		
		k=k+1;
		pA = product2(p,A,size);
		alpha = 0;
		for(i=0;i<size;i++){
			alpha = alpha + r[i]*z[i];
		}
		beta = 0;
		for(i=0;i<size;i++){
			beta = beta + pA[i]*p[i];
		}
		
		omega = alpha/beta;
		
		Ap = product(A,p,size);
		for(i = 0;i<size;i++){
			x[i] = x[i] + omega*p[i];
			r[i] = r[i] - omega*Ap[i];
		}
		
		if (norm(r,size)/norm(b,size)<itol || k == size) break;
		
		z = precond_sol(d,r,size);
		
		alpha2 = 0;
		for(i=0;i<size;i++){
			alpha2 = alpha2 + r[i]*z[i];
		}
		
		lala = alpha2/alpha;
		
		for(i=0;i<size;i++){
			p[i] = z[i] + lala*p[i];
		}
		
	}
		
		for(i = 0;i<size;i++){
		printf("x[%d] = %lf\n",i,x[i]); 
		
	}	
	printf("\n\n");	
	return(x);
}

/*****************************************************************************/
void Bi_CG(double *x,double **A,double *b,double itol,int size,double *d,double *d1,double **S){
  
  
  double iter,nb,rho,rho1,alpha,beta,omega;
  int i;
  double *z,*z1,*p,*p1,*q,*q1,*r,*r1,*Ax;
 

  r = malloc(sizeof(double)*size);
  r1 = malloc(sizeof(double)*size);
  p = malloc(sizeof(double)*size);
  p1 = malloc(sizeof(double)*size);


     
  Ax = product(A,x,size);
  for(i=0;i<size;i++){
    r[i] = b[i] - Ax[i];
    r1[i] = r[i];
  }
 
  iter = 0.0;
  rho = 0.0;
  nb = norm(b,size);
  
  
   while (1){ 
    iter = iter + 1;
    z = precond_sol(d,r,size);
    z1 = precond_sol(d1,r1,size);
    
    rho = 0;
    for(i = 0;i<size;i++){
      rho = rho + z[i]*r1[i];
    }
    
    if(fabs(rho) < DBL_EPSILON){ exit(2); }
    if(iter == 1){
      for(i=0;i<size;i++){
	p[i] = z[i];
	p1[i] = z1[i];
      }
    }
    else{
      beta = rho/rho1;
      for(i=0;i<size;i++){
	p[i] = z[i] + beta*p[i];
	p1[i] = z1[i] + beta*p1[i];
      }
    }
    
    rho1 = rho;
    q = product(A,p,size);
    q1 = product(S,p1,size);
    
    omega = 0;
    for(i = 0;i<size;i++){
      omega = omega + q[i]*p1[i];
    }
    
    if(fabs(omega) < DBL_EPSILON){ exit(2); }
    alpha = rho/omega;
    
    for(i = 0;i<size;i++){
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*q[i];
      r1[i] = r1[i] - alpha*q1[i];
    }
    
      if (norm(r,size)/norm(b,size)<itol ||iter == size) break;
  }
  printf("Mesa sti BCG kai kanw %lf epanalipseis\n",iter);
  for(i = 0;i<size;i++){
     printf("To x[%d] = %lf\n",i,x[i]);
    }
    
  free(r);
  free(r1);
  free(p);
  free(p1);

}
/**********************************************************************************/

double norm(double *vec,int size){
  double sum;
  int i;
  
  sum = 0.0;
  
  for(i=0;i<size;i++){
    sum = sum + pow(vec[i],2);
  }
  sum = sqrt(sum);
  return(sum);
  
}

/**********************************************************************************/

double *precond_sol(double *d,double *r,int size){
  int i;
  double *z;
  z = malloc(sizeof(double)*size); 
  
  for(i=0;i<size;i++){
    z[i] = d[i]*r[i];
  }
  return(z);
}


/**********************************************************************************/
double* product(double **A,double *p,int size){
  double sum_a = 0.0;
  int i,j;
  double *q;
  
   q = malloc(sizeof(double)*size);
  
  for(i=0;i<size;i++){
    sum_a=0.0;
    for(j=0;j<size;j++){
      sum_a = sum_a + A[i][j]*p[j];
    }
    q[i] = sum_a;
  }
  return(q);
}

/*************************************************************************************************************************/

double* product2(double *p,double **A,int size){
	double sum_a = 0.0;
	int i,j;
	double *q;
	
	q = malloc(sizeof(double)*size);
	
	for(i=0;i<size;i++){
		sum_a=0.0;
		for(j=0;j<size;j++){
			sum_a = sum_a + p[j]*A[j][i];
		}
		q[i] = sum_a;
	}
	return(q);
}