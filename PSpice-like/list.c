/*
 * authors :Kranta Marieta  &&  Vouronikou Mina  &&  Konstantinidou Theodosia
 */
#include "types.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "list.h"
#include <math.h>
#include <float.h>
#include "functions.h"
#include "./CSparse/Include/cs.h"
#include <ctype.h>
#include "cs.h"
#include <unistd.h>

#define PAIRS 5


struct node *root;
struct node1 *dc_root;
char str[10];
int nonz = 0;

//each node is a circuit element 
void add_list(struct elem *element){
  struct node *tmp;
  tmp = (struct node *)malloc(sizeof(struct node));
  tmp->element = element;
  tmp->nxt = root;
  root = tmp;
}


int change_nodes(int m2){
  
  struct node *tmp;
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    
    if (tmp->element->type[0] == 'L') {
      m2++;
    }
  }
  return m2;  
}

void change_name(int m2){
  
  struct node *tmp;
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    
    if (tmp->element->type[0] == 'L') {
      sprintf(tmp->element->type,"V%d",m2);
      m2--;
    }
  }
}


void clear_list(){
  struct node *tmp;
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    free(tmp->element->term);
    free(tmp->element->value);
  }
  free(tmp);
}


void print_list(){
  struct node *tmp;
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    printf("%s\n",tmp->element->type);
    printf("%lf\n",tmp->element->term[0]);
    printf("%lf\n",tmp->element->term[1]);
  }
}

void check_list(){
  struct node *tmp;
  struct node *tmp1;
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    for(tmp1 = root->nxt; tmp1!= (struct node *) 0; tmp1=tmp1->nxt){
      if ((strcmp(tmp->element->type,tmp1->element->type) == 0) && tmp!=tmp1){
	printf("%s invalid element has been redeclared\n",tmp1->element->type);
      }
    }
  }
}


int check_list2(){
  struct node *tmp;
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    if((tmp->element->type[0] == 'V') || (tmp->element->type[0] == 'R')){
      if(tmp->element->term[0]!=0 && tmp->element->term[1]!=0){
	nonz = nonz +4;
	
      }
      else if((tmp->element->term[0]!=0 && tmp->element->term[1]==0) || (tmp->element->term[0]==0 && tmp->element->term[1]!=0)){
	nonz = nonz+1;
	
      }
    }
  }
  
  return (nonz);
}
// 
/**********************************************************************************************************/
int check_list3(){
  struct node *tmp;
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    if((tmp->element->type[0] == 'L') || (tmp->element->type[0] == 'C')){
      if(tmp->element->term[0]!=0 && tmp->element->term[1]!=0){
	nonz = nonz +4;
	
      }
      else if((tmp->element->term[0]!=0 && tmp->element->term[1]==0) || (tmp->element->term[0]==0 && tmp->element->term[1]!=0)){
	nonz = nonz+1;
	
      }
    }
  }
  
  return (nonz);
}

//clears the dc_list
void clear_list2(){
  struct node1 *tmp;
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){
    free(tmp->element->type);
    free(tmp->element->tplot);
  }
  free(tmp);
}

//print list dc
void dc_print_list(){
  struct node1 *tmp;
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){
    printf("%s\n",tmp->element->type);
    printf("%lf\n",tmp->element->start);
    printf("%lf\n",tmp->element->end);
    printf("%lf\n",tmp->element->step);
    printf("%d\n",tmp->element->pos);
  }
}

// adds elements in the dc_list
void dc_add_list(struct dc_nodes *element1){
  
  struct node1 *tmp;
 
  tmp = (struct node1 *)malloc(sizeof(struct node1));
 
  tmp->element = element1;
  tmp->nxt = dc_root;
  dc_root = tmp;
  
}


/******************************************************************************************************************/

void MNAG(double **pinakasG,int size,double *b,char **vector,double max,int m2,struct dc_nodes* elem1){
  struct node *tmp;
  struct node1 *tmp2;
  int k=0;
  int i,j,temp;
  
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    if(tmp->element->type[0] == 'R'){
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	pinakasG[(int)tmp->element->term[0]-1][(int)tmp->element->term[0]-1] +=  1/tmp->element->value[0];	  
	pinakasG[(int)tmp->element->term[0]-1][(int)tmp->element->term[1]-1] += -1/tmp->element->value[0];
	pinakasG[(int)tmp->element->term[1]-1][(int)tmp->element->term[0]-1] += -1/tmp->element->value[0];
	pinakasG[(int)tmp->element->term[1]-1][(int)tmp->element->term[1]-1] +=  1/tmp->element->value[0];
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	pinakasG[(int)tmp->element->term[1]-1][(int)tmp->element->term[1]-1] +=  1/tmp->element->value[0];
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	pinakasG[(int)tmp->element->term[0]-1][(int)tmp->element->term[0]-1] += 1/tmp->element->value[0];
      }
    }
    if(tmp->element->type[0] == 'I'){
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	b[(int)tmp->element->term[0]-1] += -tmp->element->value[0];
	b[(int)tmp->element->term[1]-1] += tmp->element->value[0];
	
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = tmp->element->term[1]-1;
	    tmp2->element->pos2 = tmp->element->term[0]-1;
	  }
	}
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	b[(int)tmp->element->term[1]-1] += tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = tmp->element->term[1]-1;
	  }
	}
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	b[(int)tmp->element->term[0]-1] += -tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos2 = tmp->element->term[0]-1;
	  }
	}
      }
    }
    
    if(tmp->element->type[0] == 'V'){
      k++;
      
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	pinakasG[(int)tmp->element->term[0]-1][(int)max+k-1] +=  1;
	pinakasG[(int)tmp->element->term[1]-1][(int)max+k-1] += -1;
	pinakasG[(int)max+k-1][(int)tmp->element->term[0]-1] +=  1;
	pinakasG[(int)max+k-1][(int)tmp->element->term[1]-1] += -1;
	b[(int)max+k-1] = tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = max+k-1;
	  }
	}
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	pinakasG[(int)tmp->element->term[1]-1][(int)max+k] += -1;
	pinakasG[(int)max+k-1][(int)tmp->element->term[1]-1] += -1;
	b[(int)max+k-2] = tmp->element->value[0];
	
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = max+k-2;
	  }
	}
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	pinakasG[(int)tmp->element->term[0]-1][(int)max+k-1] += 1;
	pinakasG[(int)max+k-1][(int)tmp->element->term[0]-1] += 1;
	b[(int)max+k-1] = tmp->element->value[0];
	
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = max+k-1;
	  }
	}  
      }
    }
    
    
  }
  printf("O pinakasG A einai :\n\n");
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      printf("%.3lf\t",pinakasG[i][j]);    
    }
    printf("\n");
  }
  
  printf("To dianusma  b einai :\n\n");
  
  for(i=0;i<size;i++){
    printf("%lf\n",b[i]);
  }
  for(i=0;i<max;i++){
    
    sprintf(vector[i],"V%d",i);
  }
  for(i=max;i<max+m2;i++){
    temp = (int)i-max+1;
    sprintf(vector[i],"I%d",temp);
    
  }
  printf("To dianusma  vector einai :\n\n");
  
  for(i=0;i<max+m2;i++){
    printf("%s\n",vector[i]);
  }
}
/****************************************************************************************************************/
void MNAC(double **pinakasC,int size,double *b,char **vector,double max,int m2,struct dc_nodes* elem1){
  struct node *tmp;
  int k=0;
  int i,j;

  
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    if(tmp->element->type[0] == 'C'){
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	pinakasC[(int)tmp->element->term[0]-1][(int)tmp->element->term[0]-1] +=  tmp->element->value[0];	  
	pinakasC[(int)tmp->element->term[0]-1][(int)tmp->element->term[1]-1] += -tmp->element->value[0];
	pinakasC[(int)tmp->element->term[1]-1][(int)tmp->element->term[0]-1] += -tmp->element->value[0];
	pinakasC[(int)tmp->element->term[1]-1][(int)tmp->element->term[1]-1] +=  tmp->element->value[0];
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	pinakasC[(int)tmp->element->term[1]-1][(int)tmp->element->term[1]-1] +=  tmp->element->value[0];
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	pinakasC[(int)tmp->element->term[0]-1][(int)tmp->element->term[0]-1] += tmp->element->value[0];
      }
    }
    
    
    if(tmp->element->type[0] == 'L'){
      k++;
      pinakasC[(int)max+k-1][(int)max+k-1] = -(tmp->element->value[0]);
      
    }
    
  }
  printf("O pinakasC einai :\n\n");
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      printf("%.4lf\t",pinakasC[i][j]);    
    }
    printf("\n");
  }
  
  
}

/**********************************************************************************/
/*mapping*/

int mapa(){
  
  struct node *tmp;
  double *map;
  int k=0;
  int i=0;
  int temp2;
  
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    
    if(tmp->element->type[0] == 'M'){k=k+2;}
    if(tmp->element->type[0] == 'Q'){k=k+1;}
    k=k+2;
    
  }
  map =malloc(sizeof(double)*k);
 
  for(tmp = root,i=0; (tmp!= (struct node *) 0); tmp=tmp->nxt,i=i+2){
    
    map[i]=tmp->element->term[0];
    map[i+1]=tmp->element->term[1];
  }
  
  printf("\n\n\n");
  
  temp2=sort(map,k);
  
  free(map);
  
  return(temp2);                         //returns max 
  
}

/**********************************************************************************/
//sort the table
int sort(double *map,int k){
  int i,j,counter,temp2;
  double temp,*map2;
  struct node *tmp;
  temp2=0;
  
  for(i=0;i<k-1;i++){
    for(j=0;j<k-i-1;j++){
      if (map[j] > map[j+1]){
	temp = map[j+1];
	map[j+1] = map[j];
	map[j] = temp;	  
      }
    }
  }
  
  for(i=0;i<k;i=i+counter+1){
    temp2++;
    counter = 0;
    for(j=i+1;j<k;j++){
      if (map[i] == map[j]){ 
	counter++; 
      }
    }
  }
  
  map2 =malloc(sizeof(double)*(temp2+1));
  map2[0] =map[0];
  i = 0;
  for(j=0;j<k;j++){
    if (map[j] != map2[i]){ 
      i++;
      map2[i] = map[j];
    }
  }
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    for(i=0;i<temp2;i++){
      if(tmp->element->term[0] == map2[i]){
	tmp->element->term[0] = (double)i;
      }
      if(tmp->element->term[1] == map2[i]){
	tmp->element->term[1] = (double)i;
      }
    }
  }
  //print_list();
  free(map2);
  
  return(temp2);
}
/******************** sparse MNAG ***************************/
void s_p_a_r_s_e(int max,int m2,int flag,int flag1){
  
  
  struct node *tmp;
  struct node1 *tmp2;
  cs *A,*B,*G,*C;
  int plithos_nonz;
  int plithos_nonz2;
  int n,k,i;
  double *b,*x;
  double *d,*temp;
  
  
  plithos_nonz = check_list2();		// G
  plithos_nonz2 = check_list3();	// C
  n = (max)+m2;
  
  temp = cs_malloc(n,sizeof(double));
  x = cs_malloc(n,sizeof(double));
  d = cs_malloc(n,sizeof(double));
  b = cs_malloc(n,sizeof(double));
  A = cs_spalloc(n,n,plithos_nonz,1,1);
  B = cs_spalloc(n,n,plithos_nonz2,1,1);
  
  
  k = 0;
  
  A->nz = plithos_nonz;
  B->nz = plithos_nonz2;
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    if(tmp->element->type[0] == 'C'){
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	if( (cs_entry(B,(int)tmp->element->term[0]-1,(int)tmp->element->term[0]-1,tmp->element->value[0]) == 1)&&         
	  (cs_entry(B,(int)tmp->element->term[0]-1,(int)tmp->element->term[1]-1,-tmp->element->value[0]) == 1)&&
	  (cs_entry(B,(int)tmp->element->term[1]-1,(int)tmp->element->term[0]-1,-tmp->element->value[0]) == 1)&&
	  (cs_entry(B,(int)tmp->element->term[1]-1,(int)tmp->element->term[1]-1,tmp->element->value[0]) == 1) )
	{continue;}
	else{printf("wrong-error");break;}
      }
    }
    else if(tmp->element->type[0] == 'L'){
      k++;
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	if(cs_entry(B,(int)max+k-1,(int)max+k-1,tmp->element->value[0])){continue;}
      } 
    }
  }
  
  k = 0;
  change_name(m2);
  
  for(tmp = root; tmp!= (struct node *) 0; tmp=tmp->nxt){
    
    if(tmp->element->type[0] == 'R'){
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	if( (cs_entry(A,(int)tmp->element->term[0]-1,(int)tmp->element->term[0]-1,1/tmp->element->value[0]) == 1)&&         
	  (cs_entry(A,(int)tmp->element->term[0]-1,(int)tmp->element->term[1]-1,-1/tmp->element->value[0]) == 1)&&
	  (cs_entry(A,(int)tmp->element->term[1]-1,(int)tmp->element->term[0]-1,-1/tmp->element->value[0]) == 1)&&
	  (cs_entry(A,(int)tmp->element->term[1]-1,(int)tmp->element->term[1]-1,1/tmp->element->value[0]) == 1) )
	{continue;}
	else{printf("wrong-error");break;}
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	if(cs_entry(A,(int)tmp->element->term[1]-1,(int)tmp->element->term[1]-1,1/tmp->element->value[0]) == 1){continue;}
	else{printf("error"); break;}
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	if(cs_entry(A,(int)tmp->element->term[0]-1,(int)tmp->element->term[0]-1,1/tmp->element->value[0]) == 1){continue;}
	else{printf("error"); break;}
      }
    }
    
    if(tmp->element->type[0] == 'I'){
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	b[(int)tmp->element->term[0]-1] += -tmp->element->value[0];
	b[(int)tmp->element->term[1]-1] += tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = tmp->element->term[1]-1;
	    tmp2->element->pos2 = tmp->element->term[0]-1;
	  }
	}
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	b[(int)tmp->element->term[1]-1] += tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = tmp->element->term[1]-1;
	  }
	}
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	b[(int)tmp->element->term[0]-1] += -tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos2 = tmp->element->term[0]-1;
	  }
	}
      }
    }
    
    if(tmp->element->type[0] == 'V'){
      k++;
      if((tmp->element->term[0] != 0) && (tmp->element->term[1]!=0)){
	if( (cs_entry(A,(int)tmp->element->term[0]-1,(int)max+k-1,1) == 1)&&
	  (cs_entry(A,(int)tmp->element->term[1]-1,(int)max+k-1,-1) == 1)&&
	  (cs_entry(A,(int)max+k-1,(int)tmp->element->term[0]-1,1) == 1)&&
	  (cs_entry(A,(int)max+k-1,(int)tmp->element->term[1]-1,-1) == 1) )
	{continue;}
	else{printf("error"); break;}
	b[(int)max+k-1] = tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = max+k-1;
	  }
	}
      }
      else if((tmp->element->term[0] == 0) && (tmp->element->term[1]!=0)){
	if( (cs_entry(A,(int)tmp->element->term[1]-1,(int)max+k,-1) == 1)&&
	  (cs_entry(A,(int)max+k-1,(int)tmp->element->term[1]-1,-1) == 1) )
	{continue;}
	else{printf("error"); break;}
	b[(int)max+k-2] = tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = max+k-2;
	  }
	}
      }
      else if((tmp->element->term[0] != 0) && (tmp->element->term[1]==0)){
	if( (cs_entry(A,(int)tmp->element->term[0]-1,(int)max+k-1,1) == 1)&&
	  (cs_entry(A,(int)max+k-1,(int)tmp->element->term[0]-1,1) == 1) )
	{continue;}
	else{printf("error"); break;}
	
	b[(int)max+k-1] = tmp->element->value[0];
	for(tmp2 = dc_root; tmp2!= (struct node1 *) 0; tmp2=tmp2->nxt){
	  if(!strcmp(tmp2->element->type,tmp->element->type)){
	    tmp2->element->pos = max+k-1;
	  }
	}	 
      }
    }
    
    
  }
  
  
  G = cs_compress((cs *)A); 
  C = cs_compress((cs *)B);
  
  cs_spfree(A);
  cs_spfree(B);
  
  cs_dupl(G);
  cs_dupl(C);
  
  cs_print(G,1);
  cs_print(C,1);
  
  for(i=0;i<n;i++){
    
    temp[i] = b[i];
  }
  
  
  
  if (flag == 2){
    printf("***************************    Sparse LU decomposition    ******************************************\n");
    lu_sparse(G,plithos_nonz,n,temp);
    
  
  
  }
  else if(flag == 1){
    printf("***************************  Sparse CHOLESKY decomposition ******************************************\n");
    
    
    chol_sparse(G,plithos_nonz,n,temp);
  }
  else if(flag == 4){
    d = diag(G,n);
  printf("***************************  Sparse  CG  decomposition     ******************************************\n");

   for(i=0;i<n;i++) {x[i] = 0.0;}  
   sparse_CG(x,G,b,1e-3,n,d);

  }
  
  else if(flag == 3){
    d = diag(G,n);
    printf("***************************   Sparse Bi-CG   decomposition ******************************************\n");
    for(i=0;i<n;i++) {x[i] = 0.0;}  
    sparse_Bi_CG(x,G,b,1e-3,n,d);  
  }
  
  else if(flag == 8 || flag == 7){
     d = diag(G,n);
     for(i=0;i<n;i++) {x[i] = 0.0;}  
     sparse_CG(x,G,b,1e-3,n,d);

    if(flag == 8){
    
      printf("***************************  Sparse  BE    ******************************************\n");
      x = BE(G,C,x,b,1e-3,n,plithos_nonz,0.1);
   
    }
   else{
    
   printf("***************************    Sparse  TR    ******************************************\n");
   
   for(i=0;i<n;i++){
      if(i==n-3){temp[i]=0.1;}
      else{temp[i]=b[i];}
   }
    x = TR(G,C,x,b,temp,1e-3,n,plithos_nonz,0.1); 
  }
  
}

  if(flag1==1){
    printf("***************************  Sparse dc sweep analysis      ******************************************\n");
    dc_sweep_sparse(G,plithos_nonz,n,b,flag);
  }
  else {
    	transient_sparse(G,C,plithos_nonz,x,n,b,flag);

  }
  
}

/**********************************************************************************************************/

double *lu_sparse(cs *G,int plithos_nonz,int n,double *b){
  
  css *S;
  csn *N;
  int j;
  double *x;
  
  
  
  
  S = cs_malloc(n,sizeof(css *));
  N = cs_malloc(n,sizeof(csn *));
  x = cs_malloc(n,sizeof(double));
  
  
  S = cs_sqr(2,G,0);
  N = cs_lu(G,S,1);
  
  
  cs_ipvec(N->pinv,b,x,n);
  cs_lsolve(N->L,x);
  cs_usolve(N->U,x);
  cs_ipvec(S->q,x,b,n);
  
  for(j=0;j<n;j++){
    printf("x[%d] = %lf\n",j,b[j]);
  }
  printf("\n");
  
  
  return(b);
}
/***************************************************************************************************************/
double *chol_sparse(cs *G,int plithos_nonz,int n,double *b){
  
  css *S;
  csn *N;
  int j;
  double *x;
  
  
  S = cs_malloc(n,sizeof(css *));
  N = cs_malloc(n,sizeof(csn *));
  x = cs_malloc(n,sizeof(double));
  S = cs_schol(1,G);
  N = cs_chol(G,S);
  
  cs_ipvec(S->pinv,b,x,n);
  cs_lsolve(N->L,x);
  cs_ltsolve(N->L,x);
  cs_pvec(S->pinv,x,b,n);
  
  for(j=0;j<n;j++){
    printf("x[%d] = %lf\n",j,b[j]);
  }
  printf("\n");
  
  
  return(b);
}

/******************************************************************************************************/
double *sparse_CG(double *x,cs *A,double *b,double itol,int size,double *d){
  
  double alpha,beta,omega,lala,alpha2;
  int i,k;
  double *z,*p,*Ap,*r,*Ax,*pA;
  
  
  r = malloc(sizeof(double)*size);
  p = malloc(sizeof(double)*size);
  
  
  Ax = mul(A,x,size);
  printf("\n");
  
  for(i=0;i<size;i++){
    r[i] = b[i] - Ax[i];
  }
  
  z = precond_sol(d,r,size);
  
  for(i=0;i<size;i++) p[i] = z[i]; //copy z to p
	
	k=0;
  
  while (1){
    
    k=k+1;
    pA = mul(A,p,size);
    alpha = 0;
    for(i=0;i<size;i++){
      alpha = alpha + r[i]*z[i];
    }
    beta = 0;
    for(i=0;i<size;i++){
      beta = beta + pA[i]*p[i];
    }
    
    omega = alpha/beta;
    
    Ap = mul(A,p,size);
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
  
  free(r);
  free(p);
  
  return(x);
}

/**********************************************************************************/
double *sparse_Bi_CG(double *x,cs *A,double *b,double itol,int size,double *d){
  
  
  double iter,nb,rho,rho1,alpha,beta,omega;
  int i;
  double *z,*z1,*p,*p1,*q,*q1,*r,*r1,*Ax;
  
  
  r = malloc(sizeof(double)*size);
  r1 = malloc(sizeof(double)*size);
  p = malloc(sizeof(double)*size);
  p1 = malloc(sizeof(double)*size);
  
  
  Ax = mul(A,x,size);
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
    z1 = precond_sol(d,r1,size);
    
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
    q = mul(A,p,size);
    q1 = mul2(p1,A,size);
    
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
    
    if (norm(r,size)/nb <itol ||iter == size) break;
  }
  
  printf("Mesa sti BCG kai kanw %lf epanalipseis\n",iter);
  for(i = 0;i<size;i++){
    printf("To x[%d] = %lf\n",i,x[i]);
  }
  
  free(r);
  free(r1);
  free(p);
  free(p1);
  
  return(x);
  
}
/**********************************************************************************/
double* mul(cs *A,double *l,int size){
  int k,j;
  double *y;
  y = calloc(sizeof(double)*size,sizeof(double));
  
  for(j=0;j<size;j++){
    for(k=A->p[j];k<A->p[j+1];k++){
      y[A->i[k]] = y[A->i[k]] + A->x[k]*l[j];
    }
  }
  
  return(y);
}

/**********************************************************************************/
double* mul2(double *l,cs *A,int size){
  int k,j;
  double *y;
  y = calloc(sizeof(double)*size,sizeof(double));
  
  for(j=0;j<size;j++){
    for(k=A->p[j];k<A->p[j+1];k++){
      y[j] = y[j] + A->x[k]*l[A->i[k]];
    }
  }
  return(y);
}
/*******************************************************************************/  
double* diag(cs *A,int size){
  int k,j;
  double *d;
  int cnt;
  
  cnt=0;
  d = malloc(size*sizeof(double));
  
  for(j=0;j<size;j++){
		d[j] = 1.0;
  }
  
  for(j=0;j<size;j++){
    for(k=A->p[j];k<A->p[j+1];k++){
      if(A->i[k]==cnt){
	d[j]=A->x[k];
      }
    }
    cnt++;
  }
  
  return(d);
}  

/*************************************** dc sweep analysis **************************************************/

void dc_sweep_sparse(cs *pinakasG,int non_z,int size,double *b,int flag){
  
  double j;
  int i;//k;
  int counter=0;
  double itol = 1e-3;
  double tempb,tempc,tempa;
  double *res;
  double *d;
  double *x1;
  double *temp;
  
  struct node1 *tmp;
  
  d = malloc(size*sizeof(double));
  x1 = malloc(non_z*sizeof(double));
  res = cs_malloc(size,sizeof(double));
  
  FILE *fp2 =fopen("myfile.txt","w+");
  
  if(fp2 == NULL) {printf("error opening file\n");}
  
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){  
    counter++;
    printf("\n\n %d) DC_sweep for %s\n\n",counter,tmp->element->type);
    
    
    if(flag == 1){ 
      fprintf(fp2,"KANW DC_SWEEP SPARSE ME CHOLESKY PARAGONTOPOIHSH\n\n"); 
    }
    else if(flag==2){
      fprintf(fp2,"KANW DC_SWEEP SPARSE ME LU PARAGONTOPOIHSH\n\n"); 
    }
    else if(flag==4){
      fprintf(fp2,"KANW DC_SWEEP ME CG PARAGONTOPOIHSH\n\n"); 
      d = diag(pinakasG,size);
      for(i=0;i<size;i++){ x1[i] = 0.0; }
    }
    if(flag == 3){
      fprintf(fp2,"KANW DC_SWEEP ME Bi-CG PARAGONTOPOIHSH\n\n");
      d = diag(pinakasG,size);
      for(i=0;i<size;i++){ x1[i] = 0.0; }
    }
    
    
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){
	tempb=b[tmp->element->pos];
      }
      if(tmp->element->pos2 != -1){
	tempc=b[tmp->element->pos2];
      }
    }
    if(tmp->element->type[0] == 'V'){
      if(tmp->element->pos != -1){
	tempa=b[tmp->element->pos];
      } 
    }
    
    
    
    
    temp = b;
    
    for(j=tmp->element->start;j<=tmp->element->end;j=j+tmp->element->step){
      
     
      res = temp; 
      
      b = res;
      
      if(tmp->element->type[0] == 'I'){
	if(tmp->element->pos != -1){
	  b[tmp->element->pos]= j;
	  if(tmp->element->pos2 != -1){
	    b[tmp->element->pos2]=-j;
	  }
	}
      }
      
      if(tmp->element->type[0] == 'V'){
	if(tmp->element->pos != -1){
	  b[tmp->element->pos]=j;
	}
      }
      
      
      if(flag==1){ // * CHOLESKY *
 	fprintf(fp2,"*******************SPRARSE CHOLESKY SWEEP*******************************\n");
	
	res = chol_sparse(pinakasG,non_z,size,b);
	
      }
      if(flag==2){ // * LU *
 	fprintf(fp2,"************************ SPARSE LU SWEEP ********************************* \n");
	res = lu_sparse(pinakasG,non_z,size,b);
	
      }
      else if(flag == 3 || flag == 4){
	
	if(flag == 4){
	  fprintf(fp2,"*************************** SPARSE CG SWEEP ******************************\n");
	  x1 = sparse_CG(x1,pinakasG,b,itol,size,d);
	  
	}
	else{  
	  fprintf(fp2,"************************* SPARSE Bi-CG SWEEP *****************************\n ");
	  x1 = sparse_Bi_CG(x1,pinakasG,b,itol,size,d);
	}
	res = x1;
      }
      
      for(i=0;i<size;i++){
	fprintf(fp2,"x[%d] = %lf \n",i,res[i]);
      }
      fprintf(fp2,"\n");
        
      
    }
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){ b[tmp->element->pos]= tempb;}
      if(tmp->element->pos2 != -1){b[tmp->element->pos2]= tempc;}
  }
  if(tmp->element->type[0] == 'V'){
    if(tmp->element->pos != -1){
      b[tmp->element->pos]=tempa;
    }
}

}
fclose(fp2);
}

/**************************************************************************************************************************************************************************************************************************/

double value_for_time(struct dc_nodes *element , double t, double finalTime){
  double res;
  double pi = 3.14;
  int i;
  //An to trans_spec einai EXP//
  if(element->flag == 1){
    if( (t >= 0)&&(t < element->trans_spec[2]) ) { return(element->trans_spec[0]); }
    else if( (t >= element->trans_spec[2]) && (t < element->trans_spec[4]) ){
      res = element->trans_spec[0] + (element->trans_spec[1]-element->trans_spec[0])*( (1 - exp(-(t-element->trans_spec[2])/element->trans_spec[3])));
      return(res);
    }
    else if( (t >= element->trans_spec[4]) && (t < finalTime) ){
      res = element->trans_spec[0] + (element->trans_spec[1]-element->trans_spec[0])*( (1 - exp(-( (t-element->trans_spec[2])/element->trans_spec[3]))) -(1 - exp(-( (t-element->trans_spec[4])/element->trans_spec[5]))) ) ;
      return(res);
    }
    else{ return(-1); }
  }
  
  
  else if(element->flag == 2){
    if( (t >= 0)&&(t < element->trans_spec[2]) ) { 
      return(element->trans_spec[0]); 
    }
    else if( (t >= element->trans_spec[2]) && (t < (element->trans_spec[2] + element->trans_spec[3]) ) ){ //prepei na mpei kai to k*per kapou
     
     res = linear1(element->trans_spec[2],(element->trans_spec[2] + element->trans_spec[3]),element->trans_spec[0],element->trans_spec[1],t);
     return(res);
    }
    else if( (t >= (element->trans_spec[2] + element->trans_spec[3]) ) && ( t < (element->trans_spec[2] + element->trans_spec[3] + element->trans_spec[5] ) ) ){
      return(element->trans_spec[1]);
    }
    else if( ( t >= (element->trans_spec[2] + element->trans_spec[3] + element->trans_spec[5] ) ) && ( t < (element->trans_spec[2] + element->trans_spec[3] + element->trans_spec[5] + element->trans_spec[4]) ) ){
      
      res = linear2((element->trans_spec[2] + element->trans_spec[3] + element->trans_spec[5] ),(element->trans_spec[2] + element->trans_spec[3] + element->trans_spec[5] + element->trans_spec[4]),element->trans_spec[0],element->trans_spec[1],t);
      return(res);
    }
    else if( ( t >= (element->trans_spec[2] + element->trans_spec[3] + element->trans_spec[5] + element->trans_spec[4]) ) && ( t < element->trans_spec[2] + element->trans_spec[6])){
      return( element->trans_spec[0] );
    }
    else{return(-2);}
  }
  
  
  else if(element->flag == 3){
    if( (t >= 0) && (t < element->trans_spec[3]) ){
      res = element->trans_spec[0] + element->trans_spec[1]*sin( (2*pi*element->trans_spec[5])/360 );
      return(res);
    }
    else if( (t >= element->trans_spec[3]) && (t < finalTime) ){
      res = element->trans_spec[0] + element->trans_spec[1]*sin( 2*pi*( element->trans_spec[2]*(t-element->trans_spec[3]) + ( element->trans_spec[5]/360 ) ) )*exp(-(t-element->trans_spec[3])*element->trans_spec[4]);
      return(res);
    }
    else{ return(-3); }
  }
  
  
  else if(element->flag == 4){
    if( (t>=0) && (t<element->trans_spec[0]) ) { 
      return(element->trans_spec[1]); 
    }
    else if( (t>=element->trans_spec[2*PAIRS-2]) && (t<finalTime) ) { 
      return(element->trans_spec[2*PAIRS-1]); }
      else{
	for(i=0;i<2*PAIRS-3;i=i+2){
	  if( (t>=element->trans_spec[i]) && (t<element->trans_spec[i+2]) ){
	    if(element->trans_spec[i+1] < element->trans_spec[i+3]){
	      res = linear1(element->trans_spec[i],element->trans_spec[i+2],element->trans_spec[i+1],element->trans_spec[i+3],t);
	      return(res);
	    }
	    else{
	      res = linear2(element->trans_spec[i],element->trans_spec[i+2],element->trans_spec[i+1],element->trans_spec[i+3],t);
	      return(res);
	    }
	  }
	}
      }
  }
  
  else{ return(element->trans_spec[0]); }
  
  return(0);
}

/************************************************************************************/
double linear1(double t1, double t2, double val1, double val2,double t){
  
  
      double a,b,res;
  
  
      a = (val2-val1) / (t2-t1);
      b = val1 - a*t1;
       
      res = a*t + b;
      
      return(res);
      
      
}

/***********************************************************************************/
double linear2(double t1, double t2, double val1, double val2,double t){
  
  
      double a,b,res;
  
  
      a = (val1-val2) / (t1-t2);
      b = val2 + a*t1;
      
      res = -a*t + b;
      
      return(res);
      
      
}

/**********************************************************************************************************/
double *BE(cs *G,cs *C,double *x,double *b,double itol,int size,int plithos_nonz,double h){
  cs *A;
  double *newb;
  int i,j;
  
  
  A = cs_spalloc(size,size,plithos_nonz,1,1);
  newb = cs_malloc(size,sizeof(double));
  
  i = cs_gaxpy (C,x,newb);
  A = cs_add(G,C,1,1/h);
  
  for(i=0;i<size;i++){
    newb[i] = b[i] + 1/h*newb[i];
  }
  
  printf("LU\n");
  x = lu_sparse(A,plithos_nonz,size,b); 
  
  for(j=0;j<size;j++){
    printf("to teliko x(k+1) : x[%d] = %lf\n",j,x[j]);
  }
  printf("\n");
   
   return(x);
}

/*************************************************************************************************************/
double *TR(cs *G,cs *C,double *x,double *b,double *oldb,double itol,int size,int plithos_nonz,double h){
  cs *A,*B;
  double *newb,*d;
  int i,j;
  
  
 
  
  A = cs_spalloc(size,size,plithos_nonz,1,1);
  B = cs_spalloc(size,size,plithos_nonz,1,1);
  d = cs_malloc(size,sizeof(double));
  newb = cs_malloc(size,sizeof(double));
  
  A = cs_add(G,C,1,2/h);
  B = cs_add(G,C,1,-2/h);  
  
  newb = mul(B,x,size);
  for(i=0;i<size;i++){
    newb[i] = b[i] + oldb[i] - newb[i];
  }
  
  d = diag(A,size);
  
  x = lu_sparse(A,plithos_nonz,size,b); 
  
  
  for(j=0;j<size;j++){
    printf("to teliko x(k+1) : x[%d] = %lf\n",j,x[j]);
  }
  printf("\n");
   
  return(x);
  
}
/**********************EULER***********************************************/
double *back_euler(double **G,double **C,double *x1,double *b,int size,double h,double **A,double *tmp,double *da,double *y){
  
  int i,j;

  for(i=0;i<size;i++){
    for(j=0;j<size;j++){ 
      A[i][j] = G[i][j] + ((1/h)* C[i][j]); // A = (G+1/h *C)
    }
  }
      
    for(i=0;i<size;i++){   
	for(j=0;j<size;j++){
	      tmp[i] = tmp[i]+ C[i][j]*x1[j]; // C * x[tk-1]
	}
    
     }
	
 
  for(i=0;i<size;i++){
    da[i] = b[i] +((1/h) * tmp[i]); // da = b + 1/h ( C * x[tk-1] )
   
  }
   
     for(i=0;i<size;i++){   
	  for(j=0;j<size;j++){
	      G[i][j] = A[i][j];
	  }
	  
       }
       
	
      return(da);
  
  
}
/*****************************************************************************************************************************************/
double *trapezoidal(double **G,double **C,double *x1,double *b,double *b1,int size,double h,double **A,double **B,double *tmp,double *da,double *y){
  
 
  int i,j;
  
  
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){ 
      A[i][j] = G[i][j] + ((2/h)* C[i][j]);    // (G+2/h *C)
    }
  }
  
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){ 
      B[i][j] = G[i][j] - ((2/h)* C[i][j]);    // (G-2/h *C)
    }
  }
  
  
  for(i=0;i<size;i++){
  
    for(j=0;j<size;j++){
       tmp[i]= B[i][j] * x1[j]; // (G-2/h *C) * x[tk-1]
    }
   
  }

  
  for(i=0;i<size;i++){
    da[i] = b[i] + b1[i] -tmp[i]; // d = b + b1 (b[t-1]) -(G-2/h *C) * x[tk-1]
  }
  
       for(i=0;i<size;i++){   
	  for(j=0;j<size;j++){
	      G[i][j] = A[i][j];
	  }
	  
       }
      return(da);
  

}
/************************************************************************************/
void transient_analysis(double **pinakasG,double **pinakasC,int flag,double *x1,double *b,int size,double itol){

  
  double t,tempa,tempb,tempc;
  double **A,**B,**S,**pinakasP;
  double *da,*tmps,*y,*d,*k,*b1;
  double *data;
  int i,j,value;
  double end;
  double step;
  char type[2];
  
  struct node1 *tmp; 
  
    y = malloc(sizeof(double)*size);
    da  = calloc(sizeof(double)*size,sizeof(double)); 
    tmps = calloc(sizeof(double)*size,sizeof(double)); 
    d = malloc(sizeof(double)*size);
    b1 = calloc(sizeof(double)*size,sizeof(double));
    k = malloc(sizeof(double)*size);
    
    S = calloc(sizeof(double *)*size,sizeof(double*));
    for(i=0;i<size;i++){
      S[i]= calloc(sizeof(double)*size,sizeof(double));  
    }
    
    A = calloc(sizeof(double *)*size,sizeof(double*));
    for(i=0;i<size;i++){
      A[i]= calloc(sizeof(double)*size,sizeof(double));  
    }
    
    B = calloc(sizeof(double *)*size,sizeof(double*));
    for(i=0;i<size;i++){
      B[i]= calloc(sizeof(double)*size,sizeof(double));  
    }
    
    pinakasP = calloc(sizeof(double *)*size,sizeof(double*));
    for(i=0;i<size;i++){
      pinakasP[i]= calloc(sizeof(double)*size,sizeof(double));  
    }
    
    
    
    for(i=0;i<size;i++){
      for(j=0;j<size;j++){
	S[i][j]=pinakasG[i][j];
      }
    }
      
for(i = 0;i < size;i++){ k[i] = x1[i]; }
  
  FILE *fp3 =fopen("myNewFile.txt","w+");
 
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){  
    if(tmp->element->end != 0){
      end = tmp->element->end;
    }
    if(tmp->element->step != 0){
      step = tmp->element->step;
      type[0] = tmp->element->tplot[0];
      value = tmp->element->value_plot;
    }
    
  }

  
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){
    
    
    if( (tmp->element->type[0] == type[0]) && (atoi(&tmp->element->type[1]) == value) ){
      break;
    }
  }
    

    
    if(flag==7){fprintf(fp3,"KANW AC_SWEEP ME TH METHODO BE\n\n");}
    else if(flag==8){fprintf(fp3,"KANW AC_SWEEP ME TH METHODO TR\n\n");}

  
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){tempb=b[tmp->element->pos];}
      if(tmp->element->pos2 != -1){tempc=b[tmp->element->pos2];}
    }
   

    if(tmp->element->type[0] == 'V'){tempa=b[tmp->element->pos];}
    	for(i=0;i<size;i++){
	  if(pinakasG[i][i]!= 0){
		  d[i] = 1/pinakasG[i][i];
	  }else{
		  d[i] = 1.0;
	  }
	}
    
    for(t=0;t<=end;t=t+step){
      
      for(i=0;i<size;i++){
	for(j=0;j<size;j++){
	    pinakasP[i][j] = S[i][j];
	}
      }
      
      
      if(tmp->element->type[0] == 'V'){
	b[tmp->element->pos]=value_for_time(tmp->element,t,end); 
      }
      
      if(tmp->element->type[0] == 'I'){
	if(tmp->element->pos != -1){
 	  b[tmp->element->pos]= value_for_time(tmp->element,t,end);
	  
	  }
	  if(tmp->element->pos2 != -1){
	    b[tmp->element->pos2]=value_for_time(tmp->element,t,end);
	  }
	
      }
      	

      
      if(flag==8){
	
	printf("************** Back Euler! ********************************\n");
	
	data = back_euler(pinakasP,pinakasC,k,b,size,step,A,tmps,da,y);
	
      }
      else {				//TR
	  
	  printf("************** Trapezoidal! ********************************\n");
	  data = trapezoidal(pinakasP,pinakasC,k,b,b1,size,step,A,B,tmps,da,y);
	  
    }
      for(i=0;i<size;i++){
	    b1[i] = b[i];
	  }
	  
        k = CG(k,pinakasP,data,itol,size,d);
      
    fprintf(fp3,"Gia vhma %lf exoume %c[%d] = %lf\n",t,tmp->element->type[0],value,k[value-1]);

      if(tmp->element->type[0] == 'V'){b[tmp->element->pos]=tempa;}
      
      if(tmp->element->type[0] == 'I'){
	if(tmp->element->pos != -1){ b[tmp->element->pos]= tempb;}
	if(tmp->element->pos2 != -1){b[tmp->element->pos2]= tempc;}
     }
      
     
    }  
  
  
  fclose(fp3);
 }
 
 /***************************************************************************************************/
 
 void transient_sparse(cs *pinakasG,cs *pinakasC,int non_z,double *k,int size,double *b,int flag){
  
  double j;
  int i;//k;
  int counter=0;
  double itol = 1e-3;
  double tempb,tempc,tempa;
  double *res;
  double *d,*oldb;
  double *x1,*data;
  double *temp;
  double step,end;
  double type[2];
  double value;
  cs *A,*P;
  
  struct node1 *tmp;
  
  d = malloc(size*sizeof(double));
  oldb = malloc(sizeof(double)*size);
  x1 = malloc(non_z*sizeof(double));
  res = cs_malloc(size,sizeof(double));
  A = cs_spalloc(size,size,non_z,1,1);
  P = cs_spalloc(size,size,non_z,1,1);
  data = malloc(non_z*sizeof(double));
  
  FILE *fp2 =fopen("myfile.txt","w+");
  
  if(fp2 == NULL) {printf("error opening file\n");}
  
  
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){  
    if(tmp->element->end != 0){
      end = tmp->element->end;
    }
    if(tmp->element->step != 0){
      step = tmp->element->step;
      type[0] = tmp->element->tplot[0];
      value = tmp->element->value_plot;
    }
    
  }
  
 
  
  for(tmp = dc_root; tmp!= (struct node1 *) 0; tmp=tmp->nxt){  
    counter++;
    
    
    if(flag == 7){ 
      fprintf(fp2,"KANW TR_SWEEP SPARSE ME TRAPEZOIDAL\n\n"); 
    }
    else if(flag==8){
      fprintf(fp2,"KANW TR_SWEEP SPARSE ME BE\n\n"); 
    }
     
      d = diag(pinakasG,size);
      for(i=0;i<size;i++){ x1[i] = 0.0; }
   
    
    
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){
	tempb=b[tmp->element->pos];
      }
      if(tmp->element->pos2 != -1){
	tempc=b[tmp->element->pos2];
      }
    }
    if(tmp->element->type[0] == 'V'){
      if(tmp->element->pos != -1){
	tempa=b[tmp->element->pos];
      } 
    }
         
      A = pinakasG;
       temp = b;
      
      for(j=0;j<=end;j=j+step){
      
      res = temp;
      
      b = res;
	
	  P = A;
          
      if(tmp->element->type[0] == 'I'){
	if(tmp->element->pos != -1){
	  b[tmp->element->pos]= value_for_time(tmp->element,j,end);
	  if(tmp->element->pos2 != -1){
	    b[tmp->element->pos2]=value_for_time(tmp->element,j,end);
	  }
	}
      }
      
      if(tmp->element->type[0] == 'V'){
	if(tmp->element->pos != -1){
	  b[tmp->element->pos]=value_for_time(tmp->element,j,end);
	}
      }
      
      
	
	if(flag == 7){
	  fprintf(fp2,"*************************** SPARSE BE ******************************\n");
	  data = TR(P,pinakasC,x1,b,oldb,itol,size,non_z,step);
	  for(i=0;i<size;i++){
	    b[i] = oldb[i];
	  }
	}
	else if(flag == 8){
	  fprintf(fp2,"*************************** SPARSE TRAPEZOIDAL ******************************\n");
	  data=BE(P,pinakasC,x1,b,itol,size,non_z,step);
	}
      
	 x1 = sparse_CG(x1,P,data,itol,size,d);
      
      x1=res;
      
      for(i=0;i<size;i++){
	fprintf(fp2,"x[%d] = %lf \n",i,res[i]);
      }
      fprintf(fp2,"\n");
        
      
    }
    if(tmp->element->type[0] == 'I'){
      if(tmp->element->pos != -1){ b[tmp->element->pos]= tempb;}
      if(tmp->element->pos2 != -1){b[tmp->element->pos2]= tempc;}
  }
  if(tmp->element->type[0] == 'V'){
    if(tmp->element->pos != -1){
      b[tmp->element->pos]=tempa;
    }
}

}
fclose(fp2);
}