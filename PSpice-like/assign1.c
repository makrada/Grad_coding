  /*
   * authors :Kranta Marieta  &&  Vouronikou Mina  &&  Konstantinidou Theodosia
   *                    
   */
  
  #include "types.h"
  #include <stdio.h>
  #include <string.h>
  #include <stdlib.h>
  #include "list.h"
  #include "functions.h"
  #include <ctype.h>
  #include <math.h>
  #include "cs.h"
  
  #define PAIRS 5
  
  int main(int argc, char **argv){
    
    
    char str[20],str1[20],str2[3];
    int m2,size,i,size_of_str,j;
    double **pinakasG,**pinakasC,**A,**B;       	
    double *b,*b1,*tmp,*da;              		
    double *k, *x,*y,*x1,*y1,*d,*d1,*data;   
    double **S;
    char **vector;         				
    char c;
    int *p,*p1;
    int flag,flag1,flag2;
    double itol = 1e-3;
    double term1,term2,term3,term4,var;     		 
    int max;     					
    int cnt =0;
    FILE *fp;
    struct elem *element;
    struct dc_nodes *elem1;
    struct  dc_nodes *element1;

    /*Initialization */
    flag=0;
    flag1=2;
    max=0;
    m2=0;   			
    size = 0;
    i=0;
    
    //open file to read netlist
    fp = fopen(argv[1],"r");
    
    
    if (fp == NULL){
      printf("Cannot open file\n");
      return -1;
      
    }
    
    do{
      
      if (fscanf(fp,"%s",str) == EOF) break;
      
      size_of_str = strlen(str);
      size_of_str = size_of_str+1; 
      
      
      if (str[0] == '*'){
	do {
	  c=getc(fp);
	} while (c != '\n');
	continue;
      } 
      
      element = (struct elem *)malloc(sizeof(struct elem));
	
                                                
      if((str[0] =='V')||(str[0] =='v')||(str[0] =='I')||(str[0] =='i')||(str[0] =='R')||(str[0] =='r')||(str[0] =='L')||(str[0] =='l')||(str[0] =='C')||(str[0] =='c')){
	if((str[0] =='V')||(str[0] =='v')) m2++;
	str[0] = toupper(str[0]);
	element->type = malloc(sizeof(char)*size_of_str);
	strcpy(element->type,str);
	element->term = malloc(sizeof(double)*2);
	element->value = malloc(sizeof(double));
	fscanf(fp,"%lf",&term1);
	
	if(term1 >=0){
	  element->term[0] = term1; 
	  
	}
	else printf("Terms must be non negative!\n");
	fscanf(fp,"%lf",&term2);
	
	if(term2 >=0){
	  element->term[1] = term2;
	  
	}
	else printf("Terms must be non negative!\n");
	
	fscanf(fp,"%lf",&element->value[0]);
	
	if((str[0] =='V')||(str[0] =='v')||(str[0] =='I')||(str[0] =='i')){
	  
	 element1= (struct dc_nodes *)malloc(sizeof(struct dc_nodes));
	 element1->type = malloc(sizeof(char)*size_of_str);
	 strcpy(element1->type,str);
	 element1->tplot = malloc(sizeof(char )*5);
	  c = getc(fp);
	  if(c != '\n'){
	    fscanf(fp,"%s",str1);
	    if(!strcmp(str1,"EXP")){				
		fscanf(fp,"%s",str1);
		
		element1->trans_spec = malloc(sizeof(double)*6);
		for(i=0;i<6;i++){ fscanf(fp,"%lf",&element1->trans_spec[i]); }
		
		fscanf(fp,"%s",str1);
		element1->flag = 1;
	    }
	    else if(!strcmp(str1,"PULSE")){		
		fscanf(fp,"%s",str1);
		
		element1->trans_spec = malloc(sizeof(double)*7);
		for(i=0;i<7;i++){ fscanf(fp,"%lf",&element1->trans_spec[i]); }
		
		fscanf(fp,"%s",str1);
		
		element1->flag = 2;
	    }
	    else if(!strcmp(str1,"SIN")){		
		fscanf(fp,"%s",str1);
		
		element1->trans_spec = malloc(sizeof(double)*6);
		for(i=0;i<6;i++){ fscanf(fp,"%lf",&element1->trans_spec[i]); }
		
		fscanf(fp,"%s",str1);
		
		element1->flag = 3;
	    }
	    else if(!strcmp(str1,"PWL")){		
		element1->trans_spec = malloc(sizeof(double)*(2*PAIRS));
		
		while(getc(fp)!='\n'){
		  fscanf(fp,"%s",str1);
		  
		  fscanf(fp,"%lf",&var);
		  element1->trans_spec[cnt] = var;
		  cnt++;
		  fscanf(fp,"%lf",&var);
		  element1->trans_spec[cnt] = var;
		  cnt++;
		  
		  fscanf(fp,"%s",str1);
		}
		
		element1->flag = 4;
	    }
	  }
	  else if(c == '\n')
	  {				
	      element1->trans_spec = malloc(sizeof(double));
	      element1->trans_spec[0] = element->value[0];
	      
	      element1->flag = 5;
	  }
	  
	  element1->pos = -1;
	  element1->pos2 = -1;
	  
	   dc_add_list(element1);
         
	}

      }else if((str[0] =='Q')||(str[0] =='q')){     
	str[0] = toupper(str[0]);
	element->type = malloc(sizeof(char)*size_of_str);
	strcpy(element->type,str);
	element->term = malloc(sizeof(double)*3);
	
	fscanf(fp,"%lf",&term1);
	if(term1 >=0) element->term[0] = term1; 
	else printf("Terms must be non negative!\n");
	
	fscanf(fp,"%lf",&term2);
	if(term2 >=0) element->term[1] = term2;
	else printf("Terms must be non negative!\n");
	
	fscanf(fp,"%lf",&term3);
	if(term3 >=0) element->term[2] = term3;
	else printf("Terms must be non negative!\n");
	
	
      }else if((str[0] =='M')||( str[0] =='m')){                     
	str[0] = toupper(str[0]);
	element->type = malloc(sizeof(char)*size_of_str);
	strcpy(element->type,str);
	element->term = malloc(sizeof(double)*4);
	
	fscanf(fp,"%lf",&term1);
	if(term1 >=0) element->term[0] = term1; 
	else printf("Terms must be non negative!\n");
	
	fscanf(fp,"%lf",&term2);
	if(term2 >=0) element->term[1] = term2;
	else printf("Terms must be non negative!\n");
		    
		    fscanf(fp,"%lf",&term3);
	if(term3 >=0) element->term[2] = term3;
		    else printf("Terms must be non negative!\n");
		    
		    fscanf(fp,"%lf",&term4);
	if(term4 >=0) element->term[3] = term4;
	else printf("Terms must be non negative!\n");
	
	
	element->value = malloc(sizeof(double)*2);
	fscanf(fp,"%lf",&element->value[0]);
	fscanf(fp,"%lf",&element->value[1]);
	
      }else if((str[0] =='D')||(str[0] =='d')){ 				    
	str[0] = toupper(str[0]);
	element->type = malloc(sizeof(char)*size_of_str);
	strcpy(element->type,str);
	element->term = malloc(sizeof(double)*2);
	
	fscanf(fp,"%lf",&term1);
	if(term1 >=0) element->term[0] = term1; 
	else printf("Terms must be non negative!\n");
		    
		    fscanf(fp,"%lf",&term2);
	if (term2 >=0) element->term[1] = term2;
		    else printf("Terms must be non negative!\n");
		    
		    
		    
      }else if(str[0] == '.'){
	
	i=1;break;
      }
      
      
      // add elements
      add_list(element);
      
    }while(!feof(fp));
    
    if(i==1){
      
      do{
	if(!strcmp(str,".OPTIONS")) {
	  fscanf(fp,"%s",str);
	  
	  
	  
	  if(!strcmp(str,"SPARSE")){	    
	    flag2 = 5; 				//SPARSE
	    fscanf(fp,"%s",str);
	      
	  }
	  
	  if(!strcmp(str,"SPD")) {		
	    flag=1;				//Cholesky
	    fscanf(fp,"%s",str);
	    
	  }
	  else if(!strcmp(str,"iter")) {
	    fscanf(fp,"%s",str);
	    if(!strcmp(str,"SPD")) { 
	      flag =3;						// CG 
	      fscanf(fp,"%s",str);
	    }
	     else { 
	      flag = 4; 				 	//Bi - CG
	      fscanf(fp,"%s",str);	      
	    }  
	  }
	  else  if(!strcmp(str,"METHOD=TR")){
	    flag = 7;
	    
	  }
	  else if(!strcmp(str,"METHOD=BE")){
	    flag = 8;  
	  }
	  else{
	    flag=2;						// LU
	    fscanf(fp,"%s",str);
	    
	  }	  
	}
	if(!strcmp(str,".DC")) {
	  elem1 = (struct dc_nodes *)malloc(sizeof(struct dc_nodes));
	  flag1=1;
	  fscanf(fp,"%s",str1);
	  if((str1[0] =='V')||(str1[0] =='v')||(str1[0] =='I')||(str1[0] =='i')){
	    str1[0] = toupper(str1[0]);
	    elem1->type = malloc(sizeof(char)*size_of_str);
	    strcpy(elem1->type,str1);
	    fscanf(fp,"%lf",&elem1->start); 
	    fscanf(fp,"%lf",&elem1->end); 
	    fscanf(fp,"%lf",&elem1->step);
	    fscanf(fp,"%s",str2);
	    
	    if(!strcmp(str2,".PLOT")){
	      fscanf(fp,"%s",str2);
	      str2[0] = toupper(str2[0]);
	      elem1->tplot = malloc(sizeof(char)*size_of_str);
	      strcpy(elem1->tplot,str2);
	      fscanf(fp,"%d",&elem1->value_plot); 
	    }
	  }
	 
	  elem1->pos = -1;
	  elem1->pos2 = -1;
	  dc_add_list(elem1);
	  
	} else if(!strcmp(str,".TRAN")){
	    printf("mpikes?\n");
	    flag1=2;
	      elem1 = (struct dc_nodes *)malloc(sizeof(struct dc_nodes));
	      elem1->type = malloc(sizeof(char)*size_of_str);
	      strcpy(elem1->type,str);
	      elem1->tplot = malloc(sizeof(char )*10);
	      elem1->trans_spec = malloc(sizeof(double)*2*PAIRS);

	     fscanf(fp,"%lf",&elem1->step);
	     fscanf(fp,"%lf",&elem1->end); 
	     printf("step : %lf\n",elem1->step);
	     printf("end : %lf\n",elem1->end);
	     fscanf(fp,"%s",str2);
	    
	    if(!strcmp(str2,".PLOT")){
	      fscanf(fp,"%s",str2);
	      str2[0] = toupper(str2[0]);
	      elem1->tplot = malloc(sizeof(char)*size_of_str);
	      strcpy(elem1->tplot,str2);
	      fscanf(fp,"%d",&elem1->value_plot); 
	    }
	       elem1->pos = -1;
	      elem1->pos2 = -1;
	     dc_add_list(elem1);
	  }
	  
	
	fscanf(fp,"%s",str);
	
      }while(!feof(fp));
    }
    
    
    m2 = change_nodes(m2);
    
    max = mapa();
    max = max-1;
    size = max+m2;
    
    
    pinakasG = calloc(sizeof(double *)*size,sizeof(double));
    for(i=0;i<size;i++){
      pinakasG[i]= calloc(sizeof(double)*size,sizeof(double));  
    }
    
    pinakasC = calloc(sizeof(double *)*size,sizeof(double));
    for(i=0;i<size;i++){
      pinakasC[i]= calloc(sizeof(double)*size,sizeof(double));  
    }
    
    vector = calloc(sizeof(char *)*size,sizeof(char*));
    for(i=0;i<size;i++){
      vector[i]= calloc(sizeof(char)*100,sizeof(char));  
    }
    
    b = calloc(sizeof(double)*size,sizeof(double));
    b1 = calloc(sizeof(double)*size,sizeof(double));
    y = malloc(sizeof(double)*size);
    x = malloc(sizeof(double)*size);
    y1 = malloc(sizeof(double)*size);
    x1 = malloc(sizeof(double)*size); 
    d1 = malloc(sizeof(double)*size);
    d = malloc(sizeof(double)*size);
    da  = calloc(sizeof(double)*size,sizeof(double)); 
    tmp = calloc(sizeof(double)*size,sizeof(double)); 
    p1 = malloc(sizeof(int)*size);
    
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
    
    
    if(flag2!=5){
      
      
      printf("\n************************************     MNAC   ****************************************\n\n");
      MNAC(pinakasC,size,b,vector,max,m2,elem1);  
      
      printf("\n************************************     MNAG   ****************************************\n\n");
      change_name(m2);
      MNAG(pinakasG,size,b,vector,max,m2,elem1);
      
      
      if((flag == 1) || (flag == 2 )){
	
	if(flag == 1){							 //  .OPTIONS SPD  Cholesky
    
      printf("\n************************************ Cholesky decomposition ****************************************\n\n");
      Cholesky(pinakasG,size);
	}
	
	
	else{    							//  .OPTIONS  LU
      
      printf("\n************************************ LU decomposition **********************************************\n\n");
      p = malloc(sizeof(int)*size);
      data = malloc(sizeof(double)*size);
      p = LU(pinakasG,size);
      
      printf("\n To p einai :\n");
      for(i=0;i<size;i++){
	printf(" to [%d]  -- > [%d]\n",i,p[i]);
      } 
      
      b= reverse_b(p,b,size);
      printf("\nTo anestrammeno b einai :\n"); 
      
      for(i=0;i<size;i++){
	printf(" b[%d] = %lf\n",i,b[i]);
      }
      printf("\n\n"); 
      
	}
	y = forw_sub(pinakasG,b,size,y,flag);
	printf("\n"); 
	
	for(i=0;i<size;i++){
	  printf(" y[%d] = %lf\n",i,y[i]);
	}
	printf("\n"); 
	x= back_sub(pinakasG,y,size,x,flag);
	
	for(i=0;i<size;i++){
	  printf(" x[%d] = %lf \n",i,x[i]);
	} 
	
	
	
      }
      
      else if ((flag ==3) || (flag == 4)){
	
	k=malloc(sizeof(double)*size);
	
	for(i=0;i<size;i++){
	  
	  k[i] = 0.0;
	}
	d = malloc(sizeof(double)*size);
	for(i=0;i<size;i++){
	  if(pinakasG[i][i]!= 0){
	    d[i] = 1/pinakasG[i][i];
	  }else{
	    d[i] = 1.0;
	  }
	}
	
	
	
	if (flag == 3){
	  printf("\n************************************ CG decomposition ****************************************\n\n");
	  CG(k,pinakasG,b,itol,size,d);  
	  
	  
	}
	else{ 
	  
	  d1 = malloc(sizeof(double)*size);
	  S = calloc(sizeof(double)*size,sizeof(double));
	  for(i=0;i<size;i++){
	    S[i]= calloc(sizeof(double)*size,sizeof(double));  
	  }
	  for(i=0;i<size;i++){
	    for(j=0;j<size;j++){
	      S[j][i] = pinakasG[i][j];
	    }
	  }
	  
	  for(i=0;i<size;i++){
	    d1[i] = 1/S[i][i];
	  }
	  printf("\n************************************ Bi-CG decomposition ****************************************\n\n");
	  Bi_CG(k,pinakasG,b,itol,size,d,d1,S);
	  
	  
	}
      }

      else if(flag == 7 || flag == 8){              //BE or TR
	
	printf("\n**CG decomposition **\n\n");
	
	k=malloc(sizeof(double)*size);
	
	for(i=0;i<size;i++){
	  
	  k[i] = 0.0;
	}
	
	d = malloc(sizeof(double)*size);
	
	for(i=0;i<size;i++){
	  if(pinakasG[i][i]!= 0){
	    d[i] = 1/pinakasG[i][i];
	  }else{
	    d[i] = 1.0;
	  }
	}
	
	CG(k,pinakasG,b,itol,size,d);
	
	if(flag == 8){			//BE
	  
	  
	  printf("************** Back Euler! ********************************\n");
	  data = back_euler(pinakasG,pinakasC,x,b,size,0.1,A,tmp,da,y);
	  
	  
	}
	else {				//TR
	  
	  printf("************** Trapezoidal! ********************************\n");
	  data = trapezoidal(pinakasG,pinakasC,x,b,b1,size,0.1,A,B,tmp,da,y);
	  
	}
	
	// allagi to b egine  data
	CG(k,pinakasG,data,itol,size,d);
	
	
      }
      
    

    }else if(flag2 == 5){
      
      printf("*****************************************************************************************************\n");
      printf("************************************ This is SPARSE Analysis ****************************************\n");
      printf("*****************************************************************************************************\n\n\n\n");
	
	s_p_a_r_s_e(max,m2,flag,flag1);
	
    }
    
    fclose(fp);	
    
    if(flag1==1 && flag2!=5){
      
      printf("************************************ DC - sweep analysis  ****************************************\n\n");
      
      
      dc_sweep(pinakasG,x1,y1,b,size,flag);
      
    }
    else if(flag1==2 && flag2!=5){
// 	

// 	//transient analysis
 	transient_analysis(pinakasG,pinakasC,flag,k,b,size,itol);

       }
    
   
  
    
    //checks if there are 2 elements with the same name
    check_list();
    
    
    clear_list();
    clear_list2();
    
    return 0;
    
  }
