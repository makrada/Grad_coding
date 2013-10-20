/*
Created on 15/07/2012
Author : Marieta Kranta 
*/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "../../INCLUDE/const.h"
# include "../../INCLUDE/struct.h"
# include "../../INCLUDE/var_ext.h"

static FILE *pin_prs = NULL;
static char f[MAXSTRING];


extern int parser() ;
extern FILE* efopen();
extern FILE *jjin,*jjout,*mfp,*nfp2 ;


int **ret_pin = NULL;//is a two-dimensional matrix that stores the information provided by .pin =>ret_pin[#p-inv][#place_integer]
int ln_out = 0;// declares the number of the p-invariants

/*parsing the .pin file which contains the information that we are interested in*/
int read_PIN_file(int npl){
   
	char *buffer, *num_token;
	int  i, j, ch, k, cnt, total_size, fill;
	int *tmp, *num_pinv_buffer2;
	int N,init_N;
	int *famous_pl; // keeps in how many pinvs each place is involved
	
	buffer = NULL;
	num_token = NULL;
	famous_pl = NULL;
	i=0;
	j=0;
	ch=0;
	k=0;
	cnt=0; 
	fill =1; 
	total_size=0;
	tmp = NULL;
	num_pinv_buffer2 = NULL; // contains the number of places that are involved in each P-inv
	N = 100; // equals to Max_Token_Bound
	init_N = N;
	
	
	sprintf(f,"%spin",net_name); 
	
	
	tmp =(int *)malloc(N*sizeof(int));
	buffer = (char *)malloc(N*sizeof(char)); 
	num_pinv_buffer2 = (int *)malloc(N*sizeof(int));
	famous_pl = (int *)calloc(sizeof(int)*npl,sizeof(int));
	
	
	if(tmp == NULL){
		fprintf(stderr,"Error in malloc in tmp (read_PIN_file) \n");
	}
	
	if(buffer == NULL){
		fprintf(stderr,"Error in malloc in buffer (read_PIN_file) \n");
	}
	
	if(num_pinv_buffer2 == NULL){
		fprintf(stderr,"Error in malloc in buffer2 (read_PIN_file) \n");
	}
	
	if(famous_pl == NULL){
		fprintf(stderr,"Error in malloc in famous_pl (read_PIN_file) \n");
	}
	
	pin_prs = efopen(f,"r");
	
	if( pin_prs == NULL ){ 
		fprintf(stderr,"Can't open file %s r\n",f);
		 return -1;
	}
      
	
	while (feof(pin_prs) == 0) { //  parsing the file (.pin) put -1 in the end of each line
		fscanf(pin_prs, "%[^\n]\n", buffer);
		
				
		num_token = strtok (buffer," "); // tokenize the line
		
		while (num_token != NULL ){ 
	  
			sscanf(num_token, "%d", &ch);
			tmp[i] =(char) ch;
			num_token = strtok (NULL, " ");
			i++;
			
			if(i==N){
				N += N;
				tmp =realloc(tmp,N*sizeof(int));
			  
				if(tmp == NULL){
					fprintf(stderr,"Error in malloc in tmp (read_PIN_file) \n");
				}
			}
		}
		tmp[i] = -1; // denotes end of line
		i++;
		
	} 


	total_size = i;
	
	if(total_size > init_N){
		num_pinv_buffer2 =realloc(num_pinv_buffer2,(total_size*sizeof(int)));
	}
	
	
	for(i=0; i<total_size;i++){ // i need the index where the line changes, i.e. where each new pinv begins
		if(tmp[i]==-1){
			if(tmp[i+1] != 0){
			num_pinv_buffer2[cnt] = tmp[i+1]; // buffer2 contains the number of pinv in each line
				cnt++;
			}
		}
	}
		
	

	ln_out = tmp[0]; // number of p-inv


	// ret_pin[#p-inv][#places]
	ret_pin =(int **)calloc(sizeof(int *)*tmp[0],sizeof(int*)); 
	
	for(i=0;i<tmp[0];i++){
		ret_pin[i]=(int *)calloc(sizeof(int)*npl,sizeof(int));  
	}
/*tmp vector contains the integer which is equal to index of the place in the tool*/

         k=2; // i need to store only the places in the ret_pin
	for(j=0; j<tmp[0]; j++){ // 
	  k += 2; // skip #places , #-1, first_token
		for(i=0; i<num_pinv_buffer2[j] ;i++){
			
			ret_pin[j][tmp[k]-1] = 1; // -1 because of the vector [0...N-1]  
			k += 2;
		}
	}
	
	
// how famous is a place, for future use

	for(i=0; i<npl; i++){
		for(j=0; j<tmp[0]; j++){
			if(ret_pin[j][i] == 1){
				famous_pl[i] += 1;
			}
		}
	}
		
	
	
	fclose(pin_prs);
	free(tmp);
	free(buffer);
	free(num_pinv_buffer2);

	return 0;
}

/********************************************* NEW VERSION *******************************************************/
/*We sort the places according to the new proposed heuristic algorithm*/

int *sort_according_to_pinv(int num_pinv, int npl){

	 int min, max,order,back, ok;
	 int cnt, keep, keep_sec, map_position; //  keep = position on min, keep_sec = position of max
	 int i, j, k, tmp, itr, run; 
	 int visited_pinv, how_many, move ;
	 int  *common_pl, *new_map_sort;
	 int  *not_touched, *checked_pinv;
	 int *pntr_ret;
	 
	 
	 min = 100;
	 max = order = 0;
	 cnt = back = 0;
	 i=j=k = 0;
	 keep = 0;
	 keep_sec =0;
	 tmp =0;
	 map_position =0;
	 common_pl =NULL;
	 checked_pinv = NULL;
	 not_touched = NULL;
	 new_map_sort = NULL;
	 itr =0;
	 run =0;
	 visited_pinv=0;
	 move = 1;
	 pntr_ret = NULL;

	 common_pl = (int *)calloc(npl*sizeof(int),sizeof(int)); 
	 
	  if(common_pl == NULL){
		fprintf(stderr,"Error in malloc in common_pl (read_PIN_file) \n");
	  }
	 
	 new_map_sort = (int *)calloc(npl*sizeof(int),sizeof(int));
	 
	  if(new_map_sort == NULL){
		fprintf(stderr,"Error in malloc in new_map_sort (read_PIN_file) \n");
	  }
	  
	  not_touched = (int *)calloc(npl*sizeof(int),sizeof(int));
	  
	    if(not_touched == NULL){
		fprintf(stderr,"Error in malloc in not_touched (read_PIN_file) \n");
	    }
	  
	  checked_pinv = (int *)calloc(num_pinv*sizeof(int),sizeof(int));
	  
	    if(checked_pinv == NULL){
		fprintf(stderr,"Error in malloc in checked_pinv (read_PIN_file) \n");
	    }
	  
	  
	  
	    for(i=0; i<num_pinv; i++){
		checked_pinv[i] =-1;
	    }
	  
	  
	  
	    for(i=0; i<npl; i++){
		 not_touched[i] = -1;
		 new_map_sort[i] = -1;
	    }
	    
	pntr_ret = new_map_sort;
	 

	 
	keep = find_min(ret_pin,num_pinv,npl); 
	
	checked_pinv[keep] = keep;//  keep the first p-inv i keep OLNY the pinvs that I have dealt with and are not needed any more
		  cnt=0;
		  
	while(move){
	
	  
	 keep_sec = find_max(ret_pin, keep, num_pinv, npl,checked_pinv);
	 
	 
	if(keep_sec == -1){ //if intersection is empty, put in the vector only the places of the first p-invariant
		for(i=0; i<npl; i++){ 
		  
			if( (ret_pin[keep][i] == 1) ){
			
				itr = find_if_exists_place(new_map_sort,npl,i); // if ith place already exists i do nothing
		      
				if(!itr){ // else i store it in the vector
				      new_map_sort[cnt] = i;
				      cnt++;
				}
			}		  
		}
		
		
		
	}
	 
	else{  // if intersection is non empty
		for(i=0; i<npl; i++){ // store non common places
		  
		      if( (ret_pin[keep][i] == 1) && !(ret_pin[keep][i] & ret_pin[keep_sec][i]) ){
			
				itr = find_if_exists_place(new_map_sort,npl,i); // if the place already exists 
		      
				if(!itr){ 
				      new_map_sort[cnt] = i;
				      cnt++;
				}
		      }		  
		}
		  
		  
		for(i=0; i<npl; i++){ // store only common places between p-invariant
		
			if( (ret_pin[keep][i] & ret_pin[keep_sec][i]) ){
			      
			        back = find_if_exists_place(not_touched,npl,i);
			        
		             
				if(!back){ // if it isn't already in common places
					order = find_if_exists_place(new_map_sort,npl,i);
			       
					if(!order){ // if it isn't on the new_map_sort
					      new_map_sort[cnt] = i; // add new place
					      not_touched[i] = i; // places which I am not allowed to touch
					      cnt++;
					}
					else{ // if common places exists and i can touch it (isn't already in common )
				  
						for(j=0; j<npl; j++){
							if(new_map_sort[j] == i){ // keep the position of new common place
							break;
							}
						}

						for(tmp=j; tmp<cnt-1; tmp++){
							new_map_sort[tmp] = new_map_sort[tmp+1]; // reorder the vector
						}
					
					
						new_map_sort[cnt-1] = i;
					
					
						for(k=0; k<cnt; k++){
							for(j=k+1;j<cnt;j++){
								if(new_map_sort[k] == new_map_sort[j] ){
									fprintf(stdout,"problem\n");
									exit(0);
								}	
							}
						}
					}
				}
			}
		
		}  
		      
		
		
		for(i=0; i<npl; i++){ // store places of the second p-invariant
		
			if(ret_pin[keep_sec][i] == 1 && !(ret_pin[keep_sec][i] & ret_pin[keep][i])){
				
				itr = find_if_exists_place(new_map_sort,npl,i); 
		      
				if(!itr){ 
				      new_map_sort[cnt] = i;
				      cnt++;
				}
			}
		
		
		}
	}
			
		move =0;
			for(i=0;i<num_pinv; i++){
				if(checked_pinv[i] != -1 || cnt == npl){
					continue;
				}			
				else{
					checked_pinv[i] = i;
					keep = i; // next p-invariant to deal with
					move=1;
					break;
				}
			}
		
		}
		
	
	free(common_pl);
	free(new_map_sort);
	free(not_touched);
	free(checked_pinv);

	return pntr_ret;

}

/****************************************** FIND MIN INTERSECTION ****************************************************/
/*
 Finds the P-invariant with the minimum intersection with the others P-invariants
 and returns its position on the vector where P-invariants are stored
 
 */
int find_min(int **ret_pin, int num_pinv, int npl){

  int i,j,cnt,k;
  int *inter_vec_min;
  int min,keep;
  
  i=j=cnt=k=0;
  keep=0;
  min=100;
  
  inter_vec_min = NULL;
  
  
	  inter_vec_min = (int *)calloc(num_pinv*sizeof(int),sizeof(int));
	 
	  if(inter_vec_min == NULL){
		fprintf(stderr,"Error in malloc in inter_vec_min (read_PIN_file) \n");
	  }
	  
	 
	for(i=0; i< num_pinv; i++){
		for(j=i+1; j<num_pinv; j++){
		  k=0;
			 for(cnt=0; cnt< npl; cnt++){ 
				if(ret_pin[i][cnt] & ret_pin[j][cnt]){
					k++; // i = Pprev
				} 
			 }
			 if(k > 0){
				inter_vec_min[i]++;
			 }
			 
			if(inter_vec_min[i] < min && k>0){
				min = inter_vec_min[i]; // i shows the pinv with the minimum number of intersection with othe pinv
				keep = i;
			}
			if( k==0 ){ // if there is no intersection with the other P-invariants
				min = 0; 
				keep = i;
			}
		}
		
	}  



      free(inter_vec_min);
      return keep;

}


/******************************************** FIND MAX INTERSECTION  **************************************************/
/*

  Finds the P-invariant with the maximum intersection with the others P-invariants
 and returns its position on the vector where P-invariants are stored
 
 */


int find_max(int **ret_pin, int keep, int num_pinv, int npl, int *checked_pinv){

  int i,j,cnt,max,keep_sec,k;
  int *num_inter_max;
  
  i=j=keep_sec=cnt=k=0;
  max = -1;
  num_inter_max = NULL;
  
  
   num_inter_max =  (int *)calloc(num_pinv*sizeof(int),sizeof(int));
	 
	  if(num_inter_max == NULL){
		fprintf(stderr,"Error in malloc in num_inter_max (read_PIN_file) \n");
	  }
  

 for(i=0; i<num_pinv; i++){
			if(i == checked_pinv[i]){ // if p-invariant is already checked_pinv  
				continue;
			}
			else{ 
					k=0;
						for(cnt=0; cnt< npl; cnt++){
							if(ret_pin[keep][cnt] & ret_pin[i][cnt]){
								  k++; // i = Pprev
								  num_inter_max[i]++; // max inter #places
							} 
						}
			 
						if(num_inter_max[i] > max && k>0){
							max = num_inter_max[i];// "i" shows the pinv with the maximum number of intersection with other pinv
							keep_sec = i;
						}
						
						if(k==0){
							if(max == -1){
								keep_sec = -1; // if intersection is empty we deal only with the first p-invariant
							}
						}
			}
		
		}
		
	free(num_inter_max);
	return keep_sec;

}
 /***************************************** IF THE PLACE ALREADY EXIST DO NOTHING **************************************/
/*
 Checks if one place already exists in the new sorting vector or in the common places vector
 */
int find_if_exists_place(int *exist_in_vector, int npl, int place){

  int i;

	  for(i=0;i<npl;i++){
			if(( exist_in_vector[i] == place ) ){
				return 1;	// if exists in either new sorting vector or in common places vector		
			}		
	  }

      return 0;

}
