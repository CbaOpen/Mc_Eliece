#include <flint/fmpz.h>
#include <stdio.h>
#include <flint/fmpz_poly.h>
#include "encryptRS.h"
#include "decryptRS.h"
#include <unistd.h>
#include <time.h>
#include <math.h>

#include <math.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_poly.h>
#include "../inc/encryptRS.h"
#include "../inc/decryptRS.h"
#include "../inc/main.h"	
#include <stdlib.h>
#include <time.h>


extern fmpz_mod_poly_t reduc;
extern fmpz_mod_poly_t* tab;
extern fmpz_t* ptoi;
extern fmpz_t DEUX;    


fmpz_mat_t P,Gm,S;  //à supp
fmpz_mod_poly_t G1;


				// s  === bit par symbole   en param ou a calculer -> = ln(n+1)/ln(2) === log2(n+1)  / 
int gf;	




														//row op1, col1=row2,col op2
void multGFmat(fmpz_mat_t res,fmpz_mat_t op1,fmpz_mat_t op2,int n,int p,int m){
	fmpz_t tmp,tmp2;
	fmpz_init(tmp);
	fmpz_init(tmp2);
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			//fmpz_set_ui(tmp,0);
			//printf("%d %d \n",i,j);
			fmpz_init_set_ui(fmpz_mat_entry(res, i, j),0);
			for(int k=0;k<p;k++){
					//fmpz_set_ui(tmp2,0);
					fmpz_t o1,o2;
					fmpz_init(o1);
					fmpz_init(o2);
					fmpz_set(o1,fmpz_mat_entry(op1, i, k));
					fmpz_set(o1,ptoi[fmpz_get_ui(o1)]);    //alplha op1

					fmpz_set(o2,fmpz_mat_entry(op2, k, j));
					fmpz_set(o2,ptoi[fmpz_get_ui(o2)]);
					if(fmpz_is_zero(o1) || fmpz_is_zero(o2)){
						//printf("0 %d %d %d\n",i,j,k);
						fmpz_set_ui(o1,0);

					}
					else{
						fmpz_add(o1,o1,o2);
						fmpz_mod_ui(o1,o1,gf-1);
						tab[fmpz_get_ui(o1)]->p=gf;
						fmpz_set(o2,o1);
						fmpz_mod_poly_evaluate_fmpz(o1,tab[fmpz_get_ui(o1)],DEUX);
						//tab[fmpz_get_ui(o2)]->p=2;

					}
					fmpz_xor(fmpz_mat_entry(res, i, j),fmpz_mat_entry(res, i, j),o1);






					//fmpz_mul(tmp2,fmpz_mat_entry(op1, i, k),fmpz_mat_entry(op2, k, j));
					//fmpz_add(fmpz_mat_entry(res, i, j),tmp2);
					//fmpz_set(fmpz_mat_entry(res, i, j),tmp);
			}
		}
	}
}
void RPermut(fmpz_mat_t P,int n){
int tabP[n]	;
long int j;
int temp;
srandom(getpid()+time(NULL));
	for(int i=0;i<n;i++){
		tabP[i]=i;
	}
	
	for(int i=n-1;i>0;i--){
		j = random()%i;
		temp=tabP[i];
		tabP[i]=tabP[j];
		tabP[j]=temp;

	}

	for(int i=0;i<n;i++){
		fmpz_set_ui(fmpz_mat_entry(P, i, tabP[i]), 1);
	}

}



void get_S(fmpz_mat_t S,int k){
	
	
	flint_rand_t rand;
	flint_randinit(rand);
	fmpz_t det;
	fmpz_init(det);
	fmpz_mat_randrank(S , rand ,k , 1);
	fmpz_mat_det(det , S);
	while(!fmpz_get_ui(det)){
		fmpz_mat_randrank(S , rand ,k , 1);
		fmpz_mat_det(det , S);

	}
}


void get_G(fmpz_mat_t Gm,fmpz_mod_poly_t G,int n,int k){
	fmpz_mat_zero(Gm);
	fmpz_t coeff;
	fmpz_init(coeff);


	for(int i=0;i<k;i++){
		for(int j=i;j<fmpz_mod_poly_length(G)+i;j++){
			fmpz_mod_poly_get_coeff_fmpz(coeff,G,j-i);
			fmpz_set(fmpz_mat_entry(Gm, i,j),coeff);		
		}
	}

	fmpz_clear(coeff);
}


void keygen(fmpz_mat_t key,int n,int k){
	//verif t
	int tt=(n-k);
	//verif gf
	
	 int m = log(n+1)/log(2);
	 printf("%d",m);
	 init_tabs(m);

	 //appel fonction d'init de var global de main.c
	 get_reduc(m);
	 lookuptab();

	 //fmpz_mat_t P,Gm,S,tmpM;
	 fmpz_mat_t tmpM;
	 fmpz_mod_poly_t G;
	 fmpz_t tmp;
	 fmpz_init_set_ui(tmp,pow(2,m));
	 fmpz_mod_poly_init(G,tmp);
	 fmpz_mat_init(P,n,n);
	 fmpz_mat_init(Gm,k,n);
	 fmpz_mat_init(S,k,k);
	 fmpz_mat_init(tmpM,k,n);
	 RPermut(P,n);
	 get_S(S,k);
	 gen_poly(G,tt);
	 //à supp
	 fmpz_t d;
	 fmpz_init_set_ui(d,gf);
	 fmpz_mod_poly_init(G1,d);
	 fmpz_mod_poly_set(G1,G);

	 //-----

	 printf("G\n");
	 fmpz_mod_poly_print(G);
	 get_G(Gm,G,n,k);
	 printf("P \n");
	 fmpz_mat_print_pretty(P);
	 printf("\n S \n");
	 fmpz_mat_print_pretty(S);
	 printf("\n Gm \n");
	 fmpz_mat_print_pretty(Gm);
	 printf("\nkey \n");
	 
	 
	 //fmpz_mat_t f;
	 //fmpz_mat_init(f,k,n);
	//multGFmat(tmpM,S,Gm,k,k,n);
	//multGFmat(f,tmpM,P,k,n,n);
	fmpz_mat_mul(tmpM,S,Gm);
	fmpz_mat_mul(Gm,tmpM,P);
	 

	//fmpz_mat_mul(tmpM,Gm,P);
	 //fmpz_mat_mul(Gm,S,tmpM);


	 fmpz_mat_set(key,Gm);
	 fmpz_mat_print_pretty(key);
	 printf("\n");



}



void Gen_E(fmpz_mat_t E,int n,int t){
	int temp;
	int j;
	int tabP[n]	;

	for(int i=0;i<t;i++){
		tabP[i]=1;
	}
	for(int i=t;i<n;i++){
		tabP[i]=0;
	}
	
	for(int i=n-1;i>0;i--){
		j = random()%i;
		temp=tabP[i];
		tabP[i]=tabP[j];
		tabP[j]=temp;

	}

	for(int i=0;i<n;i++){
		if(tabP[i]==1){
			fmpz_set_ui(fmpz_mat_entry(E, 0, i),random()%gf );
			
		}
		else{
			fmpz_set_ui(fmpz_mat_entry(E, 0, i), tabP[i]);
		}
	}
	printf("\n E \n");
	fmpz_mat_print_pretty(E);

}



void Encrypt_McEliece(fmpz_mat_t res,int* m,fmpz_mat_t key,int t,int k,int n){
	fmpz_mat_t msg;
	//fmpz_mat_t c;
	fmpz_mat_t E;
	//fmpz_t tmp;
	//fmpz_init(tmp);
	fmpz_mat_init(E,1,n);
	fmpz_mat_init(msg,1,k);
	for(int i=0;i<k;i++){
		fmpz_set_ui(fmpz_mat_entry(msg, 0,i),m[i]);	
	}
	printf("\n msg \n");
	fmpz_mat_print_pretty(msg);
	//fmpz_mat_mul(res,msg,key);
	multGFmat(res,msg,key,1,k,n);
	printf("\n C af mul \n");
	fmpz_mat_print_pretty(res);
	Gen_E(E,n,t);
	//E XOR res
	printf("\n E \n");
	fmpz_mat_print_pretty(E);
	for(int i=0;i<n;i++){
		fmpz_xor(fmpz_mat_entry(res, 0, i),fmpz_mat_entry(res, 0, i),fmpz_mat_entry(E, 0, i));
		//fmpz_set(fmpz_mat_entry(res, 0, i), tmp);

	}
	//fmpz_mat_add(res,res,E);
	printf("\n C \n");
	fmpz_mat_print_pretty(res);

}






// void Decrypt_McElieceMc_Eliece(fmpz_mat_t c,fmpz_mat_t P,fmpz_mat_t S,fmpz_mat_t G){
void Decrypt_McEliece(fmpz_mat_t c,int n,int k,int t){
	 
	fmpz_t d;
	fmpz_init_set_ui(d,gf);

	 fmpz_mat_t Si,Pi,res,fres;
	 fmpz_t den;
	 fmpz_init(den);
	 fmpz_mat_init(Pi,n,n);
	 fmpz_mat_init(Si,k,k);
	 fmpz_mat_init(res,1,n);
	 fmpz_mat_init(fres,1,k);
	 //printf("p inv ? : %d",fmpz_mat_inv(Pi , den , P));
	  fmpz_mat_transpose(Pi,P);
	  printf("\n pi \n");
	  fmpz_mat_print_pretty(Pi);
	 fmpz_mat_inv(Si , den , S);
	 

	 fmpz_mat_mul(res,c,Pi);
	//multGFmat(res,c,Pi,1,7,7);
	  printf("Si \n");
	 fmpz_mat_print_pretty(Si);
	 fmpz_mat_mul(res,c,Pi);
	 printf("C inter\n");
	 fmpz_mat_print_pretty(res);


	 fmpz_mod_poly_t data,received;
	 fmpz_mod_poly_init(data,d);
	 fmpz_mod_poly_init(received,d);
	 for(int i=0;i<n;i++){
	 	fmpz_mod_poly_set_coeff_fmpz(received,i,fmpz_mat_entry(res,0,i));
	 }
	 printf("received\n");
	 fmpz_mod_poly_print(received);
	 decode(data,received,2*t);
	 printf("data\n");
	 fmpz_mod_poly_print(data);
	 printf("\n");
	 fmpz_mat_clear(res);
	 fmpz_mat_init(res,1,k);

	// for(int i=0;i<k;i++){
	// 	fmpz_mod_poly_get_coeff_fmpz(fmpz_mat_entry(res,0,i),data,i);
	//  }
	//   printf("res \n");
	//  fmpz_mat_print_pretty(res);	
	// fmpz_mat_mul(fres,res,Si);
	// //multGFmat(fres,res,Si,1,k,k);
	//  printf("m \n");
	//  fmpz_mat_print_pretty(fres);


//---------
//  TEST

	 fmpz_mod_poly_t q;
	 fmpz_mod_poly_t r;

	 fmpz_mod_poly_init(q,d);
	 fmpz_mod_poly_init(r,d);
	 // fmpz_mod_poly_init(d1,d);
	 // fmpz_mod_poly_init(d2,d);


	 // fmpz_mod_poly_set_coeff_ui(d1,0,2);
	 // fmpz_mod_poly_set_coeff_ui(d1,1,3);
	 // fmpz_mod_poly_set_coeff_ui(d1,2,1);
	 // fmpz_mod_poly_set_coeff_ui(d1,3,2);
	 // fmpz_mod_poly_set_coeff_ui(d1,4,1);
	 // fmpz_mod_poly_set_coeff_ui(d1,5,0);
	 // fmpz_mod_poly_set_coeff_ui(d1,6,3);

	 // fmpz_mod_poly_set_coeff_ui(d2,0,3);
	 // fmpz_mod_poly_set_coeff_ui(d2,1,2);
	 // fmpz_mod_poly_set_coeff_ui(d2,2,1);
	 // fmpz_mod_poly_set_coeff_ui(d2,3,3);
	 // fmpz_mod_poly_set_coeff_ui(d2,4,1);

	 division(q,r,data,G1);
	 fmpz_mod_poly_print(q);
 	for(int i=0;i<k;i++){
		fmpz_mod_poly_get_coeff_fmpz(fmpz_mat_entry(res,0,i),q,i);
	  }
	  fmpz_mat_mul(res,res,Si);
	  printf("\nmessage décoder : \n");
	  fmpz_mat_print_pretty(res);	 




}




