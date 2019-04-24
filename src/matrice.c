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
#include <stdlib.h>

void RPermut(fmpz_mat_t P,int n){
int[n] tabP;
long int j;
int temp;
	for(int i=0;i<n;i++){
		tabP[i]=i;
	}
	
	for(int i=n-1;i>=0;i--){
		j = random()%i;
		temp=tabP[i];
		tabP[i]=tabP[j];
		tabP|j]=temp;

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
	while(!fmpz_mat_det(det , S)){
		fmpz_mat_randrank(S , rand ,k , 1);

	}
}


void get_G(fmpz_mat_t Gm,fmpz_mod_poly_t G,int n){
	fmpz_mat_zero(Gm);
	fmpz_t coeff;
	fmpz_init(coeff);


	for(int i=0;i<k;i++){
		for(int j=i;j<fmpz_mod_poly_length(G)+i;j++){
			fmpz_mod_poly_get_coeff_fmpz(coeff,G,j-i)
			fmpz_set(fmpz_mat_entry(Gm, i,j),coeff);		
		}
	}

	fmpzz_clear(coeff);
}


void keygen(fmpz_mat_t key,int n,int k){
	//verif t
	int tt=(n-k);
	//verif gf
	
	 int m = ln(n+1)/ln(2);
	 //appel fonction d'init de var global de main.c

	 fmpz_mat_t P,Gm,S,tmpM;
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
	 genpoly(G,tt);
	 get_G(Gm,G,n);
	 fmpz_mat_mul(tmpM,S,Gm);
	 fmpz_mat_mul(Gm,tmpM,P);
	 fmpz_mat_set(key,Gm);


}




