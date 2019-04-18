//#include <flint/flint.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <flint/fmpz_poly.h>
#include "encryptRS.h"
#include "decryptRS.h"

#include <math.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_poly.h>
#include "../inc/encryptRS.h"
#include "../inc/decryptRS.h"
#include <stdlib.h>


fmpz_mod_poly_t reduc;
fmpz_mod_poly_t* tab;
fmpz_t* ptoi;
fmpz_t DEUX;
int gf;



void init_tabs(int m){
	fmpz_init_set_ui(DEUX,2);
	double m2 =pow(2,m);
	gf=m2;	
	//printf("%lf\n",m2);
	tab=malloc(sizeof(fmpz_mod_poly_t)*(long int)m2);
	for(int i=0;i<m2;i++){
		fmpz_mod_poly_init(tab[i],DEUX);
	}
	ptoi=malloc(sizeof(fmpz_t)*(long int)m2);
	//if(ptoi==NULL)printf("soucis dans init ptoi");
	//printf("%d\n",(int)m2);
	for(int i=0;i<(int)m2;i++){
		fmpz_init(ptoi[i]);
	
	}
	fmpz_set_ui(ptoi[1],0);
}


void get_reduc(int m){
	 FILE* fd = fopen("poly_irr.txt", "r");
	 if(fd==NULL){printf("fichier introuvable");exit(0);}
	 int tmp;
	 int tmp2;
	 int test;
	 fmpz_t mf;
	 while (tmp!=m && test){
	 	test=fscanf(fd,"%d %d",&tmp,&tmp2);
	 }
	 
	 fmpz_init_set_ui(mf,m);
	 fmpz_mod_poly_init(reduc,mf);
	 //fmpz_mod_poly_set_coeff_ui(reduc,1,1);
	 //fmpz_mod_poly_set_coeff_ui(reduc,2,1);
	 //fmpz_mod_poly_set_coeff_ui(reduc,4,1);

	 fmpz_set_ui(mf,tmp2);

	 setBinPoly(reduc,mf);
	 printf("reduc\n");
	 fmpz_mod_poly_print(reduc);
	  printf("\n");
	 fmpz_clear(mf);
}



int main(int argc, char** argv){
	// test_encryptRS();
	
	get_reduc(atoi(argv[1]));
	init_tabs(atoi(argv[1]));
	test_decode();
	return 0;
}