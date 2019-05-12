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
#include "../inc/matrice.h"
#include <stdlib.h>
	

fmpz_mod_poly_t reduc;
fmpz_mod_poly_t* tab;
fmpz_t* ptoi;
fmpz_t DEUX;    
				// s  === bit par symbole   en param ou a calculer -> = ln(n+1)/ln(2) === log2(n+1)  / 
int gf;			// gf === 2^s 											                  n = 2^s - 1			
				// k  === nb symbole dans le message  k=n-2t
				// 2t === nb de symbole de contrôle  = n-k
				// n(taille encoder) = k+2t	

/* résumé des formules 

n = 2^s -1    
s = log2(n+1)
gf(2^s)

k  = n-2t
2t = n-k
n  = k+2t

*/

//long ?					long?
int convertOctalToDecimal(int octalNumber)
{
    int decimalNumber = 0, i = 0;

    while(octalNumber != 0)
    {
        decimalNumber += (octalNumber%10) * pow(8,i);
        ++i;
        octalNumber/=10;
    }

    i = 1;

    return decimalNumber;
}



void init_tabs(int m){
	fmpz_init_set_ui(DEUX,2);
	double m2 =pow(2,m);
	gf=m2;	
	tab=malloc(sizeof(fmpz_mod_poly_t)*(long int)m2);
	for(int i=0;i<m2;i++){
		fmpz_mod_poly_init(tab[i],DEUX);
	}
	ptoi=malloc(sizeof(fmpz_t)*(long int)m2);
	if(ptoi==NULL || tab==NULL){printf("soucis dans l'init des tab");exit(0);}
	for(int i=0;i<(int)m2;i++){
		fmpz_init(ptoi[i]);
	
	}
	fmpz_set_ui(ptoi[1],0);

}

void clear_tabs(){
	for(int i=0;i<gf;i++){
		fmpz_mod_poly_clear(tab[i]);
		fmpz_clear(ptoi[i]);
	}
	fmpz_clear(DEUX);
	fmpz_mod_poly_clear(reduc);
	free(tab);
	free(ptoi);
}


void get_reduc(int m){
	 FILE* fd = fopen("poly_irr.txt", "r");
	 if(fd==NULL){printf("fichier poly_irr.txt introuvable");exit(0);}
	 int tmp;
	 int tmp2;
	 int test;
	 fmpz_t mf;
	 while (tmp!=m && test){
	 	test=fscanf(fd,"%d %d",&tmp,&tmp2);
	 }

	 fmpz_init_set_ui(mf,m);
	 fmpz_mod_poly_init(reduc,mf);
	 tmp2=convertOctalToDecimal(tmp2);
	 fmpz_set_ui(mf,tmp2);
	setBinPoly(reduc,mf);
	 fmpz_clear(mf);
	 fclose(fd);
	 
}



int main(int argc, char** argv){
	//clear_tabs();
	// get_reduc(atoi(argv[1]));
	// init_tabs(atoi(argv[1]));


	// int n=7;
	// int k=3;
	// int t=2;
	int n=15;
	int k=9;
	int t=3;
	fmpz_mat_t key;
	fmpz_mat_init(key,k,n);

	keygen(key,n,k);
	//affichage des tableaux
	printf("tab:\n");
	for(int i=0;i<8;i++){
		printf("\n%d----",i);
		fmpz_mod_poly_print(tab[i]);
		printf("  alpha 	i: ");
		fmpz_print(ptoi[i]); 
	}
	fmpz_mat_t c;
	fmpz_mat_init(c,1,n);
	//int m[3]={3,5,6};
	int m[9]={5,6,3,8,7,5,11,10,15};
	Encrypt_McEliece(c,m,key,t,k,n);
	Decrypt_McEliece(c,n,k,t);


	return 0;
}