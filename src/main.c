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
#include <string.h>
	

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
	test=fscanf(fd,"%d %d",&tmp,&tmp2);
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


void Get_msg(FILE* M,int *msg,int k){
	int tmp;
	for(int i=0;i<k;i++){
		tmp=fgetc(M);
		if(tmp==EOF){
			printf("le message ne contient pas le bon nombre de caractère\n");
			exit(1);
		}
		msg[i]=tmp;

	}
}


int main(int argc, char** argv){
	//verif arg et FILE 
	if(argc<2){
	 // help
	 exit(1); 
	}
	//mceliece -key n k 
	if(strcmp(argv[1],"-key")==0){
		int n=atoi(argv[2]);
		int k=atoi(argv[3]);
		if( (n-k)%2 != 0){
			printf("paramètre n et k incompatible \n");
			exit(1);
		}
		fmpz_mat_t key;
		fmpz_mat_init(key,k,n);
		keygen(key,n,k);
		fmpz_mat_clear(key);
		clear_tabs();
		return 0;

	}
	//mceliece -c KEYpub n k msg chiffré 
	if(strcmp(argv[1],"-c")==0){
		FILE * keyf=fopen(argv[2],"r");
		fmpz_mat_t key;
		int n=atoi(argv[3]);
		int k =atoi(argv[4]);
		if( (n-k)%2 != 0){
			printf("paramètre n et k incompatible \n");
			exit(1);
		}
		int t=(n-k)/2;
		fmpz_mat_init(key,k,n);
		fmpz_mat_fread(keyf,key);
		fclose(keyf);
		FILE* M=fopen(argv[5],"r");
		int *msg=malloc(sizeof(int)*k);
		Get_msg(M,msg,k);
		
		FILE* c=fopen(argv[6],"w");
		int m=log(n+1)/log(2);
		get_reduc(m);
		init_tabs(m);
		get_reduc(m);
		lookuptab();
		Encrypt_McEliece(c,msg,key,t,k,n);
		fclose(c);	
		fclose(M);
		clear_tabs();
		free(msg);
		fmpz_mat_clear(key);
		return 0;
	}

	//mceliece -d KEYpriv n k chiffre msg
	if(strcmp(argv[1],"-d")==0){
		FILE* keys=fopen(argv[2],"r");
		FILE* c=fopen(argv[5],"r");
		FILE* msg=fopen(argv[6],"w");
		int n=atoi(argv[3]);
		int k=atoi(argv[4]);
		if( (n-k)%2 != 0){
			printf("paramètre n et k incompatible \n");
			exit(1);
		}
		int t=(n-k)/2;
		int m=log(n+1)/log(2);
		get_reduc(m);
		init_tabs(m);
		get_reduc(m);
		lookuptab();
		Decrypt_McEliece(c,n,k,t,keys,msg);
		fclose(keys);
		fclose(c);
		fclose(msg);
		clear_tabs();
		return 0; 
	}
	//help
	return 1;

}