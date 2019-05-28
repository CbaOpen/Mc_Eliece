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
long gf;	




														//row op1, col1=row2,col op2
void multGFmat(fmpz_mat_t res,fmpz_mat_t op1,fmpz_mat_t op2,long n,long p,long m){
	fmpz_t tmp,tmp2;
	fmpz_init(tmp);
	fmpz_init(tmp2);
	fmpz_t o1,o2;
	fmpz_init(o1);
	fmpz_init(o2);
	for(long i=0;i<n;i++){
		for(long j=0;j<m;j++){
			fmpz_init_set_ui(fmpz_mat_entry(res, i, j),0);
			for(long k=0;k<p;k++){
					fmpz_set(o1,fmpz_mat_entry(op1, i, k));

					//printf("k :   %ld\no1",k);
					//fmpz_print(o1);
					//printf("\n");	
					//fmpz_print(ptoi[fmpz_get_ui(o1)]);
					fmpz_set(o1,ptoi[fmpz_get_ui(o1)]);    //alplha op1
					fmpz_set(o2,fmpz_mat_entry(op2, k, j));
					fmpz_set(o2,ptoi[fmpz_get_ui(o2)]);
					if(fmpz_is_zero(o1) || fmpz_is_zero(o2)){
						fmpz_set_ui(o1,0);

					}
					else{
						fmpz_add(o1,o1,o2);
						fmpz_mod_ui(o1,o1,gf-1);
						tab[fmpz_get_ui(o1)]->p=gf;
						fmpz_set(o2,o1);
						fmpz_mod_poly_evaluate_fmpz(o1,tab[fmpz_get_ui(o1)],DEUX);
					}
					fmpz_xor(fmpz_mat_entry(res, i, j),fmpz_mat_entry(res, i, j),o1);
			}
		}
	}
	fmpz_clear(tmp);
	fmpz_clear(tmp2);
	fmpz_clear(o1);
	fmpz_clear(o2);
}
void RPermut(fmpz_mat_t P,long n){
long tabP[n]	;
long long j;
long temp;
srandom(getpid()+time(NULL));
	for(long i=0;i<n;i++){
		tabP[i]=i;
	}
	
	for(long i=n-1;i>0;i--){
		j = random()%i;
		temp=tabP[i];
		tabP[i]=tabP[j];
		tabP[j]=temp;

	}

	for(long i=0;i<n;i++){
		fmpz_set_ui(fmpz_mat_entry(P, i, tabP[i]), 1);
	}

}



void get_S(fmpz_mat_t S,long k){
	
	srandom(getpid()+time(NULL));
	flint_rand_t rand;
	flint_randinit(rand);
	fmpz_t det;
	fmpz_init(det);
	fmpz_mat_randrank(S , rand ,k , 1);
	long r=random()%100000;
	for(long i=0;i<r;i++){
		fmpz_mat_randrank(S , rand ,k , 1);
	}
	fmpz_clear(det);
	flint_randclear(rand);
}


void get_G(fmpz_mat_t Gm,fmpz_mod_poly_t G,long n,long k){
	fmpz_mat_zero(Gm);
	fmpz_t coeff;
	fmpz_init(coeff);


	for(long i=0;i<k;i++){
		for(long j=i;j<fmpz_mod_poly_length(G)+i;j++){
			fmpz_mod_poly_get_coeff_fmpz(coeff,G,j-i);
			fmpz_set(fmpz_mat_entry(Gm, i,j),coeff);		
		}
	}

	fmpz_clear(coeff);
}


void keygen(fmpz_mat_t key,long n,long k){
	//verif t
	long tt=(n-k);
	//verif gf
	
	long m = log(n+1)/log(2);
	init_tabs(m);

	 //appel fonction d'init de var global de main.c
	get_reduc(m);
	lookuptab();

	fmpz_mat_t tmpM,P,S,Gm;
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
	


	FILE* fpriv=fopen("KEYpriv","w");
	if(fpriv == NULL){printf("le fichier de la clé privé n'as pas pu être ouvert\n");exit(1);}
	get_G(Gm,G,n,k);
	fmpz_mat_fprint(fpriv,P);
	fputs(" ",fpriv);
	fmpz_mat_fprint(fpriv,Gm);
	fputs(" ",fpriv);
	fmpz_mat_fprint(fpriv,S);
	fmpz_mat_mul(tmpM,S,Gm);
	fmpz_mat_mul(Gm,tmpM,P);
	fmpz_mat_set(key,Gm);
	FILE* fpub=fopen("KEYpub","w");
	if(fpub == NULL){printf("le fichier de la clé public n'as pas pu être ouvert\n");exit(1);}
	fmpz_mat_fprint(fpub,key);
	
	fclose(fpub);
	fclose(fpriv);
	fmpz_mod_poly_clear(G);

	fmpz_clear(tmp);
	fmpz_mat_clear(P);
	fmpz_mat_clear(Gm);
	fmpz_mat_clear(S);
	fmpz_mat_clear(tmpM);

}



void Gen_E(fmpz_mat_t E,long n,long t){
	long temp;
	long j;
	long tabP[n]	;

	for(long i=0;i<t;i++){
		tabP[i]=1;
	}
	for(long i=t;i<n;i++){
		tabP[i]=0;
	}
	
	for(long i=n-1;i>0;i--){
		j = random()%i;
		temp=tabP[i];
		tabP[i]=tabP[j];
		tabP[j]=temp;

	}

	for(long i=0;i<n;i++){
		if(tabP[i]==1){
			fmpz_set_ui(fmpz_mat_entry(E, 0, i),(random()%(gf-1))+1 );
			
		}
		else{
			fmpz_set_ui(fmpz_mat_entry(E, 0, i), tabP[i]);
		}
	}

}



void Encrypt_McEliece(FILE* c,int* m,fmpz_mat_t key,long t,long k,long n){
	fmpz_mat_t msg,res;
	fmpz_mat_t E;
	fmpz_mat_init(E,1,n);
	fmpz_mat_init(msg,1,k);
	fmpz_mat_init(res,1,n);

	for(long i=0;i<k;i++){
		fmpz_set_ui(fmpz_mat_entry(msg, 0,i),m[i]);	
	}
	multGFmat(res,msg,key,1,k,n);
	Gen_E(E,n,t);
	//E XOR res
	for(long i=0;i<n;i++){
		fmpz_xor(fmpz_mat_entry(res, 0, i),fmpz_mat_entry(res, 0, i),fmpz_mat_entry(E, 0, i));

	}
	for(long i=0;i<n;i++){
		
		fputc(fmpz_get_ui(fmpz_mat_entry(res, 0, i)),c);

	}

	fmpz_mat_clear(msg);
	fmpz_mat_clear(E);
	fmpz_mat_clear(res);
}



// void Decrypt_McElieceMc_Eliece(fmpz_mat_t c,fmpz_mat_t P,fmpz_mat_t S,fmpz_mat_t G){
void Decrypt_McEliece(FILE* chiffre,long n,long k,long t,FILE* keys,FILE* msg){
	 
	fmpz_t d;
	fmpz_init_set_ui(d,gf);	


	fmpz_mat_t P,S,Gm;
	fmpz_mat_init(P,n,n);
	fmpz_mat_init(S,k,k);
	fmpz_mat_init(Gm,k,n);

	fmpz_mat_fread(keys,P);
	fgetc(keys);
	fmpz_mat_fread(keys,Gm);
	fgetc(keys);
	fmpz_mat_fread(keys,S);

	fmpz_mat_t c;
	fmpz_mat_init(c,1,n);

	for(long i=0;i<n;i++){
		fmpz_set_ui(fmpz_mat_entry(c,0,i),fgetc(chiffre));
	}




	 fmpz_mat_t Si,Pi,res,fres;
	 fmpz_t den;
	 fmpz_init(den);
	 fmpz_mat_init(Pi,n,n);
	 fmpz_mat_init(Si,k,k);
	 fmpz_mat_init(res,1,n);
	 fmpz_mat_init(fres,1,k);
	  fmpz_mat_transpose(Pi,P);
	 fmpz_mat_inv(Si , den , S);
	 

	 fmpz_mat_mul(res,c,Pi);
	//fmpz_mat_mul(res,c,Pi);


	 fmpz_mod_poly_t data,received;
	 fmpz_mod_poly_init(data,d);
	 fmpz_mod_poly_init(received,d);
	 for(long i=0;i<n;i++){
	 	fmpz_mod_poly_set_coeff_fmpz(received,i,fmpz_mat_entry(res,0,i));
	 }
	 decode(data,received,2*t);
	 fmpz_mat_clear(res);
	 fmpz_mat_init(res,1,k);

	 fmpz_mod_poly_t q;
	 fmpz_mod_poly_t r;
	 fmpz_mod_poly_t G;
	 fmpz_mod_poly_init(q,d);
	 fmpz_mod_poly_init(r,d);
	 fmpz_mod_poly_init(G,d);

	for(long i=0;i<=2*t;i++){
	 	fmpz_mod_poly_set_coeff_fmpz(G,i,fmpz_mat_entry(Gm,0,i));
	 }
	 division(q,r,data,G);
 	for(long i=0;i<k;i++){
		fmpz_mod_poly_get_coeff_fmpz(fmpz_mat_entry(res,0,i),q,i);
	  }
	  fmpz_mat_mul(res,res,Si);
	  for(long i=0;i<k;i++){
	  	fmpz_abs(fmpz_mat_entry(res, 0, i),fmpz_mat_entry(res, 0, i));

	  	fputc(fmpz_get_ui(fmpz_mat_entry(res, 0, i)),msg);
	  }

	  fmpz_mat_clear(P);
	  fmpz_mat_clear(S);
	  fmpz_mat_clear(Gm);
	  fmpz_mat_clear(res);
	  fmpz_mat_clear(Pi);
	  fmpz_mat_clear(Si);
	  fmpz_mat_clear(c);
	  fmpz_mat_clear(fres);

	  fmpz_mod_poly_clear(q);
	  fmpz_mod_poly_clear(r);
	  fmpz_mod_poly_clear(G);
	  fmpz_mod_poly_clear(data);
	  fmpz_mod_poly_clear(received);
	  fmpz_clear(d);
	  fmpz_clear(den);
}




