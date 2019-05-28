#include "../../inc/ReedSolomon/encryptRS.h"
#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <mpf2mpfr.h>
#include <flint/fmpz_mod_poly.h>



/**
	Calcule des éléments du corps Galois GF(2^m)
*/

void lookuptab(fmpz_mod_poly_t reduc, fmpz_mod_poly_t* tab,int primitive, int length,fmpz_t* ptoi  ){
	fmpz_t d,ev;
	fmpz_t tmp;
	fmpz_init(tmp);
	fmpz_init(ev);
	fmpz_init_set_ui(d,2);
	fmpz_mod_poly_set_coeff_ui(tab[0],0,1);
	fmpz_mod_poly_set_coeff_ui(tab[1],1,1);
	
	//tab[] comporte les valeurs de a^i
	for(int i=2;i<length;i++){
		fmpz_mod_poly_mulmod(tab[i] , tab[i-1] , tab[1],reduc);

	}
	fmpz_set_ui(ptoi[1],0);

	//ptoi[] comporte les exposants  i de a^i calculer
	for(int i=1;i<length;i++){
		tab[i]->p=16;
		fmpz_mod_poly_evaluate_fmpz(ev,tab[i],d);
		fmpz_set_ui(tmp,i);
		fmpz_set(ptoi[fmpz_get_ui(ev)],tmp);

	}

	fmpz_clear(d);
	fmpz_clear(tmp);
	fmpz_clear(ev);


}


/**
	 Cette fonction permet de calculer B

*/

void findB(fmpz_mod_poly_t res,fmpz_mod_poly_t A,fmpz_mod_poly_t G,fmpz_t* ptoi,fmpz_mod_poly_t* tab){
	

	fmpz_t deux;
	fmpz_init_set_ui(deux,2);	



	fmpz_t t;
	fmpz_init(t);
	fmpz_t t2;
	fmpz_init(t2);
	fmpz_t d;
	fmpz_init_set_ui(d,16);
	fmpz_t q;
	fmpz_init_set_ui(q,15);
	fmpz_mod_poly_t tmp;
	fmpz_mod_poly_init(tmp,d);
	fmpz_mod_poly_t Gtemp;
	fmpz_mod_poly_init(Gtemp,d);
	fmpz_mod_poly_t tmp2;
	fmpz_mod_poly_init(tmp2,d);
	fmpz_mod_poly_set(res,A);



	fmpz_mod_poly_t tmp3;
	fmpz_mod_poly_init(tmp3,deux);
	fmpz_mod_poly_t tmp4;
	fmpz_mod_poly_init(tmp4,deux);


	// set coeff en ^alpha  -> G et res
	for(int i=0;i<fmpz_mod_poly_degree(res)+1;i++){
		fmpz_mod_poly_get_coeff_fmpz(t,res,i);
		fmpz_mod_poly_set_coeff_fmpz(res,i,ptoi[fmpz_get_ui(t)]);
	}
	for(int i=0;i<fmpz_mod_poly_degree(G)+1;i++){
		fmpz_mod_poly_get_coeff_fmpz(t,G,i);
		fmpz_mod_poly_set_coeff_fmpz(Gtemp,i,ptoi[fmpz_get_ui(t)]);
	}


	//boucle
	while(fmpz_mod_poly_degree(res)>fmpz_mod_poly_degree(G)-1){
		
		//decalage puissance de X
		fmpz_mod_poly_set_coeff_ui(tmp,fmpz_mod_poly_degree(res)-fmpz_mod_poly_degree(G),1);
		fmpz_mod_poly_mul(tmp2,tmp,Gtemp);


		//mult des coeff
		fmpz_mod_poly_get_coeff_fmpz(t,res,fmpz_mod_poly_degree(res));
		
		if(fmpz_get_ui(t)!=15){
			for(int i=0;i<fmpz_mod_poly_degree(tmp2)+1;i++){
				fmpz_mod_poly_get_coeff_fmpz(t2,tmp2,i);
				if(!fmpz_is_zero(t2)){
					fmpz_add(t2,t2,t);
					fmpz_mod(t2,t2,q);
					if(fmpz_get_ui(t2)==0)fmpz_set_ui(t2,15);
						fmpz_mod_poly_set_coeff_fmpz(tmp2,i,t2);	
				}
			}
		}

		//xor coeff
		for(int i=0;i<fmpz_mod_poly_degree(res)+1;i++){
			fmpz_mod_poly_get_coeff_fmpz(t,res,i);
			fmpz_mod_poly_get_coeff_fmpz(t2,tmp2,i);

			if(!fmpz_is_zero(t))fmpz_mod_poly_set(tmp3,tab[fmpz_get_ui(t)]);else fmpz_mod_poly_zero(tmp3);
			if(!fmpz_is_zero(t2))fmpz_mod_poly_set(tmp4,tab[fmpz_get_ui(t2)]);else fmpz_mod_poly_zero(tmp4);
	

			fmpz_mod_poly_add(tmp3,tmp3,tmp4);
			tmp3->p=16;
			tmp4->p=16;
		

			fmpz_mod_poly_evaluate_fmpz(t,tmp3,deux);
			if(!fmpz_is_zero(t))fmpz_set(t,ptoi[fmpz_get_ui(t)]);

			fmpz_mod_poly_set_coeff_fmpz(res,i,t);
			tmp3->p=2;
			tmp4->p=2;
		}
		
		_fmpz_mod_poly_normalise(res);
		fmpz_mod_poly_zero(tmp);
	}


	// i -> alpha^i
	for(int i=0;i<fmpz_mod_poly_degree(res)+1;i++){
		fmpz_mod_poly_get_coeff_fmpz(t,res,i);
		fmpz_mod_poly_evaluate_fmpz(t,tab[fmpz_get_ui(t)],deux);
		fmpz_mod_poly_set_coeff_fmpz(res,i,t);
	}

	fmpz_clear(d);
	fmpz_clear(deux);
	fmpz_clear(q);
	fmpz_clear(t);
	fmpz_clear(t2);

	fmpz_mod_poly_clear(Gtemp);
	fmpz_mod_poly_clear(tmp);
	fmpz_mod_poly_clear(tmp2);
	fmpz_mod_poly_clear(tmp3);
	fmpz_mod_poly_clear(tmp4);



}


/**
	Cette fonction permet de générer le polynome Générateur
*/
void gen_poly(fmpz_mod_poly_t G,fmpz_mod_poly_t* tab,fmpz_t* ptoi,int tt,int tcycle){
	fmpz_t tmp;
	fmpz_t tmp2;
	fmpz_t tmp3;
	fmpz_t ev;
	fmpz_t d;
	fmpz_t s;
	fmpz_init(tmp);
	fmpz_init(ev);
	fmpz_init(tmp2);
	fmpz_init(tmp3);
	fmpz_init_set_ui(d,2);
	fmpz_init_set_ui(s,tcycle);
	fmpz_mod_poly_set_coeff_ui(G,0,2);	
	fmpz_mod_poly_set_coeff_ui(G,1,1);
	for(int i=2;i<tt+1;i++){
		fmpz_mod_poly_set_coeff_ui(G,i,1);
		for(int j=i-1;j>0;j--){
			fmpz_mod_poly_get_coeff_fmpz(tmp,G,j);
			if(!fmpz_is_zero(tmp)){
					fmpz_mod_poly_get_coeff_fmpz(tmp,G,j-1); //tmp -> g[j-1]
					fmpz_mod_poly_get_coeff_fmpz(tmp2,G,j); // tmp2 -> g[j]
					fmpz_set(tmp2,ptoi[fmpz_get_ui(tmp2)]); //tmp2 index[g[j]]
					fmpz_add_ui(tmp2,tmp2,i); // tmp2 +i  
					fmpz_mod(tmp2,tmp2,s); 	// tmp2 mod 15		
					fmpz_mod_poly_evaluate_fmpz(ev,tab[fmpz_get_ui(tmp2)],d);
					
					fmpz_xor(tmp,tmp,ev);
					
					fmpz_mod_poly_set_coeff_fmpz(G,j,tmp);

			}else{
				fmpz_mod_poly_get_coeff_fmpz(tmp,G,j-1);
				fmpz_mod_poly_set_coeff_fmpz(G,j,tmp);				
			}

		}
					
		//       gg[0] = alpha_to[(index_of[gg[0]]+i)%nn] ; 
					fmpz_mod_poly_get_coeff_fmpz(tmp,G,0); //tmp -> g[j-1]
					 // tmp2 -> g[j]
					fmpz_set(tmp,ptoi[fmpz_get_ui(tmp)]);
					fmpz_add_ui(tmp,tmp,i); // tmp2 +i  
					fmpz_mod(tmp,tmp,s);
					fmpz_mod_poly_evaluate_fmpz(ev,tab[fmpz_get_ui(tmp)],d);
					fmpz_mod_poly_set_coeff_fmpz(G,0,ev);					



	}

	fmpz_clear(tmp);
	fmpz_clear(tmp2);
	fmpz_clear(tmp3);
	fmpz_clear(ev);
	fmpz_clear(d);
	fmpz_clear(s);
}



/**
	
*/
void setBinPoly(fmpz_mod_poly_t res,fmpz_t f){
	fmpz_mod_poly_zero(res);
	for(size_t bi=fmpz_sizeinbase(f,2);bi>0;bi--){
		fmpz_mod_poly_set_coeff_ui(res,bi,fmpz_tstbit(f, bi));
	}
	fmpz_mod_poly_set_coeff_ui(res,0,fmpz_tstbit(f, 0));
}


/**

	Cette fonctio permet de calculer le produit de 2 polynome dans le champs galois
*/

void mulPoly(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2){
	fmpz_t c1,c2 ,d,tmp,tmp2;
	fmpz_init_set_ui(d,2);
	fmpz_init(c1);
	fmpz_init(c2);
	fmpz_init(tmp);
	fmpz_init(tmp2);
	fmpz_mod_poly_t ptmp;
	fmpz_mod_poly_t ptmp2;

	fmpz_mod_poly_init(ptmp,d);
	fmpz_mod_poly_init(ptmp2,d);
	
	for(long i=0;i<fmpz_mod_poly_degree(op1)+1;i++){
		for(long j=0;j<fmpz_mod_poly_degree(op2)+1;j++){ // pour chaqu ecouple (i,j)
			fmpz_mod_poly_get_coeff_fmpz(c1,op1,i);
			fmpz_mod_poly_get_coeff_fmpz(c2,op2,j);
			setBinPoly(ptmp,c1);
			setBinPoly(ptmp2,c2);
			fmpz_mod_poly_mulmod(ptmp,ptmp,ptmp2,reduc);
			ptmp->p=gf;
			ptmp2->p=gf;
			fmpz_mod_poly_evaluate_fmpz(tmp,ptmp,d);

			fmpz_mod_poly_get_coeff_fmpz(tmp2,res,i+j);
			fmpz_xor(tmp,tmp,tmp2);
			fmpz_mod_poly_set_coeff_fmpz(res,i+j,tmp);
			ptmp->p=2;
			ptmp2->p=2;


		}
	}
	fmpz_clear(d);
	fmpz_clear(tmp);
	fmpz_clear(tmp2);
	fmpz_clear(c1);
	fmpz_clear(c2);
	fmpz_mod_poly_clear(ptmp);
	fmpz_mod_poly_clear(ptmp2);



}


/**
	Fonction de l'encodage
*/ 
void encrypt(fmpz_mod_poly_t res,fmpz_mod_poly_t m, long t,long n){    // n utile? ->contenue dans taille du poly m
	fmpz_t tmp;
	fmpz_init_set_ui(tmp,gf);
	fmpz_mod_poly_t G;
	fmpz_mod_poly_init(G,tmp);
	//lookuptab();

	//Génération du polynome générateur
	gen_poly(G,2*t);
	fmpz_mod_poly_t B;
	fmpz_mod_poly_init(B,tmp);
	fmpz_mod_poly_t A;
	fmpz_mod_poly_init(A,tmp);
	fmpz_mod_poly_t R;
	fmpz_mod_poly_init(R,tmp);
	fmpz_mod_poly_t X;
	fmpz_mod_poly_init(X,tmp);
	fmpz_mod_poly_set_coeff_ui(X, 2*t, 1);
	
	//multiplication de  M*X^{2T}
	mulPoly(B,m,X);
	// calcule de B = A*X^{2t} par G
	findB(R,B,G);
	// Calcule du C finale  C=A*X^{2t} +B
	fmpz_mod_poly_add(res,B,R);

	fmpz_mod_poly_clear(A);
	fmpz_mod_poly_clear(B);
	fmpz_mod_poly_clear(G);
	fmpz_mod_poly_clear(R);
	fmpz_mod_poly_clear(X);
	fmpz_clear(tmp);
}
