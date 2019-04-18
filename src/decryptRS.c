#include "../inc/encryptRS.h"
#include "../inc/decryptRS.h"
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
extern fmpz_mod_poly_t reduc;
extern fmpz_mod_poly_t* tab;
extern fmpz_t* ptoi;
extern fmpz_t DEUX;


void findC(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2, fmpz_t n){
	fmpz_mod_poly_t tmp;
	fmpz_mod_poly_init(tmp,n);
	fmpz_mod_poly_set_coeff_ui(tmp,6,1);
	fmpz_mod_poly_mul (res , op1 , tmp);
	fmpz_mod_poly_add (res , res, op2);	

}




void evaluation_fonction(fmpz_t res, fmpz_mod_poly_t fonction,int valeur_x){
	slong length;
	fmpz_t temp, temp_coef,val_x,d;
	fmpz_init_set_ui(d,2);
	fmpz_init(temp);
	fmpz_init(temp_coef);
	fmpz_init(val_x);

	fmpz_set(val_x,ptoi[valeur_x]);
	fmpz_mod_poly_get_coeff_fmpz ( res, fonction, 0 );
	length = fmpz_mod_poly_length ( fonction);

	for (int i=1;i<length;i++){

		fmpz_set(temp,val_x);
		fmpz_mul_ui(temp, temp, i);
		fmpz_mod_poly_get_coeff_fmpz ( temp_coef, fonction, i );
		fmpz_set(temp_coef, ptoi[fmpz_get_ui(temp_coef)]);
		if(fmpz_get_ui(temp_coef)!=0){
			fmpz_add(temp, temp, temp_coef);
			fmpz_mod_ui(temp, temp,15); /// A modifier qaund les variables deviendront globales

			fmpz_mod_poly_evaluate_fmpz(temp,tab[fmpz_get_ui(temp)],d);

			fmpz_xor(res, res, temp);
			
		}
		
		
	}

	fmpz_clear(temp);
	fmpz_clear(temp_coef);
	fmpz_clear(val_x);
	fmpz_clear(d);


}

void xor_polynome(fmpz_mod_poly_t resu, fmpz_mod_poly_t poly1_v, fmpz_mod_poly_t poly2_v){
	int deg1, deg2;
	fmpz_t temp1, temp2,n;
	fmpz_mod_poly_t poly1,poly2, res;
	fmpz_init_set_ui(n,16);
	fmpz_mod_poly_init(poly1, n);	
	fmpz_mod_poly_set(poly1,poly1_v);
	fmpz_mod_poly_init(poly2, n);	
	fmpz_mod_poly_set(poly2,poly2_v);
	fmpz_mod_poly_init(res, n);	

	fmpz_init(temp1);
	fmpz_init(temp2);

	deg1 = fmpz_mod_poly_length (poly1);
	deg2 = fmpz_mod_poly_length (poly2);

	if(deg1<=deg2){
		for (int i = 0; i < deg2; ++i)
		{ 
			fmpz_mod_poly_get_coeff_fmpz(temp1,poly1,i);
			fmpz_mod_poly_get_coeff_fmpz(temp2,poly2,i);
			fmpz_xor(temp1,temp1,temp2);
			fmpz_mod_poly_set_coeff_fmpz(res,i,temp1);
		
		}
	}
	else
		for (int i = 0; i < deg1; ++i){
			fmpz_mod_poly_get_coeff_fmpz(temp1,poly1,i);
			fmpz_mod_poly_get_coeff_fmpz(temp2,poly2,i);
			fmpz_xor(temp1,temp1,temp2);
			fmpz_mod_poly_set_coeff_fmpz(res,i,temp1);
			

		}
	fmpz_mod_poly_set(resu,res);



}

void mulPolyInDiv(fmpz_mod_poly_t res,fmpz_mod_poly_t pol2,fmpz_t valeur, int degre_poly )
{

	fmpz_mod_poly_t temp_res, poly_xor;
	fmpz_t n,d,temp;
	int deg;


	fmpz_init_set_ui(d,2);
	fmpz_init_set_ui(n,16);
	fmpz_init_set_ui(temp,16);
	fmpz_mod_poly_init(temp_res,n);
	fmpz_mod_poly_init(poly_xor,n);


	deg= fmpz_mod_poly_degree(pol2);
	
	for(int i=0; i<=deg;i++){
		fmpz_mod_poly_get_coeff_fmpz(temp,pol2,i);
		if(fmpz_get_ui(temp)!=0){
			fmpz_add(temp,ptoi[fmpz_get_ui(temp)],ptoi[fmpz_get_ui(valeur)]);
			fmpz_mod_ui(temp,temp,15);
			fmpz_mod_poly_evaluate_fmpz(temp,tab[fmpz_get_ui(temp)],d);}
		
		fmpz_mod_poly_set_coeff_fmpz(temp_res,i+degre_poly,temp);

	}


	fmpz_mod_poly_set(res,temp_res);

	fmpz_mod_poly_clear(temp_res);


}


void division(fmpz_mod_poly_t q, fmpz_mod_poly_t r, fmpz_mod_poly_t dividente, fmpz_mod_poly_t deviseur){
	int degr0, degr1,mult_deg;
	fmpz_mod_poly_t mult, poly_temp, poly_xor, r0,r1;
	fmpz_t alph,n,d,temp,temp1,temp2,mult_temp;

	fmpz_init_set_ui(d,2);

	fmpz_init(alph);
	fmpz_init(temp);
	fmpz_init(temp1);
	fmpz_init(temp2);
	fmpz_init(mult_temp);

	fmpz_init_set_ui(n,16);
	fmpz_mod_poly_init(mult, n);
	fmpz_mod_poly_init(poly_temp, n);
	fmpz_mod_poly_init(poly_xor, n);
	fmpz_mod_poly_init(r0, n);
	fmpz_mod_poly_init(r1, n);	
	fmpz_mod_poly_set(r0,dividente);
	fmpz_mod_poly_set(r1,deviseur);
	degr0 = fmpz_mod_poly_length (r0);
	degr1 = fmpz_mod_poly_length (r1);

	fmpz_mod_poly_init(q, n);
	fmpz_mod_poly_init(poly_temp, n);
	if (degr0<degr1){
		fmpz_mod_poly_set(r,r0);
		fmpz_mod_poly_zero(q);
	}

	else{

		do{
				
			fmpz_mod_poly_get_coeff_fmpz(temp1,r0,degr0-1);
			fmpz_mod_poly_get_coeff_fmpz(temp2,r1,degr1-1);
			fmpz_sub(temp,ptoi[fmpz_get_ui(temp1)],ptoi[fmpz_get_ui(temp2)]);
			fmpz_mod_ui(temp,temp,15);

			if(fmpz_get_ui(temp)==0){
				fmpz_set_ui(mult_temp,1);
				mult_deg=degr0-degr1;
				fmpz_mod_poly_set_coeff_ui(q,degr0-degr1,1);
			}
			else{
				fmpz_mod_poly_evaluate_fmpz(alph,tab[fmpz_get_ui(temp)],d);
				fmpz_set(mult_temp,alph);
				mult_deg=degr0-degr1;
				fmpz_mod_poly_set_coeff_fmpz(q,degr0-degr1,alph);

			}
			mulPolyInDiv(mult, r1,mult_temp,mult_deg);
			xor_polynome(poly_xor , mult , r0);
			degr0 = fmpz_mod_poly_length (poly_xor);			
			fmpz_mod_poly_set(r0,poly_xor);
			fmpz_mod_poly_init(poly_xor, n);
			fmpz_mod_poly_init(mult, n);


			}while(degr0>=degr1);
		fmpz_mod_poly_set(r,r0);
	}
	fmpz_clear(d);
	fmpz_clear(alph);
	fmpz_clear(temp);
	fmpz_clear(temp1);
	fmpz_clear(temp2);
	fmpz_clear(mult_temp);

	fmpz_clear(n);
	fmpz_mod_poly_clear(mult);
	fmpz_mod_poly_clear(poly_temp);
	fmpz_mod_poly_clear(poly_xor);
	fmpz_mod_poly_clear(r0);
	fmpz_mod_poly_clear(r1);	

}


void algo_euclide( fmpz_mod_poly_t localisation,fmpz_mod_poly_t amplitude , fmpz_mod_poly_t x2t_v, fmpz_mod_poly_t syndrome_v,int tt){
	fmpz_mod_poly_t poly_temp ,quotient, reste,multi,ampl,local,synd,x2t,syndrome;
	fmpz_t n;

	fmpz_init_set_ui(n,16);
	fmpz_mod_poly_init(poly_temp,n);
	fmpz_mod_poly_init(x2t,n);
	fmpz_mod_poly_init(syndrome,n);
	fmpz_mod_poly_init(ampl,n);
	fmpz_mod_poly_init(local,n);
	fmpz_mod_poly_init(synd,n);

	fmpz_mod_poly_init(quotient,n);
	fmpz_mod_poly_init(reste,n);
	fmpz_mod_poly_init(multi,n);
	fmpz_mod_poly_set(x2t,x2t_v);
	fmpz_mod_poly_set(syndrome,syndrome_v);

	fmpz_mod_poly_set_coeff_ui(local,0,0);
	fmpz_mod_poly_set_coeff_ui(ampl,0,1);
	fmpz_mod_poly_set(synd,syndrome);

	int i=2, deg_reste=0;


	while(deg_reste>(tt/2)|| i==2) {

		fmpz_mod_poly_init(reste,n);
		fmpz_mod_poly_init(quotient,n);
		fmpz_mod_poly_init(poly_temp,n);

		division(quotient,reste, x2t, synd);
		mulPoly(multi,ampl, quotient);
		xor_polynome(poly_temp ,multi, local );
		deg_reste = fmpz_mod_poly_length (reste);
		fmpz_mod_poly_init(multi,n);///////////////////////////////////////////////////
	
		fmpz_mod_poly_set(local,ampl);
		fmpz_mod_poly_set(ampl,poly_temp);
		
		fmpz_mod_poly_set(x2t,synd);
		fmpz_mod_poly_set(synd,reste);
		i++;
	}

	fmpz_mod_poly_set(localisation,poly_temp);
	fmpz_mod_poly_set(amplitude,reste);
}


bool calcul_poly_syndrome(fmpz_mod_poly_t syndrome, fmpz_mod_poly_t data, int tt){
	bool error=false;
	fmpz_t coef_sydrome,temp,d;
	fmpz_init_set_ui(d,2); 
	fmpz_init(coef_sydrome);
	fmpz_init(temp);

		for(int i=1; i<=tt;i++){
			 tab[i]->p=16;
			fmpz_mod_poly_evaluate_fmpz(temp,tab[i],d);
			evaluation_fonction(coef_sydrome, data,fmpz_get_ui(temp));
			fmpz_mod_poly_set_coeff_fmpz(syndrome, i-1,coef_sydrome);
			if(fmpz_get_ui(coef_sydrome)!=0)
				error=true;
		}

	fmpz_clear(temp);
	fmpz_clear(coef_sydrome);
	fmpz_clear(d);

	return error;

}

void derivation(fmpz_mod_poly_t res, fmpz_mod_poly_t function){
	fmpz_t temp;
	fmpz_init(temp);
	int deg=fmpz_mod_poly_length(function);

	for(int i=1; i<deg;i+=2){
		fmpz_mod_poly_get_coeff_fmpz(temp,function,i);
		fmpz_mod_poly_set_coeff_fmpz(res,i-1,temp);
	}
	fmpz_clear(temp);
}

int decode(fmpz_mod_poly_t data, fmpz_mod_poly_t received ,int nn, int tt ){
	bool error;
	int nbr_erreur=0;
	fmpz_mod_poly_t poly_temp, syndrome, localisation, amplitude,x2t, poly_position,poly_valeur, poly_racine,derive_localisation ,poly_error,test, t2, syndrome_bis, ampli_bis,ampli_mod,xt;
	fmpz_t n,temp,temp2,temp1,valeur,d;

	fmpz_init_set_ui(d,2); 
	fmpz_init_set_ui(n,nn);
	fmpz_init(temp);
	fmpz_init(temp1);
	fmpz_init(temp2);
	fmpz_init(valeur);
	fmpz_mod_poly_init(syndrome,n);
	fmpz_mod_poly_init(localisation,n);
	fmpz_mod_poly_init(amplitude,n);
	fmpz_mod_poly_init(x2t,n);
	fmpz_mod_poly_init(poly_position,n);
	fmpz_mod_poly_init(poly_valeur,n);
	fmpz_mod_poly_init(poly_racine,n);
	fmpz_mod_poly_init(derive_localisation,n);
	fmpz_mod_poly_init(poly_error,n);
	fmpz_mod_poly_init(poly_temp,n);
	fmpz_mod_poly_init(syndrome_bis,n);
	fmpz_mod_poly_init(ampli_bis,n);
	fmpz_mod_poly_init(ampli_mod,n);
	fmpz_mod_poly_init(xt,n);

	fmpz_mod_poly_init(test,n);

	fmpz_mod_poly_set_coeff_ui(test, 0, 1);
	fmpz_mod_poly_set_coeff_ui(test, 1, 9);
	fmpz_mod_poly_set_coeff_ui(test, 2,6);
		
	fmpz_mod_poly_init(t2,n);

	fmpz_mod_poly_set_coeff_ui(t2, 0, 11);
	fmpz_mod_poly_set_coeff_ui(t2, 1, 5);

	fmpz_mod_poly_set_coeff_ui(x2t, tt, 1);
	fmpz_mod_poly_set_coeff_ui(xt, tt/2, 1);


	error=calcul_poly_syndrome( syndrome, received,tt);

	fmpz_mod_poly_set(syndrome_bis,syndrome);
	if(error==false){
		fmpz_mod_poly_set(data,received);

	}
	else{ 	

		algo_euclide(localisation, amplitude ,x2t,syndrome,tt);
		mulPoly(ampli_bis,syndrome_bis,localisation);
		
		//Calcule position
		for(int i=1;i<=nn-1;i++){
			evaluation_fonction(valeur, localisation, i);
			if (fmpz_get_ui(valeur)==0)
			{
				fmpz_mod_poly_set_coeff_ui(poly_racine,nbr_erreur,i);
				fmpz_set(temp, ptoi[i]);
				fmpz_sub(temp,n,temp);
				fmpz_sub_ui(temp,temp,1);
				fmpz_mod_poly_set_coeff_fmpz(poly_position,nbr_erreur,temp);
				nbr_erreur++;
			}
			fmpz_init(valeur);
		}
		//Calcule amplitudes
		derivation(derive_localisation,localisation);

		for(int i=0;i<nbr_erreur;i++){
			fmpz_mod_poly_get_coeff_fmpz(temp,poly_racine,i);
			evaluation_fonction(temp1, derive_localisation, fmpz_get_ui(temp));
			evaluation_fonction(temp2, amplitude,  fmpz_get_ui(temp));
			fmpz_set(temp1,ptoi[fmpz_get_ui(temp1)]);
			fmpz_set(temp2,ptoi[fmpz_get_ui(temp2)]);
			fmpz_sub(temp1,temp2,temp1);
			fmpz_mod_ui(temp1,temp1,nn-1);
			fmpz_mod_poly_evaluate_fmpz(temp1,tab[fmpz_get_ui(temp1)],d);
			fmpz_mod_poly_get_coeff_fmpz(temp,poly_position,i);
			fmpz_mod_poly_set_coeff_fmpz(poly_error,fmpz_get_ui(temp),temp1);
			fmpz_init(temp);
			fmpz_init(temp1);
			fmpz_init(temp2);
			
		}
		xor_polynome(data, poly_error, received);
	}
	return 0;
}


void test_decode(){

	/********* Chiffrement ********/
	//init des tableaux
	

	fmpz_t deux;
	fmpz_init_set_ui(deux,2);

	// fmpz_mod_poly_t reduc;
	// fmpz_mod_poly_t* tab=malloc(sizeof(fmpz_mod_poly_t)*16);
	// for(int i=0;i<16;i++){
	// 	fmpz_mod_poly_init(tab[i],deux);
	// }
	// fmpz_t* ptoi=malloc(sizeof(fmpz_t)*16);
	// for(int i=0;i<16;i++){
	// fmpz_init(ptoi[i]);
	// }
	fmpz_t n;
	fmpz_t tmp;
	fmpz_init_set_ui(tmp,16);
	fmpz_init_set_ui(n,4);



	// fmpz_mod_poly_init(reduc, n);
	// fmpz_mod_poly_set_coeff_ui(reduc, 4, 1);
	// fmpz_mod_poly_set_coeff_ui(reduc, 1, 1);
	// fmpz_mod_poly_set_coeff_ui(reduc, 0, 1);
	
	fmpz_mod_poly_t G;
	fmpz_mod_poly_init(G,tmp);
	lookuptab();
	gen_poly(G,6,15);
	printf("g:\n");
	fmpz_mod_poly_print(G);
	printf("\n" );

/////----------------------------------------------

	


//affichage des tableaux
	printf("tab:\n");
	for(int i=0;i<16;i++){
		printf("\n%d ----",i);
		fmpz_mod_poly_print(tab[i]);
		printf("  alpha 	i: ");
		fmpz_print(ptoi[i]); 
	}
	

	//test B
	fmpz_mod_poly_t B;
	fmpz_mod_poly_init(B,tmp);
	fmpz_mod_poly_t A;
	fmpz_mod_poly_init(A,tmp);
	fmpz_mod_poly_t R;
	fmpz_mod_poly_init(R,tmp);
	
	fmpz_mod_poly_set_coeff_ui(A, 0, 2);
	fmpz_mod_poly_set_coeff_ui(A, 8, 8);
	fmpz_mod_poly_set_coeff_ui(A, 7, 7);
	fmpz_mod_poly_set_coeff_ui(A, 6, 6);
	fmpz_mod_poly_set_coeff_ui(A, 5, 5);
	fmpz_mod_poly_set_coeff_ui(A, 4, 4);
	fmpz_mod_poly_set_coeff_ui(A, 3, 3);
	fmpz_mod_poly_set_coeff_ui(A, 2, 2);
	fmpz_mod_poly_set_coeff_ui(A, 1, 1);
	

	fmpz_mod_poly_t X;
	fmpz_mod_poly_init(X,tmp);
	fmpz_mod_poly_set_coeff_ui(X, 6, 1);
	
	mulPoly(B,A,X);
	findB(R,B,G);

	fmpz_mod_poly_t received;
	fmpz_mod_poly_init(received,tmp);

	fmpz_mod_poly_add(received,B,R);

	printf("message clair = ");
	fmpz_mod_poly_print(A);
	printf("\n");
	printf("message chiffré = ");
	fmpz_mod_poly_print(received);
	printf("\n");
	
	/*********** Dechiffrement avec erreurs injectée ***********/
	fmpz_mod_poly_init(X,tmp);
	fmpz_mod_poly_set_coeff_ui(received, 12, 10);
	fmpz_mod_poly_set_coeff_ui(received, 11, 10);
	fmpz_mod_poly_set_coeff_ui(received, 10, 10);

	printf("\n");
	printf("message avec erreurs =");
	fmpz_mod_poly_print(received);
	printf("\n");
	decode(X,received,16,6);
	printf("Dechiffrement = ");
	fmpz_mod_poly_print(X);
	printf("\n");

	fmpz_mod_poly_clear(received);
	fmpz_mod_poly_clear(X);
	fmpz_mod_poly_clear(A);
	fmpz_mod_poly_clear(B);

}
