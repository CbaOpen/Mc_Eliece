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


void lookuptab(fmpz_mod_poly_t reduc, fmpz_mod_poly_t* tab,int primitive, int length,fmpz_t* ptoi  ){
	fmpz_t d,ev;
	fmpz_t tmp;
	fmpz_init(tmp);
	fmpz_init(ev);
	fmpz_init_set_ui(d,2);
	fmpz_mod_poly_set_coeff_ui(tab[0],0,1);
	fmpz_mod_poly_set_coeff_ui(tab[1],1,1);
	
	for(int i=2;i<length;i++){
		fmpz_mod_poly_mulmod(tab[i] , tab[i-1] , tab[1],reduc);

	}
	fmpz_set_ui(ptoi[1],0);


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




void setBinPoly(fmpz_mod_poly_t res,fmpz_t f){
	fmpz_mod_poly_zero(res);
	for(size_t bi=fmpz_sizeinbase(f,2);bi>0;bi--){
		fmpz_mod_poly_set_coeff_ui(res,bi,fmpz_tstbit(f, bi));
	}
	fmpz_mod_poly_set_coeff_ui(res,0,fmpz_tstbit(f, 0));
}



void mulPoly(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2, fmpz_mod_poly_t reduc){
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
	
	for(int i=0;i<=fmpz_mod_poly_degree(op1)+1;i++){
		for(int j=0;j<=fmpz_mod_poly_degree(op2)+1;j++){
			fmpz_mod_poly_get_coeff_fmpz(c1,op1,i);
			fmpz_mod_poly_get_coeff_fmpz(c2,op2,j);
			setBinPoly(ptmp,c1);
			setBinPoly(ptmp2,c2);
			fmpz_mod_poly_mulmod(ptmp,ptmp,ptmp2,reduc);
			ptmp->p=16;//const
			ptmp2->p=16;//const
			fmpz_mod_poly_evaluate_fmpz(tmp,ptmp,d);

			fmpz_mod_poly_get_coeff_fmpz(tmp2,res,i+j);
			fmpz_xor(tmp,tmp,tmp2);
			fmpz_mod_poly_set_coeff_fmpz(res,i+j,tmp);
			ptmp->p=2;
			ptmp2->p=2;


		}
	}

	fmpz_clear(tmp);
	fmpz_clear(tmp2);
	fmpz_clear(c1);
	fmpz_clear(c2);
	fmpz_mod_poly_clear(ptmp);
	fmpz_mod_poly_clear(ptmp2);



}

void findC(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2, fmpz_t n){
	fmpz_mod_poly_t tmp;
	fmpz_mod_poly_init(tmp,n);
	fmpz_mod_poly_set_coeff_ui(tmp,6,1);
	fmpz_mod_poly_mul (res , op1 , tmp);
	fmpz_mod_poly_add (res , res, op2);	

}




void evaluation_fonction(fmpz_t res, fmpz_mod_poly_t fonction,  fmpz_mod_poly_t* tab,fmpz_t* ptoi, int valeur_x){
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

void xor_polymone(fmpz_mod_poly_t resu, fmpz_mod_poly_t poly1_v, fmpz_mod_poly_t poly2_v){
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
void multiplication_poly(fmpz_mod_poly_t res, fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2, fmpz_mod_poly_t* tab,fmpz_t* ptoi ){
	int deg1, deg2;
//printf("toto");
	fmpz_t temp, temp1, temp2,coef,d,t,n;
	fmpz_init(temp);
	fmpz_init(coef);
	fmpz_init(t);

	fmpz_init(temp1);
	fmpz_init(temp2);

	fmpz_init_set_ui(d,2);
	fmpz_init_set_ui(n,16);

	deg1 = fmpz_mod_poly_length (pol1);
	deg2 = fmpz_mod_poly_length (pol2);
	//fmpz_mod_poly_init(res,n);

	for(int i=0;i<deg1;i++){
		for(int j=0;j<deg2;j++){
			fmpz_mod_poly_get_coeff_fmpz(temp1,pol1,i);
			fmpz_mod_poly_get_coeff_fmpz(temp2,pol2,j);
			fmpz_mod_poly_get_coeff_fmpz(coef, res, i+j);
			fmpz_sub(coef,coef,coef);
			/*if(fmpz_get_ui(temp2)==0 &&fmpz_get_ui(temp1)!=0){
				fmpz_set_ui(temp,0);
				fmpz_mod_poly_set_coeff_fmpz(res,i+j,coef);

			}

			else if (fmpz_get_ui(temp2)!=0 &&fmpz_get_ui(temp1)==0){
				fmpz_set_ui(temp,0);
				fmpz_mod_poly_set_coeff_fmpz(res,i+j,coef);

			}
			else if (fmpz_get_ui(temp2)==0 &&fmpz_get_ui(temp1)==0){
				fmpz_set_ui(temp,0);
				fmpz_mod_poly_set_coeff_fmpz(res,i+j,coef);

			}

			else{
				fmpz_add(temp,ptoi[fmpz_get_ui(temp1)],ptoi[fmpz_get_ui(temp2)]);
				fmpz_mod_ui(temp,temp,15);

				if(fmpz_get_ui(temp)==0){
					fmpz_set_ui(t,1);
					fmpz_xor(temp,t,coef);
					fmpz_mod_poly_set_coeff_fmpz(res,i+j,temp);

				}

				else{
					fmpz_mod_poly_evaluate_fmpz(temp,tab[fmpz_get_ui(temp)],d);
					fmpz_xor(temp,temp,coef);
					fmpz_mod_poly_set_coeff_fmpz(res,i+j,temp);
				}

			}*/
			if(fmpz_get_ui(temp2)!=0 &&fmpz_get_ui(temp1)!=0){
				//fmpz_add(temp,ptoi[fmpz_get_ui(temp1)],ptoi[fmpz_get_ui(temp2)]); fmpz_print(temp); printf("\n");
				fmpz_mod_ui(temp,temp,15);

		
				
					fmpz_mod_poly_evaluate_fmpz(temp,tab[fmpz_get_ui(temp)],d);
					fmpz_xor(temp,temp,coef);
					fmpz_mod_poly_set_coeff_fmpz(res,i+j,temp);
			

			}
		}
	}
	//fmpz_mod_poly_print(res);

}


void muliplication(fmpz_mod_poly_t res,fmpz_mod_poly_t pol2, fmpz_mod_poly_t* tab,fmpz_t* ptoi,fmpz_t valeur, int degre_poly )
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

void division(fmpz_mod_poly_t q, fmpz_mod_poly_t r, fmpz_mod_poly_t dividente, fmpz_mod_poly_t deviseur,fmpz_mod_poly_t* tab, fmpz_t* ptoi,fmpz_mod_poly_t reduc){
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
				//fmpz_mod_poly_set_coeff_fmpz(poly_temp,degr0-degr1,alph);
				fmpz_set_ui(mult_temp,1);
				mult_deg=degr0-degr1;
				fmpz_mod_poly_set_coeff_ui(q,degr0-degr1,1);
			}
			else{
				fmpz_mod_poly_evaluate_fmpz(alph,tab[fmpz_get_ui(temp)],d);
				//fmpz_mod_poly_set_coeff_fmpz(poly_temp,degr0-degr1,alph);
				fmpz_set(mult_temp,alph);
				mult_deg=degr0-degr1;
				fmpz_mod_poly_set_coeff_fmpz(q,degr0-degr1,alph);

			}


			
			//mulPoly(mult,poly_temp, r1, reduc);
			//multiplication_poly(mult, poly_temp, r1,tab,ptoi ); fmpz_mod_poly_print(poly_temp); printf("\n");
			muliplication(mult, r1,tab,ptoi,mult_temp,mult_deg);

			xor_polymone(poly_xor , mult , r0);

			degr0 = fmpz_mod_poly_length (poly_xor);			

			fmpz_mod_poly_set(r0,poly_xor);
			fmpz_mod_poly_init(poly_xor, n);
						fmpz_mod_poly_init(mult, n);


//printf("TRRRRRRRRRRRR   %d",degr0);	
		}while(degr0>=degr1);
		fmpz_mod_poly_set(r,r0);
		printf("yiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii\n");


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


void algo_euclide( fmpz_mod_poly_t localisation,fmpz_mod_poly_t amplitude , fmpz_mod_poly_t x2t_v, fmpz_mod_poly_t syndrome_v,int tt,fmpz_mod_poly_t* tab, fmpz_t* ptoi,fmpz_mod_poly_t reduc){
	fmpz_mod_poly_t poly_temp ,quotient, reste,multi,ampl,local,synd,x2t,syndrome;
	fmpz_t n,temp1,temp2;

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

	int i=2, deg_reste=0, deg=0;


	while(deg_reste>(tt/2)|| i==2) {
		//if(deg_reste>(tt/2) || i==2){

		//	printf(" \n ----------------------------- ");

fmpz_mod_poly_init(reste,n);
fmpz_mod_poly_init(quotient,n);
fmpz_mod_poly_init(poly_temp,n);

		division(quotient,reste, x2t, synd,tab,ptoi,reduc);
		//multiplication_poly( mult,ampl, quotient,tab,ptoi );

printf("\n \n  \n reste \n");
		fmpz_mod_poly_print(reste);
printf("\n ampl \n");

		fmpz_mod_poly_print(ampl);
		printf(" \n quo \n");
		fmpz_mod_poly_print(quotient);
		
		mulPoly(multi,ampl, quotient, reduc);
		printf(" \n mult \n");
		fmpz_mod_poly_print(multi);
		
		printf("\n ");
		fmpz_mod_poly_print(local);
		xor_polymone(poly_temp ,multi, local ); 	printf("  xooor\n ");
		fmpz_mod_poly_print(poly_temp);
			printf("\n ");
		fmpz_mod_poly_print(poly_temp);


		deg_reste = fmpz_mod_poly_length (reste);
		fmpz_mod_poly_init(multi,n);//////////////////////////////////////////////////////
	
		fmpz_mod_poly_set(local,ampl);
		fmpz_mod_poly_set(ampl,poly_temp);
		
		fmpz_mod_poly_set(x2t,synd);
		fmpz_mod_poly_set(synd,reste);
		
		
		i++;
		
	

	}

	fmpz_mod_poly_set(localisation,poly_temp);
	fmpz_mod_poly_set(amplitude,reste);

}


bool calcul_poly_syndrome(fmpz_mod_poly_t syndrome, fmpz_mod_poly_t data,fmpz_mod_poly_t* tab,fmpz_t* ptoi, int tt){
	bool error=false;
	fmpz_t coef_sydrome,temp,d;
	fmpz_init_set_ui(d,2); 
	fmpz_init(coef_sydrome);
	fmpz_init(temp);

		for(int i=1; i<=tt;i++){
			 tab[i]->p=16;
			fmpz_mod_poly_evaluate_fmpz(temp,tab[i],d);
			evaluation_fonction(coef_sydrome, data, tab,ptoi, fmpz_get_ui(temp));
			fmpz_mod_poly_set_coeff_fmpz(syndrome, i-1,coef_sydrome);
			if(fmpz_get_ui(coef_sydrome)!=0)
				error=true;
		}

	fmpz_clear(temp);
	fmpz_clear(coef_sydrome);
	fmpz_clear(d);


printf("%d \n", error);
	return error;

}

void derivation(fmpz_mod_poly_t res, fmpz_mod_poly_t function,fmpz_mod_poly_t* tab,fmpz_t* ptoi){
	int deg;

	fmpz_t temp;

	fmpz_init(temp);

	deg=fmpz_mod_poly_length(function);
printf("\n");
	fmpz_mod_poly_print(function);
 printf("mkljghkjlmùjhgfdgh  %d \n", deg);
	for(int i=1; i<deg;i+=2){
		fmpz_mod_poly_get_coeff_fmpz(temp,function,i);


		
			fmpz_mod_poly_set_coeff_fmpz(res,i-1,temp);

			fmpz_mod_poly_print(res);
			printf("\n");

		/*if(fmpz_get_ui(temp1)!=0&&fmpz_get_ui(temp2)!=0){

			fmpz_add(temp1, temp2,temp1);
			fmpz_mod_ui(temp1,temp1,15);
			fmpz_mod_poly_evaluate_fmpz(temp2,tab[fmpz_get_ui(temp1)],d);

		}*/
		//else
		//	fmpz_mod_poly_set_coeff_ui(res,i-1,0);
		
	}

		fmpz_clear(temp);



}





int decode(fmpz_mod_poly_t data, fmpz_mod_poly_t received ,fmpz_mod_poly_t* tab,fmpz_t* ptoi, int nn, int tt,fmpz_mod_poly_t reduc ){
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


	error=calcul_poly_syndrome( syndrome, received,tab, ptoi, tt);
printf("\n syndrome \n"); fmpz_mod_poly_print(syndrome);

//fmpz_mod_poly_set(syndrome_bis,syndrome);


/*fmpz_mod_poly_t AI;
fmpz_mod_poly_init(AI,n);
	fmpz_mod_poly_set_coeff_ui(AI, 5, 1);
 	fmpz_mod_poly_set_coeff_ui(AI, 4, 15);
 	fmpz_mod_poly_set_coeff_ui(AI, 3, 7);
 	fmpz_mod_poly_set_coeff_ui(AI, 2, 8);
 	fmpz_mod_poly_set_coeff_ui(AI, 0, 11);*/
fmpz_mod_poly_set(syndrome_bis,syndrome);
	if(error==false){
		fmpz_mod_poly_set(data,received);
				//printf("\n ERROr");

	}
	else{ 	

		algo_euclide(localisation, amplitude ,x2t,syndrome,tt,tab, ptoi,reduc);
		
/*
		fmpz_mod_poly_init(localisation,n);
		fmpz_mod_poly_init(amplitude,n);
fmpz_mod_poly_set_coeff_ui(localisation, 2, 6);
fmpz_mod_poly_set_coeff_ui(localisation, 1,9);
fmpz_mod_poly_set_coeff_ui(localisation, 0, 1);

fmpz_mod_poly_set_coeff_ui(amplitude, 1, 5);
fmpz_mod_poly_set_coeff_ui(amplitude,0 ,11);*/

printf("\n localisation \n"); fmpz_mod_poly_print(localisation);
		printf("\n amplitude \n"); fmpz_mod_poly_print(amplitude);





		mulPoly(ampli_bis,syndrome_bis,localisation,reduc);

		printf("\n ampli bis\n");
		fmpz_mod_poly_print(ampli_bis);

		findB(ampli_mod,ampli_bis,xt,ptoi,tab);
		printf("\n ampli mo\n");
		fmpz_mod_poly_print(ampli_mod);

		/*if (fmpz_mod_poly_degree(amplitude)==0)
			return 1;*/
		////////////////Calcule position
		for(int i=1;i<=nn-1;i++){
			evaluation_fonction(valeur, localisation, tab, ptoi, i);
		//	evaluation_fonction(valeur, test, tab, ptoi, i);
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
		printf("\n RACINE\n");
		fmpz_mod_poly_print(poly_racine	);
		///Calcule amplitudes
		derivation(derive_localisation,localisation,tab,ptoi);
		//derivation(derive_localisation,test,tab,ptoi);

		for(int i=0;i<nbr_erreur;i++){
			fmpz_mod_poly_get_coeff_fmpz(temp,poly_racine,i);
			evaluation_fonction(temp1, derive_localisation, tab,ptoi, fmpz_get_ui(temp));
			evaluation_fonction(temp2, amplitude, tab,ptoi, fmpz_get_ui(temp));
			//evaluation_fonction(temp2, t2, tab,ptoi, fmpz_get_ui(temp));


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
		printf("\n \n position \n");
		fmpz_mod_poly_print(poly_position);
printf("\n received \n ");
		fmpz_mod_poly_print(received);

		printf("\n ERROr\n");
		fmpz_mod_poly_print(poly_error);

		xor_polymone(data, poly_error, received);


	}




	return 0;
}

int main(){
////// voir à déclarer de manière global

	//init des tableaux
	fmpz_t deux;
	fmpz_init_set_ui(deux,2);
	fmpz_mod_poly_t reduc;
	fmpz_mod_poly_t* tab=malloc(sizeof(fmpz_mod_poly_t)*16);
	for(int i=0;i<16;i++){
		fmpz_mod_poly_init(tab[i],deux);
	}
	fmpz_t* ptoi=malloc(sizeof(fmpz_t)*16);
	for(int i=0;i<16;i++){
	fmpz_init(ptoi[i]);
	}

	fmpz_t n;
	fmpz_t tmp;
	fmpz_init_set_ui(tmp,16);
	fmpz_init_set_ui(n,4);



	fmpz_mod_poly_init(reduc, n);
	fmpz_mod_poly_set_coeff_ui(reduc, 4, 1);
	fmpz_mod_poly_set_coeff_ui(reduc, 1, 1);
	fmpz_mod_poly_set_coeff_ui(reduc, 0, 1);
	
	fmpz_mod_poly_t G;
	fmpz_mod_poly_init(G,tmp);
	lookuptab(reduc,tab,2, 16,ptoi);
	gen_poly(G,tab,ptoi,6,15);
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
	
	mulPoly(B,A,X,reduc);
	findB(R,B,G,ptoi,tab);



	fmpz_mod_poly_init(G,tmp);

	fmpz_mod_poly_add(G,B,R);
	
	printf("\nB\n \n\n\n");
	fmpz_mod_poly_print(G);
	
	fmpz_mod_poly_init(X,tmp);

calcul_poly_syndrome( X, G,tab, ptoi, 6);
	printf("\n");
fmpz_mod_poly_print(X);
printf("\n");
fmpz_mod_poly_set_coeff_ui(G, 10, 14);
//fmpz_mod_poly_set_coeff_ui(G, 12, 12);
//fmpz_mod_poly_set_coeff_ui(G, 11, 15);

fmpz_mod_poly_print(G);
fmpz_mod_poly_init(X,tmp);
decode(X,G,tab,ptoi, 16,6,reduc );
printf("\n");

fmpz_mod_poly_print(X);
printf("\n");


/*fmpz_mod_poly_init(G,tmp);
fmpz_mod_poly_set_coeff_ui(G, 2, 10);
fmpz_mod_poly_set_coeff_ui(G, 1, 5);
fmpz_mod_poly_set_coeff_ui(G, 0, 10); printf("\n \n");*/

//fmpz_mod_poly_print(G);
fmpz_mod_poly_init(X,tmp);
derivation(X, G,tab,ptoi);
printf("\n G\n");
fmpz_mod_poly_print(X);

fmpz_mod_poly_init(X,tmp);
fmpz_mod_poly_t Q1,Q,K;
	fmpz_mod_poly_init(Q1,tmp);

	fmpz_mod_poly_init(K,tmp);
	fmpz_mod_poly_init(Q,tmp);

 	fmpz_mod_poly_set_coeff_ui(Q1, 5, 5); 
 	fmpz_mod_poly_set_coeff_ui(Q1, 4, 2); 
 	fmpz_mod_poly_set_coeff_ui(Q1, 2, 15); 
 	fmpz_mod_poly_set_coeff_ui(Q1, 1, 11);
 	fmpz_mod_poly_set_coeff_ui(Q1, 0, 12);//12


 	fmpz_mod_poly_set_coeff_ui(Q, 0, 12);

//fmpz_mod_poly_set_coeff_ui(Q, 0, 1);
//multiplication_poly(K, Q, Q1, tab,ptoi );
/*mulPoly(K,Q, Q1, reduc);	
printf("\n multypll \n");
fmpz_mod_poly_print(K);

fmpz_t valeur;
fmpz_init_set_ui(valeur,12);
muliplication(X, Q1,tab,ptoi,valeur,0);
//mulPoly(X,Q, Q1, reduc);
printf("\n multypll \n");
fmpz_mod_poly_print(X);*/



fmpz_mod_poly_t test;
fmpz_mod_poly_init(test,tmp);
fmpz_mod_poly_set_coeff_ui(test,2,9);
fmpz_mod_poly_set_coeff_ui(test,1,5);
//fmpz_mod_poly_set_coeff_ui(test,0,1);
fmpz_t truc,truc2;
fmpz_init(truc);
evaluation_fonction(truc,test,tab,ptoi,6);
printf("\n test eval\n");
fmpz_print(truc);


// //test mult poly
// 	fmpz_mod_poly_t res,op1,op2;
// 	fmpz_mod_poly_init(res,tmp);
// 	fmpz_mod_poly_init(op1,tmp);
// 	fmpz_mod_poly_init(op2,tmp);

// 	//fmpz_mod_poly_set_coeff_ui(op1, 1, 1);
// 	fmpz_mod_poly_set_coeff_ui(op1, 0, 5);
// 	//fmpz_mod_poly_set_coeff_ui(op2, 1, 1);
// 	fmpz_mod_poly_set_coeff_ui(op2, 0, 13);
// 	mulPoly(res,op1,op2,reduc);
// 	fmpz_mod_poly_print(res);
// 	printf("\n");
// fmpz_mod_poly_t T;
// 	fmpz_mod_poly_init(T,tmp);
	
// 	fmpz_mod_poly_set_coeff_ui(T, 8, 9);
// 	fmpz_mod_poly_set_coeff_ui(T, 7, 8);
// 	fmpz_mod_poly_set_coeff_ui(T, 6, 7);
// 	fmpz_mod_poly_set_coeff_ui(T, 5, 6);
// 	fmpz_mod_poly_set_coeff_ui(T, 4, 5);
// 	fmpz_mod_poly_set_coeff_ui(T, 3, 4);
// 	fmpz_mod_poly_set_coeff_ui(T, 2, 3);
// 	fmpz_mod_poly_set_coeff_ui(T, 1, 2);
// 	fmpz_mod_poly_set_coeff_ui(T, 0, 1);

// fmpz_mod_poly_t AI;
// 	fmpz_mod_poly_init(AI,tmp);
	
	
// 	fmpz_mod_poly_set_coeff_ui(AI, 5, 1);
// 	fmpz_mod_poly_set_coeff_ui(AI, 4, 15);
// 	fmpz_mod_poly_set_coeff_ui(AI, 3, 7);
// 	fmpz_mod_poly_set_coeff_ui(AI, 2, 8);
// 	fmpz_mod_poly_set_coeff_ui(AI, 0, 11);





// fmpz_mod_poly_t BI;
// 	fmpz_mod_poly_init(BI,tmp);
// fmpz_mod_poly_t C;
// 	fmpz_mod_poly_init(C,tmp);
// fmpz_mod_poly_t D;
// 	fmpz_mod_poly_init(D,tmp);
// fmpz_mod_poly_t E;
// 	fmpz_mod_poly_init(E,tmp);
	
// 	fmpz_mod_poly_set_coeff_ui(BI, 6, 1);
	

// fmpz_mod_poly_gcd ( C , AI , BI );

// printf("C -------------------------\n");
// 		fmpz_mod_poly_print( C );printf("\n");


// printf("A\n");
// 	/*fmpz_mod_poly_print(T);
// 	findC(T,T,R,tmp);


// 	printf("JKJK------------------\n");
// 	fmpz_mod_poly_print(T);*/
// //1x5 + 15x4 + 7x3 + 8x2 + 11
	

// fmpz_mod_poly_factor_t try;
// 	fmpz_mod_poly_factor_init(try);
// 	//	fmpz_mod_poly_factor_print( try );

// 		printf("\n");
// fmpz_mod_poly_xgcd_euclidean( C,D,E ,AI ,BI);
// printf("AI\n");
// 		fmpz_mod_poly_print( AI);printf("\n");
// printf("BI\n");
// 		fmpz_mod_poly_print( BI);printf("\n");
// printf("C\n");
// 		fmpz_mod_poly_print( C );printf("\n");
// 		printf("D \n");
// 		fmpz_mod_poly_print( D );printf("\n");
// 		printf("E \n");
// 		fmpz_mod_poly_print( E );
// printf(" \n");
// /*fmpz_mod_poly_t tr;
// 	fmpz_mod_poly_init(tr,tmp-1);

// 	fmpz_mod_poly_set_coeff_ui(tr, 6, 1);
// 	fmpz_mod_poly_set_coeff_ui(tr, 5, 3);
// 	fmpz_mod_poly_set_coeff_ui(tr, 4, 1);
// 	fmpz_mod_poly_set_coeff_ui(tr, 3, 4);
// 	fmpz_mod_poly_set_coeff_ui(tr, 2, 7); 
// 	fmpz_mod_poly_set_coeff_ui(T, 1, 13);
// 	fmpz_mod_poly_set_coeff_ui(T, 0, 15);


// 		fmpz_mod_poly_factor_t factors;
// 		fmpz_mod_poly_factor_init(factors);
// 		fmpz_mod_poly_factor_berlekamp(factors , D);
// 		printf("facccccccccccccccccccccctor\n");
// 		fmpz_mod_poly_factor_print(factors);

// 		int test=5;
// 		test=nmod_poly_is_irreducible (E);

// 		printf("\n  %d   \n",test);



// printf("\nAIIIIIIIIIIIIIIIIII\n");
// 		fmpz_mod_poly_print(AI);
// fmpz_mod_poly_derivative ( tr, AI);

// printf("\nDERRRRREEEEEEEEEEEV\n");
// 		fmpz_mod_poly_print(tr);
// slong test;
// printf("\n \n");
// fmpz_mod_poly_print(AI);
// 	test= fmpz_mod_poly_length (  AI );
// printf("\n %d \n", test);*/

// fmpz_mod_poly_t resu,poly_sy;
// fmpz_t resutltat,temp,d;
// fmpz_init(resutltat);
// fmpz_init_set_ui(d,2);
// fmpz_init(temp);

// 	fmpz_mod_poly_init(resu,tmp-1);
// fmpz_mod_poly_init(poly_sy,tmp-1);
// 	fmpz_mod_poly_set_coeff_ui(resu, 14, 9); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 13, 8); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 12, 7); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 11, 1); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 10, 5); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 9, 4); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 8, 3); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 7, 2); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 6, 1); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 5, 6); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 4, 15); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 3, 15); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 2, 5); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 1, 11);
// 	fmpz_mod_poly_set_coeff_ui(resu, 0, 14);

// 	/*fmpz_mod_poly_set_coeff_ui(resu, 2, 6); 
// 	fmpz_mod_poly_set_coeff_ui(resu, 1, 9);
// 	fmpz_mod_poly_set_coeff_ui(resu, 0, 1);*/
// 			fmpz_mod_poly_evaluate_fmpz(temp,tab[4],d);

// evaluation_fonction(resutltat, D, tab,ptoi, 6);
// 			calcul_poly_syndrome(poly_sy, resu,tab,ptoi,6);
// 			printf(" ressssssssu\n");

// 			fmpz_mod_poly_print(poly_sy);
// 						printf(" ressssssssu\n");

// fmpz_print(resutltat);
// 						printf("JKKKKJK D \n");
// 						fmpz_mod_poly_print(D);


// fmpz_mod_poly_init(D,tmp);
// 						printf(" D \n");

// 	fmpz_mod_poly_set_coeff_ui(D, 0, 14);

// fmpz_mod_poly_print(D);
// 		fmpz_mod_poly_zero(D);

// printf(" ,;lkmlkmlkmlD \n");

// 	fmpz_mod_poly_set_coeff_ui(D, 0, 14);

// printf(" \n");
// printf(" \n");
// printf(" \n");
// printf(" \n");
// printf(" \n");

// fmpz_mod_poly_t Q1,R1, q, r;
// fmpz_mod_poly_init(Q1,tmp-1);
// fmpz_mod_poly_init(q,tmp-1);
// fmpz_mod_poly_init(r,tmp-1);
// fmpz_mod_poly_init(R1,tmp-1);

	
// 	/*fmpz_mod_poly_set_coeff_ui(Q1, 4, 2); 
// 	fmpz_mod_poly_set_coeff_ui(Q1, 3, 3); 
// 	fmpz_mod_poly_set_coeff_ui(Q1, 2, 6); 
// 	fmpz_mod_poly_set_coeff_ui(Q1, 1, 6);
// 	fmpz_mod_poly_set_coeff_ui(Q1, 0, 12);*/

// 	fmpz_mod_poly_set_coeff_ui(Q1, 5, 5); 
// 	fmpz_mod_poly_set_coeff_ui(Q1, 4, 2); 
// 	fmpz_mod_poly_set_coeff_ui(Q1, 2, 15); 
// 	fmpz_mod_poly_set_coeff_ui(Q1, 1, 11);
// 	fmpz_mod_poly_set_coeff_ui(Q1, 0, 12);



	
// 	fmpz_mod_poly_set_coeff_ui(R1, 3, 5); 
// 	fmpz_mod_poly_set_coeff_ui(R1, 2, 7); 
// 	fmpz_mod_poly_set_coeff_ui(R1, 1, 1);
// 	fmpz_mod_poly_set_coeff_ui(R1, 0, 2);

// 	fmpz_mod_poly_print(BI);
// printf(" \n");
// 	fmpz_mod_poly_print(AI);

// 	/*division(q, r, BI, AI, tab, ptoi, reduc);

// 	fmpz_t teste;
// 	fmpz_init(teste);
// 	fmpz_mod_poly_get_coeff_fmpz(teste,R1,15);


// 	printf(" \n Q \n");
// 	fmpz_mod_poly_print(q);
// 	printf(" \n r \n");
// 	fmpz_mod_poly_print(r);
// printf(" \n");
// 	fmpz_print(teste);*/
	
// fmpz_mod_poly_init(q,tmp-1);
// fmpz_mod_poly_init(r,tmp-1);


// printf(" \n EUCLIDE AI\n");
// 	fmpz_mod_poly_print(AI);

// algo_euclide( r,q , BI, AI,6,tab, ptoi, reduc);












// /*fmpz_t teste;
// 	fmpz_init(teste);
// 	for(int i=1;i<=15;i++){
// evaluation_fonction(teste, r, tab, ptoi, i);

// 	printf(" \n %d   ",i);
// 	fmpz_print(teste);


// 	}
// //evaluation_fonction(teste, r, tab, ptoi, 6);
// printf("\n amplitude \n");
// 	fmpz_mod_poly_print(q);
// 	printf(" \n localisation\n");
// 	fmpz_mod_poly_print(r);

// 	printf(" \n \n");
// 	//fmpz_print(teste);*/





// fmpz_mod_poly_t tr, mu;
// fmpz_mod_poly_init(tr, tmp-1);
// fmpz_mod_poly_init(mu, tmp-1);
// fmpz_mod_poly_init(r, tmp-1);






// /*fmpz_mod_poly_set_coeff_ui(mu,0,12); 
// fmpz_mod_poly_set_coeff_ui(mu,1,3);
// fmpz_mod_poly_set_coeff_ui(mu,2,10);
// fmpz_mod_poly_set_coeff_ui(mu,3,9);
// fmpz_mod_poly_set_coeff_ui(mu,4,0);
// fmpz_mod_poly_set_coeff_ui(mu,5,13);
// fmpz_mod_poly_set_coeff_ui(mu,6,11);
// fmpz_mod_poly_set_coeff_ui(mu,7,10);
// fmpz_mod_poly_set_coeff_ui(mu,8,9);
// fmpz_mod_poly_set_coeff_ui(mu,9,5);
// fmpz_mod_poly_set_coeff_ui(mu,10,7);
// fmpz_mod_poly_set_coeff_ui(mu,11,6);
// fmpz_mod_poly_set_coeff_ui(mu,12,15);
// fmpz_mod_poly_set_coeff_ui(mu,13,4);
// fmpz_mod_poly_set_coeff_ui(mu,14,3);*/



// fmpz_mod_poly_init(mu, tmp-1);
// fmpz_mod_poly_set_coeff_ui(mu,1,1);
// fmpz_mod_poly_set_coeff_ui(mu,2,2);
// fmpz_mod_poly_set_coeff_ui(mu,3,3);
// fmpz_mod_poly_set_coeff_ui(mu,4,4);
// fmpz_mod_poly_set_coeff_ui(mu,5,5);
 
// /*evaluation_fonction(teste, r, tab, ptoi, 9);

// 	printf(" \n  ttt    \n");
// 	fmpz_print(teste);
// 	evaluation_fonction(teste, q, tab, ptoi, 9);

// 	printf(" \n  ttt    \n");
// 	fmpz_print(teste);*/

// fmpz_mod_poly_init(resu, tmp-1);
// 	 decode(resu,mu,tab,ptoi, 16,6,reduc );

// 	 printf(" \n  DATA\n");

// 	fmpz_mod_poly_print(resu);
// fmpz_mod_poly_init(resu, tmp-1);

// printf(" \n  \n \n");
// calcul_poly_syndrome( resu, mu,tab, ptoi, 6);





	return 0;
}