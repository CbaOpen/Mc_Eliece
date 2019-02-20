#include "../../inc/ReedSolomon/encryptRS.h"
#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
 #include <stdio.h>
#include <math.h>
#include <flint/fq.h>
#include<flint/fq_poly.h> 

#include <flint/fmpz_mod_poly_factor.h>
//#include <flint/fmpz_mod_poly_matxx.h>
#include <flint/flint.h>
 #include <flint/fmpz.h>
 #include <flint/fmpz_poly.h>
//#include <flint/fmpz_mod_polyxx.h>
 #include <flint/fmpz_mat.h>

#include <mpf2mpfr.h>
 #include <flint/fmpz_mod_poly.h>

// uy	.h          fmpq_mat.h          fmpz_matxx.h              fmpz_polyxx.h          fq_poly_factor.h            gmpcompat.h         nmod_poly_mat.h    perm.h
// arithxx.h        fmpq_matxx.h        fmpz_mod_poly_factor.h    fmpz_vec.h             fq_poly_factor_templates.h  long_extras.h       nmod_poly_matxx.h  permxx.h
// config.h         fmpq_poly.h         fmpz_mod_poly_factorxx.h  fmpz_vecxx.h           fq_poly.h                   longlong.h          nmod_polyxx.h      profiler.h
// d_mat.h          fmpq_polyxx.h       fmpz_mod_poly.h           fmpzxx.h               fq_poly_templates.h         mpf_mat.h           nmod_vec.h         qadic.h
// double_extras.h  fmpq_vec.h          fmpz_mod_polyxx.h         fq.h                   fq_templates.h              mpfr_mat.h          nmod_vecxx.h       qadicxx.h
// d_vec.h          fmpqxx.h            fmpz_poly_factor.h        fq_mat.h               fq_vec.h                    mpfr_vec.h          NTL-interface.h    qsieve.h
// fft.h            fmpz-conversions.h  fmpz_poly_factorxx.h      fq_mat_templates.h     fq_vec_templates.h          mpf_vec.h           padic.h            templates.h
// fft_tuning.h     fmpz_factor.h       fmpz_poly.h               fq_nmod.h              fq_zech.h                   mpn_extras.h        padic_mat.h        ulong_extras.h
// flint.h          fmpz_factorxx.h     fmpz_poly_mat.h           fq_nmod_mat.h          fq_zech_mat.h               nmod_mat.h          padic_matxx.h
// flintxx          fmpz.h              fmpz_poly_matxx.h         fq_nmod_poly_factor.h  fq_zech_poly_factor.h       nmod_matxx.h        padic_poly.h
// flintxx.h        fmpz_lll.h          fmpz_poly_q.h             fq_nmod_poly.h         fq_zech_poly.h              nmod_poly_factor.h  padic_polyxx.h
// fmpq.h           fmpz_mat.h          fmpz_poly_qxx.h           fq_nmod_vec.h          fq_zech_vec.h               nmod_poly.h         padicxx.h



//void generate 

void lookuptab(fmpz_mod_poly_t reduc, fmpz_mod_poly_t* tab,int primitive, int length,fmpz_t* ptoi  ){
		fmpz_mod_poly_t poly;
		fmpz_t d,ev;
		fmpz_t tmp;
		fmpz_init(tmp);
		fmpz_init(ev);
		fmpz_init_set_ui(d,2);
			
		fmpz_mod_poly_set_coeff_ui(tab[0],0,1);

		fmpz_mod_poly_set_coeff_ui(tab[1],1,1);
	for(int i=2;i<length;i++){
		//fmpz_mod_poly_set_coeff_ui(poly,1);
		fmpz_mod_poly_mulmod(tab[i] , tab[i-1] , tab[1],reduc);
		//tab[i]->p=16;
		//fmpz_mod_poly_powmod_ui_binexp(tab[i] , poly , i, reduc)
	}
	fmpz_set_ui(ptoi[0],0);
	for(int i=1;i<length;i++){
		tab[i]->p=16;
		fmpz_mod_poly_evaluate_fmpz(ev,tab[i],d);
		//fmpz_print(ev);
		//printf("\n");
		fmpz_set_ui(tmp,i);
		//fmpz_print(tmp);
		//printf("tabi%d\n",(int)fmpz_get_d(ev));
		fmpz_set(ptoi[fmpz_get_ui(ev)],tmp);
		//fmpz_print(ptoi[(int)fmpz_get_d(ev)]);
		//printf("ev\n");
		//fmpz_mod_poly_set_coeff_ui(poly,1);
		//fmpz_mod_poly_mulmod(tab[i] , tab[i-1] , tab[1],reduc);
		//fmpz_mod_poly_powmod_ui_binexp(tab[i] , poly , i, reduc)
	}



}



void generateG(fmpz_mod_poly_t G,fmpz_mod_poly_t reduc,fmpz_mod_poly_t* tab,fmpz_mod_poly_t red){
	fmpz_t n,n2;
	fmpz_t d;
	fmpz_init_set_ui(d,2);
	fmpz_t ev;
	fmpz_init_set_ui(ev,8);
	fmpz_init_set_ui(n,16);
	fmpz_init_set_ui(n2,16);
	fmpz_t t,t2;
	fmpz_init_set_ui(t,2);
	fmpz_mod_poly_t p;
	fmpz_mod_poly_init(p,n);
	//fmpz_mod_poly_set_coeff_ui(p,0,1);
	//fmpz_mod_poly_set_coeff_ui(p,3,1);
	fmpz_t a;
	//fmpz_init_set_ui(a,2);
	//fmpz_poly_t G;
	//fmpz_mod_poly_init(G, n);
	fmpz_mod_poly_set_coeff_ui(G,1,1);
	fmpz_mod_poly_set_coeff_ui(G,0,2);

	fmpz_mod_poly_t tmp;
	//fmpz_mod_poly_init(templates,n);
	for (int i=2;i<7;i++){			
		fmpz_mod_poly_init(tmp,n);

		fmpz_mod_poly_evaluate_fmpz(a,tab[i],d);
		//fmpz_neg(a,a);
		fmpz_print(a);
		printf("a\n");
		fmpz_mod_poly_set_coeff_ui(tmp,1,1);
		fmpz_mod_poly_set_coeff_fmpz(tmp,0,a);
		printf("%d:\n",i );
		fmpz_mod_poly_print(tmp);
		printf( "\n");
		fmpz_mod_poly_mul(G , G , tmp);



		fmpz_mod_poly_clear(tmp);

	}
	fmpz_poly_t conv;
	fmpz_poly_init(conv);
	fmpz_mod_poly_t G2,G3;
	fmpz_mod_poly_init(G2,d);
	fmpz_mod_poly_init(G3,d);
	fmpz_t r;
	fmpz_init(r);
	fmpz_mod_poly_get_coeff_fmpz(r, G   , 3	);
	fmpz_print(r);
	for(size_t bi=fmpz_sizeinbase(r,2);bi>0;bi--){
		fmpz_mod_poly_set_coeff_ui(G2,bi,fmpz_tstbit(r, bi));
	}
	fmpz_mod_poly_set_coeff_ui(G2,0,fmpz_tstbit(r, 0));
	//fmpz_init_set_ui(r,fmpz_mod_poly_lead(G)[0]);
	 //fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, constfmpz_poly_t g)
	//fmpz_mod_poly_set_fmpz(G2,r);	

	//fmpz_mod_poly_set(p,G);	
	 //fmpz_mod_poly_get_fmpz_poly(G3,G3);
	
	//fmpz_mod_poly_get_fmpz_poly( conv,G3);
	//fmpz_mod_poly_set_fmpz_poly(G3,conv);

	fmpz_mod_poly_rem_basecase(G3, G2 , red);
    _fmpz_mod_poly_normalise(G);
	printf("G2\n");
	fmpz_mod_poly_print(G);
	printf("\n");

	fmpz_mod_poly_get_fmpz_poly( conv,G3);
	fmpz_mod_poly_set_fmpz_poly(G3,conv);
	printf("G3\n");
	fmpz_mod_poly_print(G3);
	printf("\n");
	//G->p=500000;
	fmpz_mod_poly_evaluate_fmpz(a,G,ev);
	fmpz_print(a);
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
		//fmpz_set(t,ptoi[fmpz_get_ui(t)]);
		fmpz_mod_poly_set_coeff_fmpz(res,i,ptoi[fmpz_get_ui(t)]);
	}
		for(int i=0;i<fmpz_mod_poly_degree(G)+1;i++){
		fmpz_mod_poly_get_coeff_fmpz(t,G,i);
		//fmpz_set(t,ptoi[fmpz_get_ui(t)]);
		fmpz_mod_poly_set_coeff_fmpz(Gtemp,i,ptoi[fmpz_get_ui(t)]);
	}


	//boucle
	while(fmpz_mod_poly_degree(res)>fmpz_mod_poly_degree(G)-1){
			printf(" %ld  ca tourne\n", fmpz_mod_poly_degree(res));
		
		//decalage puissance de X
		fmpz_mod_poly_set_coeff_ui(tmp,fmpz_mod_poly_degree(res)-fmpz_mod_poly_degree(G),1);
		fmpz_mod_poly_mul(tmp2,tmp,Gtemp);
		// printf("res\n");
		// fmpz_mod_poly_print(tmp2);

		//mult des coeff
		fmpz_mod_poly_get_coeff_fmpz(t,res,fmpz_mod_poly_degree(res));
		printf("\n t: ");fmpz_print(t);printf("\n");
		if(fmpz_get_ui(t)!=15){
			for(int i=0;i<fmpz_mod_poly_degree(tmp2)+1;i++){
				fmpz_mod_poly_get_coeff_fmpz(t2,tmp2,i);
				if(!fmpz_is_zero(t2)){
					fmpz_add(t2,t2,t);
					fmpz_mod(t2,t2,q);
					if(fmpz_get_ui(t2)==0)fmpz_set_ui(t2,15);
					fmpz_mod_poly_set_coeff_fmpz(tmp2,i,t2);	
				//printf("t2 : \n");fmpz_print(t2);
			}}
		}

		// printf("tmp2\n");
		// fmpz_mod_poly_print(tmp2);
	// 	//xor coeff
		for(int i=0;i<fmpz_mod_poly_degree(res)+1;i++){
			fmpz_mod_poly_get_coeff_fmpz(t,res,i);
			fmpz_mod_poly_get_coeff_fmpz(t2,tmp2,i);

			if(!fmpz_is_zero(t))fmpz_mod_poly_set(tmp3,tab[fmpz_get_ui(t)]);else fmpz_mod_poly_zero(tmp3);
			// printf("\ntmp2\n");
			// fmpz_mod_poly_print(tmp2);
			// printf("\n");
			if(!fmpz_is_zero(t2))fmpz_mod_poly_set(tmp4,tab[fmpz_get_ui(t2)]);else fmpz_mod_poly_zero(tmp4);
			// printf("\ntmp3p\n");
			// fmpz_mod_poly_print(tmp3);
			// printf("\n");

			fmpz_mod_poly_add(tmp3,tmp3,tmp4);
			tmp3->p=16;
			tmp4->p=16;
			// printf("\ntmp3\n");
			// fmpz_mod_poly_print(tmp3);
			// printf("\n");

			fmpz_mod_poly_evaluate_fmpz(t,tmp3,deux);
			if(!fmpz_is_zero(t))fmpz_set(t,ptoi[fmpz_get_ui(t)]);

			fmpz_mod_poly_set_coeff_fmpz(res,i,t);

			
	    //pour chaque puissance xor des poly des alpha correspondant (=== addition mod 2) puis retrouver leur puissance avec ptoi
			//fmpz_xor(t,t,t2);

			//fmpz_mod_poly_set_coeff_fmpz(res,i,t);
			tmp3->p=2;
			tmp4->p=2;
		}
			fmpz_mod_poly_print(res);
		_fmpz_mod_poly_normalise(res);
		fmpz_mod_poly_zero(tmp);
	}

	for(int i=0;i<fmpz_mod_poly_degree(res)+1;i++){
		fmpz_mod_poly_get_coeff_fmpz(t,res,i);
		//fmpz_set(t,ptoi[fmpz_get_ui(t)]);
		fmpz_mod_poly_evaluate_fmpz(t,tab[fmpz_get_ui(t)],deux);
		fmpz_mod_poly_set_coeff_fmpz(res,i,t);
	}



}




// for i = 1 step 1 until 2t do
//      begin
//      for j = 2t step -1 until 1 do
//            begin
//            g[j] = mult(g[j], alpha[i]) XOR g[j-1]
//            end
//      g[0] = mult(g[0]), alpha[i])
//      end



// 1 -> 15
// mod 15 0->14
// 	 

void genG(fmpz_mod_poly_t G,fmpz_mod_poly_t* tab,fmpz_t* ptoi){
	fmpz_t tmp;
	fmpz_t ev;
	fmpz_t d;
	fmpz_t s;
	fmpz_init_set_ui(d,2);
	fmpz_init_set_ui(s,15);
	fmpz_mod_poly_set_coeff_ui(G,0,1);	
	fmpz_init(tmp);
	fmpz_init(ev);
	for (int i=1;i<7;i++){
		for (int j=6;j>0;j--){
			fmpz_mod_poly_evaluate_fmpz(ev,tab[i],d);
			fmpz_mod_poly_get_coeff_fmpz(tmp,G,j); 
			fmpz_mul(tmp,tmp,ev);
			//if (!fmpz_is_zero(tmp)){
			
				fmpz_mod(tmp,tmp,s);
				fmpz_set(tmp,ptoi[fmpz_get_ui(tmp)]);
				fmpz_mod_poly_get_coeff_fmpz(ev,G,j);
				fmpz_xor(tmp,tmp,ev);
			
			//}
			fmpz_mod_poly_set_coeff_fmpz(G,j,tmp);
		}

		fmpz_mod_poly_get_coeff_fmpz(tmp,G,0);
		fmpz_mod_poly_evaluate_fmpz(ev,tab[i],d);
		fmpz_mul(tmp,tmp,ev);
		//if (!fmpz_is_zero(tmp)){
			fmpz_mod(tmp,tmp,s);
			fmpz_set(tmp,ptoi[fmpz_get_ui(tmp)]);
		//}
		
		fmpz_mod_poly_set_coeff_fmpz(G,0,tmp);

printf("%d\n",i);
fmpz_mod_poly_print(G);
printf("\n");
	}
}


// void gen_poly()
// /* Obtain the generator polynomial of the tt-error correcting, length
//   nn=(2**mm -1) Reed Solomon code  from the product of (X+alpha**i), i=1..2*tt
// */
//  {
//    register int i,j ;

//    gg[0] = 2 ;    /* primitive element alpha = 2  for GF(2**mm)  */
//    gg[1] = 1 ;    /* g(x) = (X+alpha) initially */
//    for (i=2; i<=nn-kk; i++)
//     { gg[i] = 1 ;
//       for (j=i-1; j>0; j--)
//         if (gg[j] != 0)  gg[j] = gg[j-1]^ alpha_to[(index_of[gg[j]]+i)%nn] ;
//         else gg[j] = gg[j-1] ;
//       gg[0] = alpha_to[(index_of[gg[0]]+i)%nn] ;     /* gg[0] can never be zero */
//     }
//    /* convert gg[] to index form for quicker encoding */
//    for (i=0; i<=nn-kk; i++)  gg[i] = index_of[gg[i]] ;
//  }
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
					fmpz_print(tmp);
					fmpz_mod_poly_set_coeff_fmpz(G,j,tmp);

			}else{
				fmpz_mod_poly_get_coeff_fmpz(tmp,G,j-1);
				fmpz_mod_poly_set_coeff_fmpz(G,j,tmp);				
			}

		}
					fmpz_mod_poly_print(G);
					printf("\n");
		//       gg[0] = alpha_to[(index_of[gg[0]]+i)%nn] ; 
					fmpz_mod_poly_get_coeff_fmpz(tmp,G,0); //tmp -> g[j-1]
					//fmpz_mod_poly_get_coeff_fmpz(tmp2,G,j); // tmp2 -> g[j]
					fmpz_set(tmp,ptoi[fmpz_get_ui(tmp)]);
					fmpz_add_ui(tmp,tmp,i); // tmp2 +i  
					fmpz_mod(tmp,tmp,s);
					fmpz_mod_poly_evaluate_fmpz(ev,tab[fmpz_get_ui(tmp)],d);
					fmpz_mod_poly_set_coeff_fmpz(G,0,ev);					



	}
}


void g_poly(fmpz_mod_poly_t G,fmpz_mod_poly_t* tab,fmpz_t* ptoi){

}


// void mul_poly(fmpz_mod_poly_t f,fmpz_mod_poly_t p1,fmpz_mod_poly_t p2,fmpz_mod_poly_t* tab,fmpz_t* ptoi){
// 	fmpz_mod_poly_t* res=malloc(sizeof(fmpz_mod_poly_t)*fmpz_mod_poly_degree(p1)+fmpz_mod_poly_degree(p2)+1);
// 	fmpz_mod_poly_t* tmp=malloc(sizeof(fmpz_mod_poly_t)* fmpz_mod_poly_length(p1));
// 	fmpz_mod_poly_t* tmp2=malloc(sizeof(fmpz_mod_poly_t)* fmpz_mod_poly_length(p2));
// 	fmpz_t temp;
// 	fmpz_init(temp);

// 	for (int i=0;i<fmpz_mod_poly_degree(p1)+fmpz_mod_poly_degree(p2)+1;i++){
// 		fmpz_mod_poly_init(res[i]);
// 	}


// 	for (int i=0;i<fmpz_mod_poly_length(p1)+1;i++){
		 
// 		 fmpz_mod_poly_init(tmp[i]);
// 		 fmpz_mod_poly_get_coeff_fmpz(temp,p1,i);
// 		fmpz_mod_poly_set(tmp[i],tab[(int)fmpz_get_d(ptoi[(int)fmpz_get_d(temp)])]);
// 	}
// 	for (int i=0;i<fmpz_mod_poly_length(p2)+1;i++){
// 		 fmpz_mod_poly_init(tmp2[i]);
// 		 fmpz_mod_poly_get_coeff_fmpz(temp,p2,i);
// 		fmpz_mod_poly_set(tmp2[i],tab[(int)fmpz_get_d(ptoi[(int)fmpz_get_d(temp)])]);
// 	}

// 	for (int i=0;i<fmpz_mod_poly_length(p1)+1;i++){
// 		for (int i=0;i<fmpz_mod_poly_length(p2)+1;i++){
				
// 		}
// 	}



// }

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
	
	for(int i=0;i<fmpz_mod_poly_degree(op1)+1;i++){
		for(int j=0;j<fmpz_mod_poly_degree(op2)+1;j++){
			fmpz_mod_poly_get_coeff_fmpz(c1,op1,i);
			fmpz_mod_poly_get_coeff_fmpz(c2,op2,j);
			

		
			setBinPoly(ptmp,c1);
			setBinPoly(ptmp2,c2);
			// fmpz_mod_poly_print(ptmp);
			// printf("\n");
			// fmpz_mod_poly_print(ptmp2);
			// printf("\n");

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



}

void divPoly(fmpz_poly_t res,fmpz_poly_t A,fmpz_poly_t div){
	
	fmpz_t s;
	fmpz_init_set_ui(s,16);

	fmpz_mod_poly_t reduc;
	fmpz_mod_poly_init(reduc, s);
	fmpz_mod_poly_set_coeff_ui(reduc, 4, 1);
	fmpz_mod_poly_set_coeff_ui(reduc, 3, 1);
	fmpz_mod_poly_set_coeff_ui(reduc, 0, 1);




	fmpz_t dA;
	fmpz_t ddiv;
	fmpz_t tmp;
	fmpz_t cA;
	fmpz_t cdiv;
	fmpz_t q;
	fmpz_poly_t ccA;
	fmpz_poly_init(ccA);
	fmpz_poly_set(ccA,A);
	fmpz_poly_t tmp2;

	
	


	fmpz_mod_poly_t tmpm1;
	fmpz_mod_poly_t tmpm2;
	fmpz_mod_poly_t tmpm3;
	fmpz_mod_poly_init(tmpm1,s);
	fmpz_mod_poly_init(tmpm2,s);
	fmpz_mod_poly_init(tmpm3,s);
	//boucle
	while (fmpz_poly_degree(ccA)>=fmpz_poly_degree(div)){
		fmpz_init(dA);fmpz_init(ddiv);fmpz_init(tmp);fmpz_init(cA);fmpz_init(cdiv);
		fmpz_init_set_ui(dA,fmpz_poly_degree(ccA));
		fmpz_init_set_ui(ddiv,fmpz_poly_degree(div));
		fmpz_sub(tmp,dA,ddiv);
		fmpz_poly_get_coeff_fmpz(cA,ccA, fmpz_poly_degree(ccA));
		fmpz_poly_get_coeff_fmpz(cdiv,div, fmpz_poly_degree(div));

		fmpz_fdiv_q(q,cA,cdiv); 
		
		fmpz_poly_init(tmp2);
			
		fmpz_poly_set_coeff_fmpz(tmp2,fmpz_get_ui(tmp),q);

		fmpz_mod_poly_set_fmpz_poly(tmpm1,div);
		fmpz_mod_poly_set_fmpz_poly(tmpm2,tmp2);		
		mulPoly(tmpm3,tmpm1,tmpm2,reduc);

		fmpz_mod_poly_get_fmpz_poly( res,tmpm3);


		
		fmpz_poly_add(res,res,ccA);
		
		fmpz_poly_set_coeff_ui(res, fmpz_poly_degree(res),0);
		fmpz_clear(dA);fmpz_clear(ddiv);fmpz_clear(tmp);fmpz_clear(cA);fmpz_clear(cdiv);fmpz_poly_clear(tmp2);
		fmpz_poly_set(ccA,res);
		printf("hey\n");
	}
}








void main(int argc, char** argv){

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
	fmpz_t n2;
	fmpz_t tmp;
	fmpz_t tmp2;
	fmpz_init_set_ui(tmp2,5000000);
	

	fmpz_init_set_ui(tmp,16);
	fmpz_init_set_ui(n,4);

	fmpz_init(n2);
	fmpz_mul(n2,n,n);


	//fmpz_init_set_ui(n, 4);
	fmpz_mod_poly_init(reduc, n);
	fmpz_mod_poly_set_coeff_ui(reduc, 4, 1);
	fmpz_mod_poly_set_coeff_ui(reduc, 3, 1);
	fmpz_mod_poly_set_coeff_ui(reduc, 0, 1);
	fmpz_mod_poly_t G;
	fmpz_mod_poly_init(G,tmp);
	lookuptab(reduc,tab,2, 16,ptoi);
	//generateG(G,tab[4],tab,reduc);
	//genG(G,tab,ptoi);
	gen_poly(G,tab,ptoi,6,15);
	printf("g:\n");
	fmpz_mod_poly_print(G);
	printf("\n" );
	//fmpz_mod_poly_print(reduc);

	printf("tab:\n");
	for(int i=0;i<16;i++){
		printf("\n%d----",i);
		fmpz_mod_poly_print(tab[i]);
		printf("  alpha 	i: ");
		fmpz_print(ptoi[i]); 
	}


	fmpz_mod_poly_t B;
	fmpz_mod_poly_init(B,tmp);
	fmpz_mod_poly_t A;
	fmpz_mod_poly_init(A,tmp);
	fmpz_mod_poly_t R;
	fmpz_mod_poly_init(R,tmp);
	fmpz_mod_poly_set_coeff_ui(A, 8, 9);
	fmpz_mod_poly_set_coeff_ui(A, 7, 8);
	fmpz_mod_poly_set_coeff_ui(A, 6, 7);
	fmpz_mod_poly_set_coeff_ui(A, 5, 6);
	fmpz_mod_poly_set_coeff_ui(A, 4, 5);
	fmpz_mod_poly_set_coeff_ui(A, 3, 4);
	fmpz_mod_poly_set_coeff_ui(A, 2, 3);
	fmpz_mod_poly_set_coeff_ui(A, 1, 2);
	fmpz_mod_poly_set_coeff_ui(A, 0, 1);

	fmpz_mod_poly_t X;
	fmpz_mod_poly_init(X,tmp);
	fmpz_mod_poly_set_coeff_ui(X, 6, 1);
	mulPoly(B,A,X,reduc);

	//fmpz_mod_poly_set_fmpz_poly(G3,conv);
	//divPoly(X,B,conv);
	findB(R,B,G,ptoi,tab);
	//fmpz_mod_poly_get_fmpz_poly( conv,reduc);
	//fmpz_poly_rem_basecase(X,X,conv);
	//fmpz_mod_polyxxoly_mul(B,A,X);
	printf("\nB\n");
	fmpz_mod_poly_print(R);
	printf("\n");
	// fmpz_mod_poly_mul(A,A,X);
	// fmpz_mod_poly_add(A,A,B);
	// fmpz_mod_poly_rem_basecase(X,A,G);
	// fmpz_mod_poly_print(X);




//test mult poly
	fmpz_mod_poly_t res,op1,op2;
	fmpz_mod_poly_init(res,tmp);
	fmpz_mod_poly_init(op1,tmp);
	fmpz_mod_poly_init(op2,tmp);

	//fmpz_mod_poly_set_coeff_ui(op1, 1, 1);
	fmpz_mod_poly_set_coeff_ui(op1, 0, 5);
	//fmpz_mod_poly_set_coeff_ui(op2, 1, 1);
	fmpz_mod_poly_set_coeff_ui(op2, 0, 13);
	mulPoly(res,op1,op2,reduc);
	fmpz_mod_poly_print(res);







// fmpz_t a;

// fmpz_init(a);



// 	fq_ctx_t  GF16;

//  fq_ctx_init_modulus(GF16 , reduc , "gen");
// fq_poly_t Gt,tmpt;
// fq_poly_init(Gt ,  GF16);


// fq_poly_set_coeff(Gt,1,1,GF16);
// 	fq_poly_set_coeff(Gt,0,2,GF16);
// for (int i=2;i<7;i++){			
// 		fq_init(tmpt,GF16);

// 		fmpz_mod_poly_evaluate_fmpz(a,tab[i],deux);
// 		//fmpz_neg(a,a);
// 		fmpz_print(a);
// 		printf("a\n");
// 		fq_poly_set_coeff(tmpt,1,1,GF16);
// 		fq_poly_set_coeff_fmpz(tmpt,0,a);
// 		printf("%d:\n",i );
// 		fq_print(tmpt,GF16);
// 		printf( "\n");
// 		fq_mul(Gt , Gt , tmpt,GF16);



// 		fq_clear(tmpt,GF16);

// 	}
// 	fq_print(Gt,GF16);















 
}