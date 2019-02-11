#include "../../inc/ReedSolomon/encryptRS.h"
#include <stdlib.h>
#include <gmp.h>
#include "/usr/local/include/flint/flint.h"
#include "/usr/local/include/flint/fmpz_mod_poly.h" // A REVOIR comment éviter de mettre tout le chemin?

fmpz_poly_t p_irred, p_gen, p_b, p_data, p_transm,p_temp,p_temp_mult,p_alpha;
int alpha=2;

fmpz_poly_init(p_irred);
fmpz_poly_init(p_gen);
fmpz_poly_init(p_b);
fmpz_poly_init(p_data);
fmpz_poly_init(p_transm);
fmpz_poly_init(p_temp);
fmpz_poly_init(p_temp_mult);
fmpz_poly_init(p_alpha);

fmpz_poly_set_coeff_ui(p_irred, 4, 1);
fmpz_poly_set_coeff_ui(p_irred, 1,  1);
fmpz_poly_set_coeff_ui(p_irred, 0,  1);

fmpz_poly_set_coeff_ui(p_gen, 0, 1);
fmpz_poly_set_coeff_ui(p_temp, 1, 1);

fmpz_poly_set_coeff_ui(p_temp_mult, 1, 1);

void generateG(fmpz_mod_poly_t res){

	for (int i =1;i<=t; i++){
		fmpz_poly_set_coeff_ui(p_alpha, 0,  temp);
		fmpz_poly_add( p_temp, p_temp ,p_alpha ); // (x-a^i)
		

		temp=*temp;

	}

	fmpz_poly_get_nmod_poly

}

void findB(fmpz_mod_poly_t res,fmpz_mod_poly_t Generator){

}





main(int argc, char** argv){
	printf("toto");

}