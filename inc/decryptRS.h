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

#ifndef DECRYPTRS_H
#define DECRYPTRS_H


void evaluation_fonction(fmpz_t res, fmpz_mod_poly_t fonction, long valeur_x);
void xor_polynome(fmpz_mod_poly_t resu, fmpz_mod_poly_t poly1_v, fmpz_mod_poly_t poly2_v);
void division(fmpz_mod_poly_t q, fmpz_mod_poly_t r, fmpz_mod_poly_t dividente, fmpz_mod_poly_t deviseur);
void algo_euclide( fmpz_mod_poly_t localisation,fmpz_mod_poly_t amplitude , fmpz_mod_poly_t x2t_v, fmpz_mod_poly_t syndrome_v,long tt);
bool calcul_poly_syndrome(fmpz_mod_poly_t syndrome, fmpz_mod_poly_t data,long tt);
void derivation(fmpz_mod_poly_t res, fmpz_mod_poly_t function);
void decode(fmpz_mod_poly_t data, fmpz_mod_poly_t received ,long tt );
#endif