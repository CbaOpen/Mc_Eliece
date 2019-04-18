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

void findC(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2, fmpz_t n);
void evaluation_fonction(fmpz_t res, fmpz_mod_poly_t fonction, int valeur_x);
void xor_polynome(fmpz_mod_poly_t resu, fmpz_mod_poly_t poly1_v, fmpz_mod_poly_t poly2_v);
void division(fmpz_mod_poly_t q, fmpz_mod_poly_t r, fmpz_mod_poly_t dividente, fmpz_mod_poly_t deviseur);
void algo_euclide( fmpz_mod_poly_t localisation,fmpz_mod_poly_t amplitude , fmpz_mod_poly_t x2t_v, fmpz_mod_poly_t syndrome_v,int tt);
bool calcul_poly_syndrome(fmpz_mod_poly_t syndrome, fmpz_mod_poly_t data,int tt);
void derivation(fmpz_mod_poly_t res, fmpz_mod_poly_t function);
int decode(fmpz_mod_poly_t data, fmpz_mod_poly_t received , int nn, int tt );
void test_decode();
#endif