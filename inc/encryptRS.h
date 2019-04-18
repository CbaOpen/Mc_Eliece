#include <gmp.h>
#include <mpfr.h>
#include <flint/flint.h>
#include <stdlib.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/flintxx.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <mpf2mpfr.h>
#include <flint/fmpz_mod_poly.h>

#ifndef ENCRYPTRS_H
#define ENCRYPTRS_H

void lookuptab(fmpz_mod_poly_t reduc, fmpz_mod_poly_t* tab,int primitive, int length,fmpz_t* ptoi);
void findB(fmpz_mod_poly_t res,fmpz_mod_poly_t A,fmpz_mod_poly_t G,fmpz_t* ptoi,fmpz_mod_poly_t* tab);
void gen_poly(fmpz_mod_poly_t G,fmpz_mod_poly_t* tab,fmpz_t* ptoi,int tt,int tcycle);
void setBinPoly(fmpz_mod_poly_t res,fmpz_t f);
void mulPoly(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2, fmpz_mod_poly_t reduc);
void test_encryptRS();

#endif