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

void lookuptab();
void findB(fmpz_mod_poly_t res,fmpz_mod_poly_t A,fmpz_mod_poly_t G);
void gen_poly(fmpz_mod_poly_t G,long tt);
void setBinPoly(fmpz_mod_poly_t res,fmpz_t f);
void mulPoly(fmpz_mod_poly_t res,fmpz_mod_poly_t op1, fmpz_mod_poly_t op2);
void encrypt(fmpz_mod_poly_t res,	fmpz_mod_poly_t m, long t,long n);

#endif