#ifndef MCELIECE_H
#define MCELIECE_H


void RPermut(fmpz_mat_t P,long n);
void get_S(fmpz_mat_t S,long k);
void get_G(fmpz_mat_t Gm,fmpz_mod_poly_t G,long n,long k);
void multGFmat(fmpz_mat_t res,fmpz_mat_t op1,fmpz_mat_t op2,long n,long p,long m);
void keygen(fmpz_mat_t key,long n,long k);
void Encrypt_McEliece(FILE* c,int* m,fmpz_mat_t key,long t,long k,long n);
void Decrypt_McEliece(FILE* chiffre,long n,long k,long t,FILE* keys,FILE* msg);
#endif