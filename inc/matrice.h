#ifndef MATRICE_H
#define MATRICE_H


void RPermut(fmpz_mat_t P,int n);
void get_S(fmpz_mat_t S,int k);
void get_G(fmpz_mat_t Gm,fmpz_mod_poly_t G,int n,int k);
void multGFmat(fmpz_mat_t res,fmpz_mat_t op1,fmpz_mat_t op2,int n,int p,int m);
void keygen(fmpz_mat_t key,int n,int k);
void Encrypt_McEliece(FILE* c,int* m,fmpz_mat_t key,int t,int k,int n);
void Decrypt_McEliece(FILE* chiffre,int n,int k,int t,FILE* keys,FILE* msg);
#endif