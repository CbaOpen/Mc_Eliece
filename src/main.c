#include <flint/fmpz.h>
#include <stdio.h>
#include <flint/fmpz_poly.h>
#include <math.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod_poly.h>
#include "../inc/encryptRS.h"
#include "../inc/decryptRS.h"
#include "../inc/mceliece.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
	

fmpz_mod_poly_t reduc;
fmpz_mod_poly_t* tab;
fmpz_t* ptoi;
fmpz_t DEUX;    
long gf;					


/* convertit un entier en base 8 en entier en base 10 */
long convertOctalToDecimal(long octalNumber)
{
    long decimalNumber = 0, i = 0;

    while(octalNumber != 0)
    {
        decimalNumber += (octalNumber%10) * pow(8,i);
        ++i;
        octalNumber/=10;
    }

    i = 1;

    return decimalNumber;
}


/*initialise la mémoire des tables */
void init_tabs(long m){
	fmpz_init_set_ui(DEUX,2);
	double m2 =pow(2,m);
	gf=m2;	
	tab=malloc(sizeof(fmpz_mod_poly_t)*(long long)m2);
	for(long i=0;i<m2;i++){
		fmpz_mod_poly_init(tab[i],DEUX);
	}
	ptoi=malloc(sizeof(fmpz_t)*(long long)m2);
	if(ptoi==NULL || tab==NULL){printf("problème dans l'initialisation des tableaux");exit(0);}
	for(long i=0;i<(long)m2;i++){
		fmpz_init(ptoi[i]);
	
	}
	fmpz_set_ui(ptoi[1],0);

}

/* libère la mémoire des tables */
void clear_tabs(){
	for(long i=0;i<gf;i++){
		fmpz_mod_poly_clear(tab[i]);
		fmpz_clear(ptoi[i]);
	}
	fmpz_clear(DEUX);
	fmpz_mod_poly_clear(reduc);
	free(tab);
	free(ptoi);
}

/* récupère le polynome irreductible */
void get_reduc(long m){
	FILE* fd = fopen("poly_irr.txt", "r");
	if(fd==NULL){printf("fichier poly_irr.txt introuvable");exit(0);}
	long tmp;
	long tmp2;
	long test;
	fmpz_t mf;
	test=fscanf(fd,"%ld %ld",&tmp,&tmp2);
	while (tmp!=m && test){
		test=fscanf(fd,"%ld %ld",&tmp,&tmp2);
	}

	fmpz_init_set_ui(mf,m);
	fmpz_mod_poly_init(reduc,mf);
	tmp2=convertOctalToDecimal(tmp2);
	fmpz_set_ui(mf,tmp2);
	setBinPoly(reduc,mf);
	fmpz_clear(mf);
	fclose(fd);
	 
}

/* récupère le message à chiffré */
void Get_msg(FILE* M,int *msg,long k){
	long tmp;
	for(long i=0;i<k;i++){
		tmp=fgetc(M);
		if(tmp==EOF){	
			printf("le message ne contient pas le bon nombre de caractère\n");
			exit(1);
		}
		msg[i]=tmp;

	}
}

/* affich l'aide */
void help(){
	printf("les paramètres de sécurité sont n et k, n - k doit être pair\n");
	printf("Pour générer les clés : -key n k\n un fichier KEYpub contenant la clé public et un fichier KEYpriv contenant la clé privé sont générés\n\n");
	printf("Pour chiffré : -c KEYpub n k [msg] [chiffe]\n");
	printf("KEYpub, chemin vers le fichier de la clé publique\n");
	printf("[msg], chemin vers le fichier contenant le message à chiffrer\n");
	printf("[chiffre], chemin vers le fichier dans lequel sera écrit le message chiffré\n\n");
	printf("pour déchiffrer : -d KEYpriv n k [chiffre] [msg]\n");
	printf("KEYpriv, chemin vers le fichier contenant la clé privé\n");
	printf("[chiffre], chemin vers le fichier contenant le message chiffré\n");
	printf("[msg], chemin vers le fichier dans lequel sera écrit le message déchiffré\n\n");
	printf("n et k sont les paramètres de securité\n");

	exit(1);
}


/* analyse des temps, non utilisé par l'application */
void analyse_RS(){
	clock_t enc;
	clock_t dec;
	fmpz_t gff;

	long n = pow(2,8)-1;
	fmpz_mod_poly_t ms,c,res;
	long m=log(n+1)/log(2);
	get_reduc(m);
	init_tabs(m);
	get_reduc(m);
	lookuptab();
	fmpz_init_set_ui(gff,gf);
	fmpz_mod_poly_init(ms,gff);
	fmpz_mod_poly_init(res,gff);
	FILE* ResRS=fopen("ResRS2_t","w");
	for(long i=0;i<85;i++){
		fmpz_mod_poly_set_coeff_ui(ms,i,67);
	}
	for(long k=1;k<86;k+=1){
		if((n-k)%2==0){
			fmpz_mod_poly_init(c,gff);
			long t=(n-85)/2;
			printf("k %ld\n",k);
			//deb=clock();
			printf("encode\n");
			encrypt(c,ms,t,n);
			for(long i=0;i<k;i++){
				fmpz_mod_poly_set_coeff_ui(c,i,1);
			}
			enc=clock();
			printf("decode\n");
			decode(res,c,2*t);	
			dec=clock();
			fprintf(ResRS, "%ld %ld %LF \n",n,k,((long double)((dec-enc))));
			fmpz_mod_poly_clear(c);
			}
		}

}


/* analyse des temps, non utilisé par  l'application */
void analyse_mce(){
	    long n=pow(2,12)-1;	
		printf("n %ld\n",n);
		clock_t deb;
		clock_t gen;
		clock_t chiffr;
		clock_t dechiffr;
		FILE* ressim=fopen("ResSimuln_2_12.data","w");
		for(long k=1;k<n-1;k+=1000){
			if((n-k)%2==0){
				printf("k: %ld\n",k);
				fmpz_mat_t key;
				fmpz_mat_init(key,k,n);		
				deb=clock();
				keygen(key,n,k);
				gen=clock();
				fmpz_mat_clear(key);
				clear_tabs();
				//chiffre
				printf("chiffre\n");
				FILE * keyf=fopen("KEYpub","r");
				long t=(n-k)/2;
				fmpz_mat_init(key,k,n);
				fmpz_mat_fread(keyf,key);
				fclose(keyf);
				int *msg=malloc(sizeof(int)*k);
				for(long v=0;v<k;v++){
					msg[v]=67;
				} 
				FILE* c=fopen("chiffre","w");
				long m=log(n+1)/log(2);
				get_reduc(m);
				init_tabs(m);
				get_reduc(m);
				lookuptab();
				Encrypt_McEliece(c,msg,key,t,k,n);
				fclose(c);	
				clear_tabs();
				free(msg);
				fmpz_mat_clear(key);
				chiffr=clock();

				//dechiffr
				printf("dechiffre\n");	
				FILE* keys=fopen("KEYpriv","r");
				c=fopen("chiffre","r");
				FILE* mesg=fopen("msg","w");
				t=(n-k)/2;
				m=log(n+1)/log(2);
				get_reduc(m);
				init_tabs(m);
				get_reduc(m);
				lookuptab();
				Decrypt_McEliece(c,n,k,t,keys,mesg);
				fclose(keys);
				fclose(c);
				fclose(mesg);
				clear_tabs();
				dechiffr=clock();
				fprintf(ressim, "%ld %ld %LF %Lf %LF\n",n,k,((long double)((gen-deb)/CLOCKS_PER_SEC)),((long double)((chiffr-gen)/CLOCKS_PER_SEC)),((long double)((dechiffr-chiffr)/CLOCKS_PER_SEC)));	
			}
		}

}

	


int main(int argc, char** argv){
	if(argc<2){
	 help(); 
	}
	//mceliece -key n k 
	if(strcmp(argv[1],"-key")==0){
		if(argc!=4){help();}
		long n=	atol(argv[2]);
		long k=atol(argv[3]);
		if( (n-k)%2 != 0){
			printf("paramètre n et k incompatible \n");
			exit(1);
		}
		fmpz_mat_t key;
		fmpz_mat_init(key,k,n);
		keygen(key,n,k);
		fmpz_mat_clear(key);
		clear_tabs();
		return 0;

	}
	//mceliece -c KEYpub n k msg chiffré 
	if(strcmp(argv[1],"-c")==0){
		if(argc!=7){printf("nani");help();}
		FILE * keyf=fopen(argv[2],"r");
		if(keyf==NULL){printf("erreur lors de l'ouverture du fichier contenant la clé\n");exit(1);}
		fmpz_mat_t key;
		long n=atol(argv[3]);
		long k =atol(argv[4]);
		if( (n-k)%2 != 0){
			printf("paramètre n et k incompatible \n");
			exit(1);
		}
		long t=(n-k)/2;
		fmpz_mat_init(key,k,n);
		fmpz_mat_fread(keyf,key);
		fclose(keyf);
		FILE* M=fopen(argv[5],"r");
		if(M==NULL){printf("erreur lors de l'ouverture du fichier message\n");exit(1);}
		int *msg=malloc(sizeof(long)*k);
		Get_msg(M,msg,k);
		
		FILE* c=fopen(argv[6],"w");
		if(c==NULL){printf("erreur lors de l'ouverture du fichier chiffré\n");exit(1);}
		long m=log(n+1)/log(2);
		get_reduc(m);
		init_tabs(m);
		get_reduc(m);
		lookuptab();
		Encrypt_McEliece(c,msg,key,t,k,n);
		fclose(c);	
		fclose(M);
		clear_tabs();
		free(msg);
		fmpz_mat_clear(key);
		return 0;
	}

	//mceliece -d KEYpriv n k chiffre msg
	if(strcmp(argv[1],"-d")==0){
		if(argc!=7){help();}
		FILE* keys=fopen(argv[2],"r");
		if(keys==NULL){printf("erreur lors de l'ouverture du fichier contenant la clé\n");exit(1);}
		FILE* c=fopen(argv[5],"r");
		if(c==NULL){printf("erreur lors de l'ouverture du fichier chiffré\n");exit(1);}
		FILE* msg=fopen(argv[6],"w");
		if(msg==NULL){printf("erreur lors de l'ouverture du message\n");exit(1);}
		long n=atol(argv[3]);
		long k=atol(argv[4]);
		if( (n-k)%2 != 0){
			printf("paramètre n et k incompatible \n");
			exit(1);
		}
		long t=(n-k)/2;
		long m=log(n+1)/log(2);
		get_reduc(m);
		init_tabs(m);
		get_reduc(m);
		lookuptab();
		Decrypt_McEliece(c,n,k,t,keys,msg);
		fclose(keys);
		fclose(c);
		fclose(msg);
		clear_tabs();
		return 0; 
	}
	help();
	return 1;


}






