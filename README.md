# Cryptosystème de McEliece basé sur les codes de Reed-Solomon
Projet réalisé par Claire Baskevitch, Ziadath Babaedjou et Tristan Bessac - Mai 2019

Le cryptosystème de McEliece est un chiffrement à clé publique utilisant les codes correcteurs d'erreurs, ici les codes de Reed-Solomon, pour sécuriser un message. Ce programme implémente ce cryptosystème à l'aide de la librairie Flint.

Librairies nécessaires
================================================================

Flint est une librairie pour le langage C qui permet de faire de la théorie des nombres. Vous trouverez une archive pour son installation ici : http://www.flintlib.org/downloads.html.
Les bibliothèques GMP (version 5.1.1 minimum) et MPFR (version 3.0.0 minimum) sont également nécéssaires à l'installation et l'utilisation de Flint.

 Commandes
 ===============================================================

### - Compilation

    make compil
  
   L'exécutable se trouve dans le dossier bin/
### - Génération des clés

    mceliece -key n k

n et k étant les paramètres de sécurité.

les fichier KEYpub et KEYpriv contenant respectivement la clé publique et la clé privé sont créés dans le dossier courant.

### - Chiffrement

    mceliece -c KEYpub n k [msg] [chiffré]

n , k : paramètres de sécurité.

KEYpub : chemin vers la clé publique.

[msg] chemin vers le message à chiffrer.

[chiffré] chemin vers le fichier où sera inscrit le message chiffré. 

### - Déchiffrement

    mceliece -d KEYpriv n k [chiffre] [msg]

n , k : paramètres de sécurité.    

KEYpiv : chemin vers la clé privé.

[chiffre] : chemin vers le message chiffré.

[msg] : chemin vers le fichier où sera inscrit le message à déchiffrer.
