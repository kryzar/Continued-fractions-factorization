La fonction main prend un argument parmi la liste suivante (à passer 
en ligne de commande) :

F7_lp_eas : trouver un facteur de F7 avec la large prime variation
            et la early abort strategy (30 s environ) ;
F7_lp     : trouver un facteur de F7 avec la large prime variation
            (100 s environ) ;
F7_eas    : trouver un facteur de F7 avec la early abort strategy
            (60 s environ) ;
F7        : trouver un facteur de F7 sans variante (300 s environ).

Par exemple, pour factoriser F7 :
$ ./a.out F7

Si l'on ne donne pas l'un des arguments ci-dessus, on peut donner en argument
un nombre entier x.
Si x <= 150 : trouver un facteur d'un nombre aléatoire de x bits produit de
			  deux premiers (avec la early abort strategy et la large prime
			  variation) ;
autrement   : trouver un facteur de x (avec la early abort strategy et la
			  large prime variation).

Par exemple, pour trouver un facteur de 5394826801 :
$ ./a.out 5394826801
Par exemple, pour factoriser un nombre aléatoire de 100 bits :
$ ./a.out 100

Attention : la fonction main ne vérifie pas que l'argument est de la 
bonne forme. 
