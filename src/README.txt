La fonction main prend un argument parmi la liste suivante (à passer 
en ligne de commande) :

F7_lp_eas : pour trouver un facteur de F7 avec la large prime variation
            et la early abort strategy (30 s environ)
F7_lp     : pour trouver un facteur de F7 avec la large prime variation
            (100 s environ)
F7_eas    : pour trouver un facteur de F7 avec la early abort strategy
            (60 s environ)
F7        : pour trouver un facteur de F7 sans variante (300 s environ)

un entier nb_bits <= 150 : pour trouver un facteur d'un nombre aléatoire
                           de nb_bits produit de deux premiers (avec la
                           early abort strategy et la large prime variation)
un entier N en base 10   : pour trouver un facteur de N (avec la early
                           abort strategy et la large prime variation)

Attention : la fonction main ne vérifie pas que l'argument est de la 
bonne forme. 
