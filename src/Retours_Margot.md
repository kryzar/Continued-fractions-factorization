# Remarques générales
- Attention à la longueur de tes lignes ! C'est pénible les sauts indigestes
  faits automatiquement. Les faire à la main assure que tout le monde aura le
  même visuel partout et que ce visuel est maitrisé.
- Fais attention aussi aux alignements, regarde comment je fais. Je pense
  notamment aux opérateurs binaires comme l'opérateur =. Regarde aussi les
  signatures de fonctions sur plusieurs lignes.
- Vois les choses comme ça : une fonction = une tâche. Si tu fais deux choses
  différentes au sein d'une fonction (sauf les grosses fonctions), tu
  compliques l'architecture.
- Noms des macros en majuscules.
- Dans une boucle for t'es pas censée déclarer l'incrément en dehors. Faut
  faire for (size_t i = 0, …).
- Pour les accolades après parenthèses, ) { est mieux que ){. Idem, } else {
  est mieux que }else{ (même si je crois que la convention C est pas d'accord).
- Les cadres faits avec des /* */ sont faits pour des titres, non des longues
  explications. Dans ce cas, tu fais un cadre de tire et tu commentes
  l'explication dessous.
- Pas de guillemets autour des noms de variables !
- Sois plus concise sur les commentaires inline.
- Dans les signatures de fonctions, vaut mieux dire "Create …" que "This
  functions create …".
- Pour les struct, c'est pas mieux de mettre -> (sans espaces) plutôt que  ->
  (avec espaces) ? Si t'es d'accord que c'est mieux sans espaces, change.

# main.c
1. Créer une fonction find_factor qui contient les allocations ; les appels aux
   étapes A, B, B ; la libération ; et plus généralement tout l'algorithme. Ce
   n'est pas à main de faire ça (même si tu précises une section genre "Looking
   for a factor"). Tu peux donner en argument de cette fonction un booléen lpv
   qui indique si oui ou non tu fais la large prime variation. Dans ce cas,
   plutôt que de commenter et décommenter, tu mets tout dans une boucle if.
2. Créer une fonction print_results qui prend en argument les résultats et
   affiche tout, comme dans la deuxième section "Testing" de ton main. Tu peux
   éventuellement rajouter un booléen printing dans la fonction find_factor qui
   permet d'appeler ou non print_results.
3. Mettre les tests dans des fonctions et fichiers dédiés.
4. Commenter la fonction main.

# step_A.h
1. Ou bien tu écris ta macro en majuscules, ou bien tu la remplaces par une
   fonction inline, qui est sûrement mieux. Voir
   https://www.geeksforgeeks.org/inline-function-in-c/.
1. Documenter Params. Pourquoi on la créée, à quoi sert-elle ? La description
   de ses composantes ne suffit pas. C'est comme si tu analysais l'architecture
   d'une maison en ne parlant que de ses murs individuels.

# step_A.c
1. Documenter init_mpz_array et free_mpz_array ? Pas sûr que ça vaille le coup.

# step_B.h
1. Choisis entre déclarer la struct et le typedef séparément ou en même temps
   et modifie Data_exp_vect et Params en conséquence.
2. Renomme Data_exp_vect en exp_vect_data.

# step_B.c
1. Peut-être un mot sur ce que fait mp_bitcnt_t.
2. Remplace is_qn_factorisable par is_Qn_factorisable.
3. Dans is_qn_factorisable, pourquoi tu passes Q_temp en paramètre ? De ce que
   je comprends c'est le Qn divisé à chaque fois par ses diviseurs premiers (et
   leurs puissances). Cette variable n'a donc pas lieu d'exister en dehors du
   scope de la fonction. Si elle fait bien ça, vire la des paramètres. Je
   propose aussi de la renommé Qn_divided.
4. La documentation de init_hist_vects n'est pas précise : 'computes' ça veut
   dire quoi, en quoi nb_AQp sera nécessaire ?
## Fonction init_exp_vect
1. La documentation de init_exp_vect est pas archi claire je trouve. Tu
   pourrais mettre plus de détails sur ce qu'elle fait et comment elle le fait.
   Autre problème dans la signature, tu dis 'voir step_B.h' mais quand tu vas
   voir step_B.h tout ce que ça te dit c'est que ça sert à faire ce que tu
   veux. Au final t'es pas plus avancé sur pourquoi ça existe, pourquoi c'est
   utile.
2. Dans init_exp_vect le nom de variable D n'est pas terrible, change le non ?
3. Pour le premier if de init_exp_vect tu as des commentaires avant et un
   commentaire en inline dans le if. Regroupe les et keep it short.
4. Dans init_exp_vect, si je comprends bien tu ajoutes des valeurs à
   D.reduced_fb_indexes, au sens où tu mets des valeurs à la case
   nb_reduced_fb_indexes. Il faudrait préciser quelque part qu'au départ la
   taille du tableau reduced_fb_indexes est énorme. Là on peut croire que c'est
   un truc genre liste chaînée.
## Fonction create_AQ_pairs
1. Dans la signature, expliquer pourquoi on peut des fois trouver un facteur
   non trivial de Qn.
2. Tu peux pas mettre la partie expension de la fonction dans une fonction
   dédiée, quitte à lui redonner plein d'arguments ? Je parle des quatre
   premiers blocs de la boucle while. Pour alléger la notation tu peux créer
   une struct (avec typedef of course) pour toutes ces variables.
3. Rajoute un bref paragraphe dans la signature sur ce que fait la boucle.

# step_C.h
1. Renommer find_A_Q en calculate_A_Q ?

# step_C.c
1. Rajouter en commentaire que la raison pour laquelle on commence la
   numérotation de msb_indexes à 1 est technique (retour de mpz_sizeinbase).

# Misc
- Pourquoi ne pas remplacer tous les "subscript" par "indexe" ?
