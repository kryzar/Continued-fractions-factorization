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
