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
2. Mettre les tests dans des fonctions et fichiers dédiés.
3. Commenter la fonction main.
