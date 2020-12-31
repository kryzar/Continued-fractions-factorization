# Fractions continues
## Intuition
- c'est une expression de la forme truc
- fini ou infini ; si c'est fini alors c'est juste un nombre rationnel « normal
  ». Si c'est infini il y a du travail à faire pour définir proprement les
  notions d'égalité et tout
- à la base motivation pour approcher des irrationnels par des fractions, voir
  les approximations de π
- notez qu'il n'y a que des 1 sur les numérateurs, sinon on appelle ça frac.
  continue généralisée 
- partant d'un réel on lui associe une fraction continue comme ça
- Note: dans le cas d'un rationnel ce procéssus se déduit de l'algorithme
  d'Euclide
- ça s'arrête jamais si le réel de départ est irrationel
- exemple de π
## Formalisation
- formellement on peut définir une fraction continue (finie ou infinie) comme
  une suite qui vérifie ça et ça
- on la note comme ça
- on définit désormais les réduites, on se place dans Q et l'on se donne une
  famille X0, X1, … une fmaille d'indéterminées indexée par les entiers. On
  définit alors une famille de réduites associées à la famille
- définition réduites
- Note : l'utilisation des indéterminées permet de ne pas se soucier des
  divisions par zéro
- les réduites sont donc des éléments de Q(X_1, …, X_n-1)
- Exemple : premières réduites
- à toute fraction continue on (finie ou infinie) peut canoniquement associer
  une famille de réduites (séparer les cas finis et infinis)
- ces réduites sont dans tous les cas des rationnels. si le nombre de départ
  est rationnel, il est égal à sa dernière réduite. Si le nombre de départ est
  irationnel il n'est égal à aucune de ses réduites mais on n'a l'égalité
  fondamentale suivante (utile dans le programme)
- il convient maintenant de donner un sens à l'égalité entre une fraction
  continue. Si la fc est finie, c'est quand r = [a_1, …, a_n]. Si la fraction
  continue est infinie ou peut dire ça
- définition, on dit qu'une fraction continue infinie et un réel sont égaux si
  la suite des réduites converge verts le réel. Dans ce cas on note machin =
  truc
- on note comme ça
- Exemple : nombre d'or
- à quoi est égale une fraction continue ?
- si la fraction continue est finie, elle est évidemment égale au
  rationnel [a_1, …, a_n]. L'application qui à une fraction continue associe un
  rationnel n'est pas bijective car non injective, le rationnel [a_1, …, a_n]
  est égal à la fraction continue (autre).
- on a alors le résultat suivant : tout irationnel est égal qui lui est
  canoniquement associée. Mieux on a une bijection blabla
## Fraction continue d'un irrationel quadratique
- on s'intéresse ici aux fractions continues d'irrationnels quadratiques et
  montrons que les fractions continues d'irrationnels quadratiques sont
  sujettes à des phénomènes de périodicité
- déf. nombre réel algébrique sur Q de degré 2
- Note : forcément irrationel
- notation d'une période pour une fraction continue
- on a le théorème fondamental suivant
- Note : il est de coutûme de l'associer à un résultat sur les quadratiques
  irationnels réduits mais ce n'est pas utile ici
- Intéressons nous plus précisément au cas de \sqrt(N), qui sera utile pour
  notre algorithme — c'est bien un irrationel quadratique (non réduit)
- déjà, on peut renforcer le th. précédent avec le th. Legendre
- exemples avec sqrt(14)
- en suite, on peut lui associer une autre suite, en vertu du lemme suivant :
- lemme. [prop. 2.5.2] qui dit que truc est encore un irr. quad.
- en reprenant l'algo précédent et l'identité bidule on a zeta = [q0, …, qn-1,
  zetan], où zetan est encore un irr. quad.. Il s'écrit donc zetan = (sqrt(N) +
  Pn) / Qn. La suite des Qn est éminemment importante. Elle vérifie les
  propriétés suivantes qui sont utiles dans l'algo : 
  * Pn < truc, Qn < machin
  * identité 1) et donc congruence 2) du M-B
