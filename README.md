# TPs-Optimisation

La première partie de ce TP-projet concerne les problèmes d’optimisation sans contraintes. On étudie la méthode de Newton et sa globalisation par l’algorithme des régions de confiance. La résolution du sous-problème des régions de confiance sera réalisée de deux façons, soit à l’aide du point de Cauchy, soit par l’algorithme du Gradient Conjugué Tronqué. La seconde partie du projet exploite la partie précédente pour résoudre des problèmes d’optimisation avec contraintes par l’algorithme du Lagrangien augmenté.

# Optimisation sans contraintes

Dans cette partie, on s’intéresse à la résolution du problème

\min _{x \in \mathbb{R}^{n}} f(x)min 
x∈R 
n
 
​	
 f(x)

où la fonction ff est de classe C^{2}C 
2
  sur R^{n}R 
n
  . On cherche donc à exploiter l’information fournie par ses dérivées première et seconde, que l’on représente en tout point x par le vecteur gradient \nabla f (x) \in R^{n}∇f(x)∈R 
n
  et la matrice Hessienne \nabla^{2} f (x) \in R^{n\times n}∇ 
2
 f(x)∈R 
n×n
 .}
