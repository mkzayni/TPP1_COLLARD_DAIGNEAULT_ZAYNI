# TPP1_COLLARD_DAIGNEAULT_ZAYNI


Fait : 
 - (Audrey) Définition d'un maillage simple permettant de vérifier les calculs géométriques/vectoriels
 - (Audrey) Fonction find_centroids() permettant de calculer les centres de face pour tous les éléments
 - (Audrey) Fonction permettant de calculer l'aire des éléments (TRI ou QUAD)
 - (Audrey) Solver jusqu'au terme du direct gradient
 - (Audrey) Implémentation du cas 1  (Versteeg 4.2) (sans Sdc so far)
 - (Audrey) Gérer la méthode itérative pour que ça fonctionne
 - (Karim) Implémenter les modules de post-processing (voir énoncé) - Comparaison Num-MMS/Comp sur des coupes ** A verifier si t'es d'accord avec la méthode **/Comparaison pour différents maillages
 - (Karim) Comparer les cas avec et sans terme de cross-diffusion 
 - (Karim) Calculer le computational time avec full matrix et sparse matrix (A verifier sur la taille de notre maillage, mon pc semble mourir à partir de 400)

En cours :
- (Audrey) À VÉRIFIER Calcul du cross-diffusion term (gradϕ•η) par la méthode des moindre carrés (voir LAP3)
- (Audrey) WIP Résoudre le cas 2 (Oberkampf & Roy 6.3.4) => Bug avec la MMS et/ou terme source


 -(Karim) Faire diagramme UML (à vérifier si j'ai eu le bon UML)

À faire : 
 - Terminer l'implémentation complète
 - Tester et vérifier que ça fonctionne pour tous les types de maillages (semble ok)
