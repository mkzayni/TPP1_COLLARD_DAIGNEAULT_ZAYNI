# TPP1_COLLARD_DAIGNEAULT_ZAYNI


Fait : 
 - (Audrey) Définition d'un maillage simple permettant de vérifier les calculs géométriques/vectoriels
 - (Audrey) Fonction find_centroids() permettant de calculer les centres de face pour tous les éléments
 - (Audrey) Fonction permettant de calculer l'aire des éléments (TRI ou QUAD)
 - (Audrey) Solver jusqu'au terme du direct gradient
 - (Audrey) Implémentation du cas 1  (Versteeg 4.2) (sans Sdc so far)


À faire : 
 - Calcul du cross-diffusion term (gradϕ•η) par la méthode des moindre carrés (voir LAP3)
 - Gérer la méthode itérative pour que ça fonctionne
 - Terminer l'implémentation complète
 - Tester et vérifier que ça fonctionne pour tous les types de maillages (semble ok)
 - Résoudre le cas 2 (Oberkampf & Roy 6.3.4)
 - Implémenter les modules de post-processing (voir énoncé)
 - Comparer les cas avec et sans terme de cross-diffusion
 - Calculer le computational time avec full matrix et sparse matrix
 - Faire diagramme UML