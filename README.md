# TPP1_COLLARD_DAIGNEAULT_ZAYNI

Maillage de test : 

![image](maillage_test.png)


Maillage choisi : 

![image](maillage.png)

Fait : 
 - (Audrey) Définition d'un maillage simple permettant de vérifier les calculs géométriques/vectoriels
 - (Audrey) Fonction find_centroids() permettant de calculer les centres de face pour tous les éléments
 - (Audrey) Fonction permettant de calculer l'aire des éléements (TRI ou QUAD)
 - (Audrey) Fonction permettant de calculer la longueur des faces (pas nécessaire???)
 - (Audrey) Solver jusqu'au terme du direct gradient



À faire : 
 - Calcul du cross-diffusion term (gradϕ•η) par la méthode des moindre carrés (voir LAP3)
 - Gérer la méthode itérative pour que ça fonctionne
 - Terminer l'implémentation complête
 - Tester et vérifier que ça fonctionne pour tous les types de maillages
 - Résoudre le cas 1 (Versteeg 4.2)
 - Résoudre le cas 2 (Oberkampf & Roy 6.3.4)
 - Implémenter les modules de post-processing (voir énoncé)
 - Comparer les cas avec et sans terme de cross-diffusion
 - Calculer le computational time avec full matrix et sparse matrix
 - Faire diagramme UML