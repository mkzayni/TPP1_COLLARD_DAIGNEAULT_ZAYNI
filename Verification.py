"""


MEC6616 Aérodynamique Numérique

@author: Mohamad Karim ZAYNI
Matricule: 2167132 

@author: Mohamed Dhia SLAMA
Matricule:2115178
    
    
    

"""


#----------------------------------------------------------------------------#
#                                 MEC6616                                    #
#                               LABORATOIRE 3                                #
#               SLAMA Mohamed Dhia, ZAYNI Mohamad Karim                      #
#----------------------------------------------------------------------------#

#%% NOTES D'UTILISATION

"""
Cette classe gère la vérification de notre programme

"""

#%% Importation des librairies

import sympy as sp
import numpy as np
import pyvista as pv
import pyvistaqt as pvQt
from meshGenerator import MeshGenerator
from meshConnectivity import MeshConnectivity
from meshPlotter import MeshPlotter
import Validate_Mesh as VM
import matplotlib.pyplot as plt
import sys


#%% Définition de la classe

class Verification():
    
    # Initilisation de la classe
    def __init__(self, _ListePas,_Prob):
        self._ListePas = _ListePas #Liste des pas pour la vérification list[int]
        self._Prob=_Prob #Problème ou géométrie Traitée [Validate_Mesh]
        
        #Vérifier qu'on a minimum 3 éléments dans notre liste de pas
        if len(_ListePas) < 3:
            print("(ERREUR CRITIQUE) : Les étapes de vérifications nécessitent au moins 3 pas de temps.")
            print("(SOLUTION 1) : Révisez")
            print("(ARRÊT DU PROGRAMME")
            sys.exit()
            
    # Accesseurs des attributs
    def get_ListePas(self):
        return self._ListePas
    
    def get_Prob(self):
        return self._Prob
    

    # Modificateurs des attributs        
    def set_ListePas(self, _New):
        self._ListePas = _New
        
    def set_Prob(self,_New):
        self._Prob=_New
    
        
    # Définition des méthodes
    # Calcul de l'ordre standard de convergence : ln(e2h/eh)/ln(r)
    def OrdreStandard(self, _x, _y):
        _x = np.log(_x)
        _y = np.log(_y)
        # Construction de la matrice des X, pour alimenter la fonction existante de Numpy
        Matrice = np.vstack([_x, np.ones(len(_x))]).T
        # Construction de la matrice des Y, pour alimenter la fonction existante de Numpy
        Colonne = _y[:, np.newaxis]
        # Sortir la liste de deux coefficents, en première place dans les réponses, puis en sortir le premier coefficent
        # (pente de la droite de régression fittée par méthode des moindres carrés)
        Ordre = np.linalg.lstsq(Matrice, Colonne, rcond=None)[0][0]
        return Ordre

    # Normes d'erreur
    def Normes(self, _U1, _U2, _Methode):
        VolumeMaille = 1 # Car on parle en nombre de pas
        VolumeDomaine = np.amin([len(_U1), len(_U2)]) # Nombre total de pas du plus petit vecteurs
        # Redimensionner les vecteurs sur la taille du plus petit, ne rien faire si Ratio = 1
        Ratio = int((np.amax([len(_U1), len(_U2)])-1)/(np.amin([len(_U1), len(_U2)])-1))
        if len(_U2) > len(_U1):
            _U2 = _U2[::Ratio]
        elif len(_U2) < len(_U1):
            _U1 = _U1[::Ratio]
        # Calculer la norme sur les vecteurs redimensionnés
        if _Methode == "L1":
            L = np.sum(np.absolute(_U2 - _U1)*VolumeMaille)/VolumeDomaine
        elif _Methode == "L2":
            L = np.sqrt(np.sum(np.square(_U2 - _U1)*VolumeMaille)/VolumeDomaine)
        elif _Methode == "Li":
            L = np.amax(np.absolute(_U2 - _U1))
        else:
            self.ErrorMessage("Norme")
        return L
        
    #Calcul de l'erreur + Affichage
    def Executer(self, _Methode,phi,CL_Face,CL_type,CL_Val):
       
        #Récupération des données
        Pas=self.get_ListePas()
        Liste_dx=[]
        Geometrie=self.get_Prob()
        Erreurs=[]
        
        for i in range(len(Pas)):
            
            #Changement de maillage + Résolution
            Liste_dx.append(1/Pas[i])
            mesher = MeshGenerator()
            mesh_parameters = {'mesh_type': 'TRI',
                                'Nx': Pas[i],
                                'Ny': Pas[i]
                                }
            mesh_obj = mesher.rectangle([0.0, 1.0, 0.0, 1.0], mesh_parameters)

            conec = MeshConnectivity(mesh_obj)
            conec.compute_connectivity()
            
            Geometrie.set_Mesh(mesh_obj)
            
            phi,CL_Face,CL_type,CL_Val=Geometrie.Conditions()
            Solution_Num=Geometrie.Gradient_LS(phi,CL_Face,CL_type,CL_Val)
            
            #Solution Exacte
            Solution_Ex=Geometrie.Solution_exacte()
          
            #Calcul de l'Erreur selon la norme choisie
            Erreur=self.Normes(Solution_Num[:,0], Solution_Ex[:,0],_Methode)
            Erreurs.append(Erreur)
            
        #Graphe: Erreur + solution dernière maillage
        Figure1, (Conv) = plt.subplots(1,1, figsize=(8,8), dpi=300)
            
        Figure1.suptitle("VÉRIFICATION DE CODE")
        
        # Calcul de l'ordre
        p=self.OrdreStandard(Liste_dx, Erreurs)
        
        
        #Calcul 1/dx
        for i in range(len(Liste_dx)):
            Liste_dx[i]=1.0/Liste_dx[i]
            
        Conv.loglog(Liste_dx, Erreurs)
        Conv.scatter(Liste_dx, Erreurs)
        Conv.set_xlabel("ln(1/dx)")
        Conv.set_ylabel("ln( norme "+_Methode+" )")
        Conv.grid()
        
        Conv.set_title("Diagramme de convergence p̂ = %.2f" %p)
        Figure1.tight_layout()
        
        # Enregistrer
        plt.savefig("Verif_Code.png")

        