"""
Date :    8 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Mohamad Karim Zayni ()
Utilité : TPP1 - Méthode des volumes finis avec diffusion
"""
# IMPORT STATEMENTS
import numpy as np
import pyvista as pv
import pyvistaqt as pvQt
from meshGenerator import MeshGenerator
from meshPlotter import MeshPlotter
from solver import MethodeVolumesFinisDiffusion
from case import Case
from post_traitement import PostTraitement
import matplotlib.pyplot as plt

if __name__ == '__main__':
    #-----------------------------------------  Cas 1 - Versteeg 4.2 ------------------------------------------------#
    post_traitement = PostTraitement('Versteeg 4.2')

    # Données du problème
    L = 0.02     # m
    gamma = 0.5  # W/mK
    q = 1000000. # W/m³ (doit être à 1000000 pour le 4.2)
    TA = 100     # °C
    TB = 200     # °C

    for facteur in [1]:  # Ajouter des facteurs pour modifier le niveau de rafinement
        # Création du maillage pour la conception du solver
        mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': facteur*100,
                           'Ny': facteur*5
                           }
        bcdata = (['DIRICHLET', TA], ['NEUMANN', 0], ['DIRICHLET', TB], ['NEUMANN', 0])

        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, 0.5*L], mesh_parameters)
        plotter = MeshPlotter()

        # Initialisation du cas
        cas = Case(mesh_obj, gamma, source_term=q)
        cas.compute_mesh_and_connectivity()
        cas.set_bc(bcdata)

        cross_diffusion = False
        solver = MethodeVolumesFinisDiffusion(cas, cross_diffusion)
        solver.solve()
        solution, area = cas.get_solutions()

        """
        # Question 2 : Reconstruction du gradient par moindres carrées
        phi = lambda x, y: np.sin(x) + np.cos(y)
        dphi_analytical = [lambda x: np.cos(x), lambda y: -np.sin(y)]

        # Conditions de Dirichlet ou Neumann
        bcdata = (['DIRICHLET'], ['DIRICHLET'], ['DIRICHLET'], ['DIRICHLET'], ['NEUMANN'])
        cas.set_bc(bcdata)
        solver = MethodeVolumesFinisDiffusion(cas, phi, dphi_analytical)
        
        

        post_traitement.set_data(cas)

    post_traitement.genere_graphiques()"""


    # Affichage de champ scalaire avec pyvista du dernier maillage
    nodes, elements = plotter.prepare_data_for_pyvista(cas.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Température'] = solution

    pl = pvQt.BackgroundPlotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Température', cmap="RdBu")
    pl.show()
    plt.plot(1,1)  # Petit hack pour Audrey, sinon la fenetre de PyVista se ferme x)
    plt.show()