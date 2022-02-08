"""
Date :    2 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Zacharie Prigent (2096785)
Utilité : LAF3 : Tests de maillage et reconstruction du gradient par moindres carrées
"""
# IMPORT STATEMENTS
import numpy as np
import pyvista as pv
from meshGenerator import MeshGenerator
from meshPlotter import MeshPlotter
from solver import GradientMoindresCarres
from unitTests import UnitTest
from case import Case
from post_traitement import PostTraitement
import matplotlib.pyplot as plt

if __name__ == '__main__':
    post_traitement = PostTraitement('Q2')

    for facteur in [1, 2, 3, 4]:
        # Création du maillage (Cercle dans un rectangle non structuré)
        print(f'\nMaillage {facteur}')
        mesh_parameters = {'mesh_type': 'MIX',
                           'lc_rectangle': facteur*0.05,
                           'lc_circle': facteur*0.01
                           }
        rayon = 0.25
        mesher = MeshGenerator()
        mesh_obj = mesher.circle([0.0, 2.0, 0.0, 1.0], rayon, mesh_parameters)

        # Initialisation du cas étudié
        cas = Case(mesh_obj, nb_trou=1)
        cas.compute_mesh_and_connectivity()


        # Question 1 : Tests de base sur les maillages
        unittest = UnitTest(cas)
        unittest.test_relation_euler()
        unittest.test_divergence()


        # Question 2 : Reconstruction du gradient par moindres carrées
        phi = lambda x, y: np.sin(x) + np.cos(y)
        dphi_analytical = [lambda x: np.cos(x), lambda y: -np.sin(y)]

        # Conditions de Dirichlet ou Neumann
        bcdata = (['DIRICHLET'], ['DIRICHLET'], ['DIRICHLET'], ['DIRICHLET'], ['NEUMANN'])
        cas.set_bc(bcdata)
        solver = GradientMoindresCarres(cas, phi, dphi_analytical)
        solver.solve()
        solution, analytical, area = cas.get_solutions()

        post_traitement.set_data(cas)

    post_traitement.genere_graphiques()
    plt.show()

    # Affichage de champ scalaire avec pyvista du dernier maillage
    plotter = MeshPlotter()
    nodes, elements = plotter.prepare_data_for_pyvista(cas.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Gradient'] = solution

    pl = pv.Plotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Gradient', cmap="RdBu")
    pl.show()