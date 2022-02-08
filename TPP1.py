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
    post_traitement = PostTraitement('Q2')

    for facteur in [5]:
        # Création du maillage pour la conception du solver (rectangle simple)
        mesh_parameters = {'mesh_type': 'TRI',
                           'Nx': 2,
                           'Ny': 2
                           }
        bcdata = (['DIRICHLET', 1], ['DIRICHLET', 0], ['DIRICHLET', 5], ['DIRICHLET', 0],)

        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, 2.0, 0.0, 1.0], mesh_parameters)
        plotter = MeshPlotter()
        cas = Case(mesh_obj)
        cas.compute_mesh_and_connectivity()
        cas.set_bc(bcdata)
        solver = MethodeVolumesFinisDiffusion(cas)



        """
        # Création du maillage (Cercle dans un rectangle non structuré)
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


        # Question 2 : Reconstruction du gradient par moindres carrées
        phi = lambda x, y: np.sin(x) + np.cos(y)
        dphi_analytical = [lambda x: np.cos(x), lambda y: -np.sin(y)]

        # Conditions de Dirichlet ou Neumann
        bcdata = (['DIRICHLET'], ['DIRICHLET'], ['DIRICHLET'], ['DIRICHLET'], ['NEUMANN'])
        cas.set_bc(bcdata)
        solver = MethodeVolumesFinisDiffusion(cas, phi, dphi_analytical)
        solver.solve()
        solution, analytical, area = cas.get_solutions()

        post_traitement.set_data(cas)

    post_traitement.genere_graphiques()"""


    # Affichage de champ scalaire avec pyvista du dernier maillage
    plotter.plot_mesh(mesh_obj, label_points=True, label_elements=True, label_faces=True)

    """plotter = MeshPlotter()
    nodes, elements = plotter.prepare_data_for_pyvista(cas.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Gradient'] = solution

    pl = pvQt.BackgroundPlotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Gradient', cmap="RdBu")
    pl.show()"""
    plt.show()