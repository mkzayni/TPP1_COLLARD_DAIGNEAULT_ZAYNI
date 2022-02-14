"""
Date :    8 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Mohamad Karim Zayni (2167132)
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
import sympy as sp

if __name__ == '__main__':
    #-----------------------------------------  Cas 1 - Versteeg 4.2 ------------------------------------------------#
    cas1_nom = 'Cas 1 - Versteeg 4.2'
    post_traitement1 = PostTraitement(cas1_nom)

    # Données du problème
    L = 0.02     # m
    k = 0.5  # W/mK
    def q(x, y): return 1000000  # W/m³
    def TA(x, y): return 100     # °C
    def TB(x, y): return 200     # °C
    def dphidn(x, y): return 0   # W/m²

    def solution_function(x, y):
        T = ((TB(x,y) - TA(x,y))/L + q(x,y)/(2*k)*(L-x))*x + TA(x, y)
        return T

    for facteur in [5, 10, 15, 20]:  # Niveau de rafinement
        # Création du maillage pour la conception du solver
        mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': facteur*5,
                           'Ny': facteur*5
                           }
        bcdata = (['DIRICHLET', TA], ['NEUMANN', dphidn], ['DIRICHLET', TB], ['NEUMANN', dphidn])
        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, 0.5*L], mesh_parameters)

        cas1 = Case(mesh_obj, gamma=k, source_term=q, analytical_function=solution_function)
        cas1.compute_mesh_and_connectivity()
        cas1.set_bc(bcdata)

        solver = MethodeVolumesFinisDiffusion(cas1, cross_diffusion=True)
        solver.solve(matrix_type="DENSE")
        solver.solve(matrix_type="SPARSE")
        solution, analytical = cas1.get_solutions()

        post_traitement1.set_data(cas1)

    # %% Ordre de convergence
    post_traitement1.show_error()

    # %% Temps de résolution avec matrice sparse et dense
    post_traitement1.show_time(title=f"Comparaison Temps de calculs pour  pour {cas1_nom}",
                                save_path="Comparaison_temps_CPU.png")

    # %% Solution Analytique et Numérique à Nx = Ny = 100
    post_traitement1.show_solutions(mesh=-1,
                                    title=f"Comparaison Solution Exacte avec la solution numérique pour {cas1_nom}",
                                    save_path="Comp_Solutions_NUM_EX_1.png")

    # %% Solution Analytique et Numérique sur des coupes à Nx = Ny = 10 (maillage 0)
    post_traitement1.show_plan_solutions(mesh=2,
                                         title=f"Effet du maillage sur la solution pour {cas1_nom}",
                                         save_path="Comp_Coupes_1.png",
                                         X_Coupe=L/2, Y_Coupe=0.5*L/2)

    # %% Comparaison Pour deux maillages différents
    # Maillage non Structuré
    mesh_parameters = {'mesh_type': 'TRI',
                           'Nx': 50,
                           'Ny': 50
                       }
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, 0.5 * L], mesh_parameters)
    cas1 = Case(mesh_obj, gamma=k, source_term=q, analytical_function=solution_function)
    cas1.compute_mesh_and_connectivity()
    cas1.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas1, cross_diffusion=True)
    solver.solve(matrix_type="SPARSE")

    post_traitement1.set_data(cas1)
    post_traitement1.show_mesh_differences(mesh1=[-1, "Maillage non structuré"],
                                           mesh2=[1,  "Maillage structuré"],
                                           title=f"Effet du maillage sur la solution pour {cas1_nom}",
                                           save_path="Comp_Maillages_1.png")

    # %% Comparaison Solution sans le cross diffusion et avec
    mesh_parameters = {'mesh_type': 'QUAD',
                       'lc': L/6
                       }
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, 0.5 * L], mesh_parameters)

    cas1 = Case(mesh_obj, gamma=k, source_term=q, analytical_function=solution_function)
    cas1.compute_mesh_and_connectivity()
    cas1.set_bc(bcdata)

    # Avec cross-diffusion
    solver = MethodeVolumesFinisDiffusion(cas1, cross_diffusion=True)
    solver.solve(matrix_type="SPARSE")
    post_traitement1.set_data(cas1)

    # Sans cross-diffusion
    solver = MethodeVolumesFinisDiffusion(cas1, cross_diffusion=False)
    solver.solve(matrix_type="SPARSE")
    post_traitement1.set_data(cas1)

    post_traitement1.show_mesh_differences(mesh1=[-2, "Effet du terme Cross-Diffusion"],
                                           mesh2=[-1, "Terme Cross Diffusion n'est pas Inclus"],
                                           title=f"Effet du terme Cross-Diffusion pour {cas1_nom}",
                                           save_path="Comp_CD_1.png")


    # --------------------------------------  Cas 2 - Oberkampf 6.4.1 ------------------------------------------------#
    cas2_nom = "Cas 2 - Oberkampf 6.4.1"
    post_traitement2 = PostTraitement(cas2_nom)

    # Données du problème
    L = 5  # m
    H = 3  # m
    T0 = 400  # K
    Tx = 45  # K
    Ty = 35  # K
    Txy = 27.5  # K
    ax = 1/3.
    ay = 1/4.
    axy = 1/2.

    # Récupérer les MMS et les dériver
    x,y = sp.symbols('x y')
    T_MMS=T0 + Tx*sp.cos(ax*np.pi*x)+Ty*sp.sin(ay*np.pi*y)+Txy*sp.sin(axy*np.pi*x*y)
    f_T_MMS = sp.lambdify([x, y], T_MMS, "numpy")
    source = sp.diff(T_MMS, x, 2)+sp.diff(T_MMS, y, 2)
    f_source = sp.lambdify([x, y], source, "numpy")

    # Le terme source ne fonctionne pas !
    def q(x,y):
        source = -f_source(x,y)
        return source
    
    def MMS(x,y):
        T = f_T_MMS(x,y)
        return T


    for facteur in [5, 10, 15, 20]:  # Ajouter des facteurs pour modifier le niveau de rafinement
        # Création du maillage pour la conception du solver
        mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': facteur*5,
                           'Ny': facteur*5
                           }
        bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)

        cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
        cas2.compute_mesh_and_connectivity()
        cas2.set_bc(bcdata)

        solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
        solver.solve(matrix_type="DENSE")
        solver.solve(matrix_type="SPARSE")
        solution, analytical = cas2.get_solutions()

        post_traitement2.set_data(cas2)

    # %% Ordre de convergence
    post_traitement2.show_error()

    # %% Temps de résolution avec matrice sparse et dense
    post_traitement2.show_time(title=f"Comparaison Temps de calculs pour  pour {cas2_nom}",
                               save_path="Comparaison_temps_CPU.png")

    # %% Solution Analytique et Numérique à Nx = Ny = 100
    post_traitement2.show_solutions(mesh=-1,
                                    title=f"Comparaison Solution Exacte avec la solution numérique pour {cas2_nom}",
                                    save_path="Comp_Solutions_NUM_EX_2.png")

    # %% Solution Analytique et Numérique sur des coupes à Nx = Ny = 10 (maillage 0)
    post_traitement2.show_plan_solutions(mesh=2,
                                         title=f"Effet du maillage sur la solution pour {cas2_nom}",
                                         save_path="Comp_Coupes_2.png",
                                         X_Coupe=2.5, Y_Coupe=1.5)

    # %% Comparaison Pour deux maillages différents
    # Maillage non Structuré
    mesh_parameters = {'mesh_type': 'QUAD',
                       'lc': 0.085
                       }
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    solver.solve(matrix_type="SPARSE")

    post_traitement2.set_data(cas2)
    post_traitement2.show_mesh_differences(mesh1=[-1, "Maillage non structuré"],
                                           mesh2=[1,  "Maillage structuré"],
                                           title=f"Effet du maillage sur la solution pour {cas2_nom}",
                                           save_path="Comp_Maillages_2.png")

    # %% Comparaison Solution sans le cross diffusion et avec
    mesh_parameters = {'mesh_type': 'QUAD',
                       'lc': 0.1
                       }
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)

    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    # Avec cross-diffusion
    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    solver.solve(matrix_type="SPARSE")
    post_traitement2.set_data(cas2)

    # Sans cross-diffusion
    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=False)
    solver.solve(matrix_type="SPARSE")
    post_traitement2.set_data(cas2)

    post_traitement2.show_mesh_differences(mesh1=[-2, "Effet du terme Cross-Diffusion"],
                                           mesh2=[-1, "Terme Cross Diffusion n'est pas Inclus"],
                                           title=f"Effet du terme Cross-Diffusion pour {cas2_nom}",
                                           save_path="Comp_CD_2.png")


    """# Affichage de champ scalaire avec pyvista du dernier maillage
    plotter = MeshPlotter()
    nodes, elements = plotter.prepare_data_for_pyvista(cas2.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Température'] = cas2.get_solutions()[0]

    pl = pvQt.BackgroundPlotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Température', cmap="RdBu")
    pl.show()"""
    plt.show()
