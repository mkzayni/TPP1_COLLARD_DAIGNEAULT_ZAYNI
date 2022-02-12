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
    """Post_traitement = PostTraitement('Versteeg 4.2')

    # Données du problème
    L = 0.02     # m
    gamma = 0.5  # W/mK
    def q(x, y): return 1000000  # W/m³
    def TA(x, y): return 100     # °C
    def TB(x, y): return 200     # °C
    def dphidn(x, y): return 0   # W/m²

    for facteur in [1]:  # Ajouter des facteurs pour modifier le niveau de rafinement
        # Création du maillage pour la conception du solver
        mesh_parameters = {'mesh_type': 'MIX',
                           'Nx': facteur*50,
                           'Ny': facteur*10
                           }
        bcdata = (['DIRICHLET', TA], ['NEUMANN', dphidn], ['DIRICHLET', TB], ['NEUMANN', dphidn])

        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, 0.5*L], mesh_parameters)
        plotter = MeshPlotter()

        # Initialisation du cas
        cas = Case(mesh_obj, gamma, source_term=q)
        cas.compute_mesh_and_connectivity()
        cas.set_bc(bcdata)

        solver = MethodeVolumesFinisDiffusion(cas, cross_diffusion=True)
        solver.solve()
        solution, area = cas.get_solutions()
 

        post_traitement.set_data(cas)

    post_traitement.genere_graphiques()


    # Affichage de champ scalaire avec pyvista du dernier maillage
    nodes, elements = plotter.prepare_data_for_pyvista(cas.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Température'] = solution

    pl = pvQt.BackgroundPlotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Température', cmap="RdBu")"""


    # --------------------------------------  Cas 2 - Oberkampf 6.4.1 ------------------------------------------------#
    post_traitement = PostTraitement('Oberkampf 6.4.1')

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

    # L'équation manufacturée à un /L... mais ça ne donne pas les bons résultats...
    def MMS(x, y):
        T = T0 + Tx*np.cos(ax*np.pi*x) + Ty*np.sin(ay*np.pi*y) + Txy*np.sin(axy*np.pi*x*y)
        return T

    """# Test de la MMS
    GRIDx = np.linspace(0,5,500)
    GRIDy = np.linspace(0,3,300)
    z = np.zeros([len(GRIDy),len(GRIDx)])

    for i in range(len(GRIDx)):
        for j in range(len(GRIDy)):
            z[j,i] = MMS(GRIDx[i], GRIDy[j])

    c = plt.pcolor(GRIDx, GRIDy, z)
    plt.colorbar(c)"""

    # Le terme source ne fonctionne pas !
    def q(x,y): return MMS(x,y)

    for facteur in [1]:  # Ajouter des facteurs pour modifier le niveau de rafinement
        # Création du maillage pour la conception du solver
        mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': 75,
                           'Ny': 75
                           }
        bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])

        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
        plotter = MeshPlotter()

        # Initialisation du cas
        cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_solution=MMS)
        cas2.compute_mesh_and_connectivity()
        cas2.set_bc(bcdata)

        solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
        solver.solve()
        solution, area = cas2.get_solutions()

        #post_traitement.set_data(cas2)

    #post_traitement.genere_graphiques()

    # Affichage de champ scalaire avec pyvista du dernier maillage
    nodes, elements = plotter.prepare_data_for_pyvista(cas2.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Température'] = solution

    pl = pvQt.BackgroundPlotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Température', cmap="RdBu")
    pl.show()
    plt.plot(1, 1)  # Petit hack pour Audrey, sinon la fenetre de PyVista se ferme x)
    plt.show()