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
from timeit import default_timer as timer

if __name__ == '__main__':
    #%% Cas 1 - Versteeg 4.2
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

    #%% Cas 2 - Oberkampf 6.4.1
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
    # def MMS(x, y):
    #     T = T0 + Tx*np.cos(ax*np.pi*x) + Ty*np.sin(ay*np.pi*y) + Txy*np.sin(axy*np.pi*x*y)
    #     return T
    
    # Récupérer les MMS et les dériver
    x,y = sp.symbols('x y')
    T_MMS=T0 + Tx*sp.cos((ax*np.pi*x)/L)+Ty*sp.sin((ay*np.pi*y)/H)+Txy*sp.sin((axy*np.pi*x*y)/(L*H))
    f_T_MMS = sp.lambdify([x,y], T_MMS, "numpy")
    source = sp.diff(sp.diff(T_MMS,x),x)+sp.diff(sp.diff(T_MMS,y),y)
    f_source = sp.lambdify([x,y], source, "numpy")

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
    def q(x,y):
        source=f_source(x,y)
        return source
    
    def MMS(x,y):
        T=f_T_MMS(x,y)
        return T
    def Solution_Exacte(Centres):
        
        Sol_Ex=[]
        for i in range(len(Centres)):
            Sol_Ex.append(MMS(Centres[i,0],Centres[i,1]))
        
        Sol_Ex=np.array(Sol_Ex)
        return Sol_Ex
    

    for facteur in [1,2,4]:  # Ajouter des facteurs pour modifier le niveau de rafinement
        # Création du maillage pour la conception du solver
        mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': facteur*20,
                           'Ny': facteur*20
                           }
        bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
        plotter = MeshPlotter()

        # Initialisation du cas
        cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
        cas2.compute_mesh_and_connectivity()
        cas2.set_bc(bcdata)

        solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
        solver.solve("SPARSE")
        solution, area, analytical = cas2.get_solutions()

        post_traitement.set_data(cas2)

    post_traitement.genere_graphiques()

    # Affichage de champ scalaire avec pyvista du dernier maillage
    nodes, elements = plotter.prepare_data_for_pyvista(cas2.get_mesh())
    pv_mesh = pv.PolyData(nodes, elements)
    pv_mesh['Température'] = solution

    pl = pvQt.BackgroundPlotter()
    pl.add_mesh(pv_mesh, show_edges=True, scalars='Température', cmap="RdBu")
    pl.show()
    plt.plot(1, 1)  # Petit hack pour Audrey, sinon la fenetre de PyVista se ferme x)
    plt.show()
    
    #%% Solution Analytique et Numérique
    #Solution Avec le terme Cross Diffusion
    mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': 100,
                           'Ny': 100
                           }
    bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
    plotter = MeshPlotter()
    
    
    # Solution Numérique
    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    Centres=solver.centroids
    solver.solve("SPARSE")
    solutionNUM, area, analytical = cas2.get_solutions()
    
    #Solution Analytique Avec MMS
    solutionEX=np.zeros(len(Centres))
    for i in range(len(Centres)):
        solutionEX[i]=MMS(Centres[i,0],Centres[i,1])
        
    Figure1, (NUM,EX) = plt.subplots(1,2, figsize=(16,8), dpi=300)
        
    Figure1.suptitle("Comparaison Solution Exacte avec la solution Numérique")
    
    NUM.tricontourf(Centres[:,0],Centres[:,1],solutionNUM)
    NUM.set_xlabel("L(m)")
    NUM.set_ylabel("H(m)")
    NUM.set_title("Solution Numérique" )
    
    EX.tricontourf(Centres[:,0],Centres[:,1],solutionEX)
    EX.set_xlabel("L(m)")
    EX.set_ylabel("H(m)")
    EX.set_title("Solution Analytique MMS" )
    
    Figure1.tight_layout()
    
    # Enregistrer
    plt.savefig("Comp_Solutions_NUM_EX.png")
    #%% Solution Analytique et Numérique sur des coupes
    
    #Chercher l'indice des éléments à un X donné
    def Coupe_X(Coordonnees,X,Solution):
        Elements_ds_coupe=[]
        Solution_coupe=[]
        eps=10e-6 #Précision
        for i in range(len(Coordonnees)):
            if(np.abs(Coordonnees[i,0]-X)<eps):
                Elements_ds_coupe.append(Coordonnees[i,:])
                Solution_coupe.append(Solution[i])
        Elements_ds_coupe=np.array(Elements_ds_coupe)
        Solution_coupe=np.array(Solution_coupe)
        return Elements_ds_coupe,Solution_coupe
    
    def Coupe_Y(Coordonnees,Y,Solution):
        Elements_ds_coupe=[]
        Solution_coupe=[]
        eps=10e-6 #Précision
        for i in range(len(Coordonnees)):
            if(np.abs(Coordonnees[i,1]-Y)<eps):
                Elements_ds_coupe.append(Coordonnees[i,:])
                Solution_coupe.append(Solution[i])
        Elements_ds_coupe=np.array(Elements_ds_coupe)
        Solution_coupe=np.array(Solution_coupe)
        return Elements_ds_coupe,Solution_coupe
    
    #Test sur un Maillage QUAD
    #Solution Avec le terme Cross Diffusion
    mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': 10,
                           'Ny': 10
                           }
    bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
    plotter = MeshPlotter()
    
    # Solution Numérique
    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    Centres=solver.centroids
    solver.solve("SPARSE")
    solutionNUM, area, analytical = cas2.get_solutions()
    
    #Coupe sur X
    X_Coupe=4.25
    Elem_ds_coupeX,Solution_coupeX=Coupe_X(Centres,X_Coupe,solutionNUM)
    SolutionEX_coupeX=Solution_Exacte(Elem_ds_coupeX)
    
    #Coupe sur Y
    Y_Coupe=2.25
    Elem_ds_coupeY,Solution_coupeY=Coupe_Y(Centres,Y_Coupe,solutionNUM)
    SolutionEX_coupeY=Solution_Exacte(Elem_ds_coupeY)
    
    Figure1, (COUPEX,COUPEY) = plt.subplots(1,2, figsize=(16,8), dpi=300)
        
    Figure1.suptitle("Effet du maillage sur la solution")
    
    COUPEX.plot(Solution_coupeX,Elem_ds_coupeX[:,1],label="Solution Numérique")
    COUPEX.plot(SolutionEX_coupeX,Elem_ds_coupeX[:,1],label="Solution MMS")
    COUPEX.set_xlabel("Température")
    COUPEX.set_ylabel("Y(m)")
    COUPEX.set_title("Solution dans une coupe à X donnée" )
    
    COUPEY.plot(Elem_ds_coupeY[:,1],Solution_coupeY,label="Solution Numérique")
    COUPEY.plot(Elem_ds_coupeY[:,1],SolutionEX_coupeY,label="Solution MMS")
    COUPEY.set_xlabel("X(m)")
    COUPEY.set_ylabel("Température")
    COUPEY.set_title("Solution dans une coupe à Y donnée" )
    
    Figure1.tight_layout()
    
    # Enregistrer
    plt.savefig("Comp_Coupes.png")
    
    
    #%% Comparaison Pour deux maillages différents
    #Solution Avec le terme Cross Diffusion
    mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': 50,
                           'Ny': 50
                           }
    bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
    plotter = MeshPlotter()

    # Maillage QUad Structuré
    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    Centres=solver.centroids
    solver.solve("SPARSE")
    solutionST, area, analytical = cas2.get_solutions()
    
    #Maillage non Structuré
    mesh_parameters = {'mesh_type': 'QUAD',
                    'lc': 0.2
                    }
    bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
    plotter = MeshPlotter()

    # Maillage QUad Structuré
    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    CentresNS=solver.centroids
    solver.solve("SPARSE")
    solutionNS, area, analytical = cas2.get_solutions()
    
    Figure1, (NS,ST) = plt.subplots(1,2, figsize=(16,8), dpi=300)
        
    Figure1.suptitle("Effet du maillage sur la solution")
    
    NS.tricontourf(CentresNS[:,0],CentresNS[:,1],solutionNS)
    NS.set_xlabel("L(m)")
    NS.set_ylabel("H(m)")
    NS.set_title("Maillage non structuré" )
    
    ST.tricontourf(Centres[:,0],Centres[:,1],solutionST)
    ST.set_xlabel("L(m)")
    ST.set_ylabel("H(m)")
    ST.set_title("Maillage structuré" )
    
    Figure1.tight_layout()
    
    # Enregistrer
    plt.savefig("Comp_Maillages.png")
    
    #%% Comparaison Solution sans le cross diffusion et avec
    
    #Solution Avec le terme Cross Diffusion
    mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': 50,
                           'Ny': 50
                           }
    bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
    mesher = MeshGenerator()
    mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
    plotter = MeshPlotter()

    # Initialisation du cas
    cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
    cas2.compute_mesh_and_connectivity()
    cas2.set_bc(bcdata)

    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
    Centres=solver.centroids
    solver.solve("SPARSE")
    solutionAvec, area, analytical = cas2.get_solutions()
    
    solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=False)
    solver.solve("SPARSE")
    solutionSans, area, analytical = cas2.get_solutions()
     
    Figure1, (Sans,Avec) = plt.subplots(1,2, figsize=(16,8), dpi=300)
        
    Figure1.suptitle("Effet du terme Cross-Diffusion")
    
    Avec.tricontourf(Centres[:,0],Centres[:,1],solutionAvec)
    Avec.set_xlabel("L(m)")
    Avec.set_ylabel("H(m)")
    Avec.set_title("Terme Cross Diffusion Inclus" )
    
    Sans.tricontourf(Centres[:,0],Centres[:,1],solutionSans)
    Sans.set_xlabel("L(m)")
    Sans.set_ylabel("H(m)")
    Sans.set_title("Terme Cross Diffusion n'est pas Inclus" )
    
    Figure1.tight_layout()
    
    # Enregistrer
    plt.savefig("Comp_CD.png")
    
    
    #%% Mesure de temps de calcul
    #Définition du cas sur lequel on va tester
    # Création du maillage pour la conception du solver
    taille=[10,50,100,200]
    Time_Sparse=np.zeros(len(taille))
    Time_Dense=np.zeros(len(taille))
    for i in range(len(taille)):
        mesh_parameters = {'mesh_type': 'QUAD',
                           'Nx': taille[i],
                           'Ny': taille[i]
                           }
        bcdata = (['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS], ['DIRICHLET', MMS])
        mesher = MeshGenerator()
        mesh_obj = mesher.rectangle([0.0, L, 0.0, H], mesh_parameters)
        plotter = MeshPlotter()
    
        # Initialisation du cas
        cas2 = Case(mesh_obj, gamma=1, source_term=q, analytical_function=MMS)
        cas2.compute_mesh_and_connectivity()
        cas2.set_bc(bcdata)
    
        solver = MethodeVolumesFinisDiffusion(cas2, cross_diffusion=True)
        
        #Résolution Sparse
        start = timer()
        solver.solve("SPARSE")
        end = timer()
        Time_Sparse[i]=end-start
        
        #Résolution Dense
        start = timer()
        solver.solve("DENSE")
        end = timer()
        Time_Dense[i]=end-start
        
        #Affiche résultat
        Figure1, (Comp) = plt.subplots(1,1, figsize=(8,8), dpi=300)
            
        Figure1.suptitle("Comparaison Temps de calculs")

        Comp.plot(taille,Time_Sparse,label="Méthode Sparse")
        Comp.plot(taille,Time_Dense,label="Méthode Dense")
        Comp.set_xlabel("Taille de la matrice (Nombre d'éléments)")
        Comp.set_ylabel("Temps de résolution (s)")
        Comp.grid()
        Comp.legend()
        
        Figure1.tight_layout()
        
        # Enregistrer
        plt.savefig("Comparaison_temps_CPU.png")
    