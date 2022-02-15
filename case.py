"""
Date :    8 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Mohamad Karim Zayni (2167132)
Utilité : TPP1 - Méthode des volumes finis avec diffusion

"""

from meshConnectivity import MeshConnectivity

# Cas étudié regroupant le maillage, les conditions frontières, l'utilisation d'un solver
# et la solution calculée.
class Case:
    """
    Preparer les données d'entrée pour l'exemple à traiter.

    Parmeters
    ---------
    mesh_obj: mesh
        Maillage du problème
    
    gamma: float
        Coefficient de diffusion
    
    source_term: function
        Fonction qui calcule le terme source sur le maillage
        
    analytical_function: function
        Fonction qui calcule la solution analytique sur le maillage


    Attributes
    ----------
    time: Dictionnaire de Str/float
        time de calcul pour les méthodes de résolution.

    """
    
    def __init__(self, mesh_obj, gamma, source_term, analytical_function=0):
        self.mesh_obj = mesh_obj                        # Maillage de la géométrie
        self.gamma = gamma                              # Coefficient diffusif
        self.source_term = source_term                  # Terme source
        self.analytical_function = analytical_function  # Fonction analytique
        self.time = {}

    # Exécute la connectivité avec le maillage généré.
    def compute_mesh_and_connectivity(self):
        """
        Exécute la connectivité avec le maillage généré.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.conec = MeshConnectivity(self.mesh_obj, verbose=False)
        self.conec.compute_connectivity()

    # Applique les conditions frontières
    def set_bc(self, bcdata):
        """
        Modifier les conditions aux frontières.

        Parameters
        ----------
        bcdata: liste de float et str.
            Type + valeur des conditions aux limites imposées. 

        Returns
        -------
        None

        """
        self.bcdata = bcdata

    # Enregistre la solution numérique et analytique, l'aire et la position du centroides des éléments
    def set_solution(self, solution, analytical_solution, area, centroid):
        """
        Enregistre la solution numérique et analytique, l'aire et la position du centroides des éléments.

        Parameters
        ----------
        solution: ndarray float.
            Solution numérique. 
        
        analytical_solution: ndarray float.
            Solution Exacte. 
            
        area: ndarray float
            Surface des éléments
            
        centroid: ndarray float
            positions des centres des éléments
        
        Returns
        -------
        None

        """
        self.solution = solution
        self.analytical_solution = analytical_solution
        self.area = area
        self.centroid = centroid

    # Enregistre les temps de résolution (Dense et/ou Sparse)
    def set_resolution_time(self, matrix_type, time):
        """
        Enregistre les temps de résolution (Dense et/ou Sparse)

        Parameters
        ----------
        matrix_type: str.
            type de matrice Dense ou Sparse. 
        
        time: float.
            Temps de résolution. 
        
        Returns
        -------
        None

        """
        self.time[matrix_type] = time

    # Retourne le maillage du cas étudié
    def get_mesh(self):
        """
        Retourne le maillage du cas étudié

        Parameters
        ----------
        None
        
        Returns
        -------
        mesh: Maillage étudié

        """
        return self.mesh_obj

    # Retourne les conditions frontières
    def get_bc(self):
        """
        Retourne les conditions frontières

        Parameters
        ----------
        None
        
        Returns
        -------
        bcdata: liste de float et str.
            Type + valeur des conditions aux limites imposées. 

        """
        return self.bcdata

    # Retourne le coefficient diffusif
    def get_gamma(self):
        """
        Retourne le coefficient diffusif

        Parameters
        ----------
        None
        
        Returns
        -------
        gamma: float.
            Coefficient de diffusion. 
        """
        return self.gamma

    # Retourne la fonction de la solution analytique
    def get_analytical_function(self):
        """
        Retourne la fonction de la solution analytique

        Parameters
        ----------
        None
        
        Returns
        -------
        analytical_solution: function.
            Fonction de calcul de la solution analytique 
        """
        return self.analytical_function

    # Retourne la solution numérique et analytique
    def get_solutions(self):
        """
        Retourne la solution numérique et analytique

        Parameters
        ----------
        None
        
        Returns
        -------
        solution, analytical_solution: ndarray float.
            solution numérique et analytique
        """
        return self.solution, self.analytical_solution

    # Retourne l'aire et la position des centroides des éléments
    def get_areas_and_centroids(self):
        """
        Retourne l'aire et la position des centroides des éléments

        Parameters
        ----------
        None
        
        Returns
        -------
        area,centroid: ndarray float.
            surface des éléments et la position de leur centre
        """
        return self.area, self.centroid

    # Retourne le/les temps de résolution (Dense et/ou Sparse)
    def get_times(self):
        """
        Retourne le/les temps de résolution (Dense et/ou Sparse)

        Parameters
        ----------
        None
        
        Returns
        -------
        time: dictionnaire str/float.
            Temps de calcul correspondant à la méthode utilisé
        """
        return self.time
