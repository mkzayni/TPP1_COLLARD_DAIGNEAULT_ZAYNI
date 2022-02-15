"""
Date :    8 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Mohamad Karim Zayni (2167132)
Utilité : TPP1 - Méthode des volumes finis avec diffusion

"""

from meshConnectivity import MeshConnectivity

# Cas étudié regroupant le maillage, les conditions frontières, l'utilisation d'un solver
# et la solution calculée.
class Case:
    def __init__(self, mesh_obj, gamma, source_term, analytical_function=0):
        self.mesh_obj = mesh_obj                        # Maillage de la géométrie
        self.gamma = gamma                              # Coefficient diffusif
        self.source_term = source_term                  # Terme source
        self.analytical_function = analytical_function  # Fonction analytique
        self.time = {}

    # Exécute la connectivité avec le maillage généré.
    def compute_mesh_and_connectivity(self):
        self.conec = MeshConnectivity(self.mesh_obj, verbose=False)
        self.conec.compute_connectivity()

    # Applique les conditions frontières
    def set_bc(self, bcdata):
        self.bcdata = bcdata

    # Enregistre la solution numérique et analytique, l'aire et la position du centroides des éléments
    def set_solution(self, solution, analytical_solution, area, centroid):
        self.solution = solution
        self.analytical_solution = analytical_solution
        self.area = area
        self.centroid = centroid

    # Enregistre les temps de résolution (Dense et/ou Sparse)
    def set_resolution_time(self, matrix_type, time):
        self.time[matrix_type] = time

    # Retourne le maillage du cas étudié
    def get_mesh(self):
        return self.mesh_obj

    # Retourne les conditions frontières
    def get_bc(self):
        return self.bcdata

    # Retourne le coefficient diffusif
    def get_gamma(self):
        return self.gamma

    # Retourne la fonction de la solution analytique
    def get_analytical_function(self):
        return self.analytical_function

    # Retourne la solution numérique et analytique
    def get_solutions(self):
        return self.solution, self.analytical_solution

    # Retourne l'aire et la position des centroides des éléments
    def get_areas_and_centroids(self):
        return self.area, self.centroid

    # Retourne le/les temps de résolution (Dense et/ou Sparse)
    def get_times(self):
        return self.time
