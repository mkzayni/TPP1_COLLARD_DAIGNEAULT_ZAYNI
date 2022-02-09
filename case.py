from meshConnectivity import MeshConnectivity

# Cas étudié regroupant le maillage, les conditions frontières, l'utilisation d'un solver
# et la solution calculée.
class Case:
    def __init__(self, mesh_obj, gamma, source_term):
        self.mesh_obj = mesh_obj        # Maillage de la géométrie
        self.gamma = gamma              # Coefficient diffusif
        self.source_term = source_term  # Terme source

    # Exécute la connectivité avec le maillage généré.
    def compute_mesh_and_connectivity(self):
        self.conec = MeshConnectivity(self.mesh_obj, verbose=False)
        self.conec.compute_connectivity()

    # Applique les conditions frontières
    def set_bc(self, bcdata):
        self.bcdata = bcdata

    # Enregistre la solution numérique et analytique puis l'aire des cellules
    def set_solution(self, solution, area):
        self.solution = solution
        self.area = area

    # Permet d'obtenir le maillage du cas étudié
    def get_mesh(self):
        return self.mesh_obj

    # Permet d'obtenir les conditions frontières
    def get_bc(self):
        return self.bcdata

    # Permet d'obtenir le coefficient diffusif
    def get_gamma(self):
        return self.gamma

    # Retourne la solution numérique et analytique puis l'aire des cellules
    def get_solutions(self):
        return self.solution, self.area
