import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve


# Solveur utilisant la méthode des volumes finis avec diffusion seulement
class MethodeVolumesFinisDiffusion:
    def __init__(self, case, cross_diffusion):
        self.case = case                        # Cas à résoudre
        self.mesh_obj = case.get_mesh()         # Maillage du cas
        self.bcdata = case.get_bc()             # Conditions frontières
        self.cross_diffusion = cross_diffusion  # Si le terme de cross diffusion est activé

        self.centroids = np.zeros([self.mesh_obj.get_number_of_elements(), 2])  # Array des centroides
        self.volumes = np.zeros([self.mesh_obj.get_number_of_elements()])       # Array de l'aire des elements

        self.preprocessing()

    # Effectue les calculs relatifs au maillage préalablement à l'utilisation du solver
    def preprocessing(self):
        # Détermine les centroides de l'élément par la moyenne des coordonnées
        def find_centroid(i_elem):
            nodes = self.mesh_obj.get_element_to_nodes(i_elem)
            for i_node in nodes:
                x, y = self.mesh_obj.get_node_to_xycoord(i_node)[0], self.mesh_obj.get_node_to_xycoord(i_node)[1]
                self.centroids[i_elem][0] += x
                self.centroids[i_elem][1] += y

            self.centroids[i_elem] /= len(nodes)

        # Calcul l'aire de l'élément par une méthode utilisant les déterminants
        def compute_area(i_elem):
            nodes = self.mesh_obj.get_element_to_nodes(i_elem)
            area_matrices = [np.zeros([2, 2]) for i in range(len(nodes))]
            for i in range(len(nodes)):
                x, y = self.mesh_obj.get_node_to_xycoord(nodes[i])[0], self.mesh_obj.get_node_to_xycoord(nodes[i])[1]
                area_matrices[i][:][0] = [x, y]
                area_matrices[i - 1][:][1] = [x, y]

            self.volumes[i_elem] = np.sum([np.linalg.det(area_matrices[i]) for i in range(len(nodes))]) / 2.

        # Calculs pour les éléments (centres d'élément et aires)
        for i_elem in range(self.mesh_obj.get_number_of_elements()):
            find_centroid(i_elem)
            compute_area(i_elem)

    def solve(self):
        # Itinitialisation des matrices et du terme de cross-diffusion
        NELEM = self.mesh_obj.get_number_of_elements()
        A = np.zeros((NELEM, NELEM))
        B = np.zeros(NELEM)
        Sdc = 0  # Cross-diffusion term reste nul si False

        # Calcule les distances et vecteurs nécessaires selon les coordonnées fournies
        def compute_lengths_and_unit_vectors(pta, ptb, ptA, ptP):
            (xa, ya), (xb, yb), (xA, yA), (xP, yP) = pta, ptb, ptA, ptP

            # Détermination des distances
            dx, dy = (xb - xa), (yb - ya)
            dA = np.sqrt(dx ** 2 + dy ** 2)
            dKSI = np.sqrt((xA - xP) ** 2 + (yA - yP) ** 2)

            # Détermination des vecteurs
            n = np.array([dy / dA, -dx / dA])
            eKSI = np.array([(xA - xP) / dKSI, (yA - yP) / dKSI])
            eETA = np.array([dx / dA, dy / dA])

            return dA, dKSI, n, eKSI, eETA

        # Parcours les faces sur les conditions frontières et remplis la matrice A et le vecteur B
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces()):
            tag = self.mesh_obj.get_boundary_face_to_tag(i_face)   # Numéro de la frontière de la face
            bc_type, bc_value = self.bcdata[tag]                   # Condition frontière (Dirichlet ou Neumann)
            nodes = self.mesh_obj.get_face_to_nodes(i_face)        # Noeuds de la face
            element = self.mesh_obj.get_face_to_elements(i_face)[0]  # Élément de la face

            # Détermine la position du centre de la face
            ptA = ((self.mesh_obj.get_node_to_xycoord(nodes[0])[0] + self.mesh_obj.get_node_to_xycoord(nodes[1])[0])/2,
                  (self.mesh_obj.get_node_to_xycoord(nodes[0])[1] + self.mesh_obj.get_node_to_xycoord(nodes[1])[1])/2)

            dA, dKSI, n, eKSI, eETA = \
                compute_lengths_and_unit_vectors(pta=self.mesh_obj.get_node_to_xycoord(nodes[0]),
                                                 ptb=self.mesh_obj.get_node_to_xycoord(nodes[1]),
                                                 ptA=ptA,
                                                 ptP=self.centroids[element])

            # Calcule les projections de vecteurs unitaires
            PNKSI = np.dot(n, eKSI)       # Projection de n sur ξ
            PKSIETA = np.dot(eKSI, eETA)  # Projection de ξ sur η

            if bc_type == "DIRICHLET":
                D = (1 / PNKSI) * self.case.get_gamma() * (dA / dKSI)  # Direct gradient term

                ###### Calculer le cross-diffusion term ici pour dirichlet ici...
                # Calcule le terme correction de cross-diffusion si activé
                if self.cross_diffusion is True:
                    Sdc = 0

                A[element, element] += D
                B[element] += D * bc_value + Sdc
            elif bc_type == "NEUMANN":
                B[element] += self.case.get_gamma() * bc_value * dA

        # Parcours les faces internes et remplis la matrice A et le vecteur B
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces(), self.mesh_obj.get_number_of_faces()):
            nodes = self.mesh_obj.get_face_to_nodes(i_face)
            elements = self.mesh_obj.get_face_to_elements(i_face)

            dA, dKSI, n, eKSI, eETA = \
                compute_lengths_and_unit_vectors(pta=self.mesh_obj.get_node_to_xycoord(nodes[0]),
                                                 ptb=self.mesh_obj.get_node_to_xycoord(nodes[1]),
                                                 ptA=self.centroids[elements[1]],
                                                 ptP=self.centroids[elements[0]])

            # Calcule les projections de vecteurs unitaires
            PNKSI = np.dot(n, eKSI)       # Projection de n sur ξ
            PKSIETA = np.dot(eKSI, eETA)  # Projection de ξ sur η

            D = (1 / PNKSI) * self.case.get_gamma() * (dA / dKSI)  # Direct gradient term

            ###### Calculer le cross-diffusion term ici pour dirichlet ici...
            # Calcule le terme correction de cross-diffusion si activé
            if self.cross_diffusion is True:
                Sdc = 0

            # Remplissage de la matrice et du vecteur
            A[elements[0], elements[0]] += D
            A[elements[1], elements[1]] += D
            A[elements[0], elements[1]] -= D
            A[elements[1], elements[0]] -= D

            B[elements[0]] += Sdc
            B[elements[1]] -= Sdc

        # Ajout de la contribution du terme source sur les éléments
        B = np.add(B, self.volumes*self.case.source_term)

        # Résolution du problème
        PHI = linsolve.spsolve(sps.csr_matrix(A, dtype=np.float64), B)

        self.case.set_solution(PHI, self.volumes)
