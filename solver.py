import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve


# Solveur utilisant la méthode des volumes finis avec diffusion seulement
class MethodeVolumesFinisDiffusion:
    def __init__(self, case, cross_diffusion):
        self.case = case  # Cas à résoudre
        self.mesh_obj = case.get_mesh()  # Maillage du cas
        self.bcdata = case.get_bc()  # Conditions frontières
        self.cross_diffusion = cross_diffusion  # Si le terme de cross diffusion est activé

        self.centroids = np.zeros([self.mesh_obj.get_number_of_elements(), 2])  # Array des centroides
        self.volumes = np.zeros([self.mesh_obj.get_number_of_elements()])       # Array de l'aire des elements

        self.preprocessing()

    # Effectue les calculs relatifs au maillage préalablement à l'utilisation du solver
    def preprocessing(self):
        # Détermine les centroides et l'aire de l'élément par les déterminants
        def compute_centroid_and_volume(i_elem):
            nodes = self.mesh_obj.get_element_to_nodes(i_elem)
            area_matrices = [np.zeros([2, 2]) for i in range(len(nodes))]
            for i in range(len(nodes)):
                x, y = self.mesh_obj.get_node_to_xycoord(nodes[i])[0], self.mesh_obj.get_node_to_xycoord(nodes[i])[1]
                area_matrices[i][:, 0] =   [x, y]
                area_matrices[i-1][:, 1] = [x, y]

            # Calcule l'aire de l'élément
            self.volumes[i_elem] = np.sum([np.linalg.det(area_matrices[i]) for i in range(len(nodes))]) / 2
            cx = (np.sum(
                [np.sum(area_matrices[i][0, :]) * np.linalg.det(area_matrices[i]) for i in range(len(nodes))]) /
                  (6 * self.volumes[i_elem]))
            cy = (np.sum(
                [np.sum(area_matrices[i][1, :]) * np.linalg.det(area_matrices[i]) for i in range(len(nodes))]) /
                  (6 * self.volumes[i_elem]))

            self.centroids[i_elem] = [cx, cy]

            """N_nodes=len(nodes)
            #Cas Triangle
            #if(N_nodes==3):
            for i_node in nodes:
                x, y = self.mesh_obj.get_node_to_xycoord(i_node)[0], self.mesh_obj.get_node_to_xycoord(i_node)[1]
                self.centroids[i_elem][0] += x
                self.centroids[i_elem][1] += y

            self.centroids[i_elem] /= len(nodes)

            #Cas QUAD
            elif(N_nodes>3):
                #Calcul de A
                A=0
                for i_node in nodes-1:
                    x1, y1 = self.mesh_obj.get_node_to_xycoord(i_node)[0], self.mesh_obj.get_node_to_xycoord(i_node)[1]
                    x2, y2 = self.mesh_obj.get_node_to_xycoord(i_node+1)[0], self.mesh_obj.get_node_to_xycoord(i_node+1)[1]
                    A+=0.5*(x1*y2-x2*y1)
                #Calcul coordonnes Centroide
                for i_node in nodes-1:
                    x1, y1 = self.mesh_obj.get_node_to_xycoord(i_node)[0], self.mesh_obj.get_node_to_xycoord(i_node)[1]
                    x2, y2 = self.mesh_obj.get_node_to_xycoord(i_node+1)[0], self.mesh_obj.get_node_to_xycoord(i_node+1)[1]
                    self.centroids[i_elem][0] += (1.0/6.0*A)*(x1+x2)*(x1*y2-x2*y1)
                    self.centroids[i_elem][1] += (1.0/6.0*A)*(y1+y2)*(x1*y2-x2*y1)
                print("centroide")
                print(self.centroids[i_elem])
                """

        # Calculs pour les éléments (centres d'élément et aires)
        for i_elem in range(self.mesh_obj.get_number_of_elements()):
            compute_centroid_and_volume(i_elem)

    def solve(self):
        # Itinitialisation des matrices et du terme de cross-diffusion
        NELEM = self.mesh_obj.get_number_of_elements()
        A = np.zeros((NELEM, NELEM))
        B = np.zeros(NELEM)
        PHI = np.zeros(NELEM)
        PHI_EX = np.zeros(NELEM)
        GRAD = np.zeros((NELEM, 2))
        gamma = self.case.get_gamma()
        Sdc = 0  # Cross-diffusion term reste nul si False
        it = 0

        solver_moindrescarres = GradientMoindresCarres(self.case)
        solver_moindrescarres.set_centroids_and_volumes(self.centroids, self.volumes)

        if self.cross_diffusion is True:
            it = 3

        for i in range(it):
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
                tag = self.mesh_obj.get_boundary_face_to_tag(i_face)     # Numéro de la frontière de la face
                bc_type, bc_value = self.bcdata[tag]                     # Condition frontière (Dirichlet ou Neumann)
                nodes = self.mesh_obj.get_face_to_nodes(i_face)          # Noeuds de la face
                element = self.mesh_obj.get_face_to_elements(i_face)[0]  # Élément de la face

                # Détermine la position du centre de la face
                (xa, ya) = (
                (self.mesh_obj.get_node_to_xycoord(nodes[0])[0] + self.mesh_obj.get_node_to_xycoord(nodes[1])[0]) / 2,
                (self.mesh_obj.get_node_to_xycoord(nodes[0])[1] + self.mesh_obj.get_node_to_xycoord(nodes[1])[1]) / 2)

                dA, dKSI, n, eKSI, eETA = \
                    compute_lengths_and_unit_vectors(pta=self.mesh_obj.get_node_to_xycoord(nodes[0]),
                                                     ptb=self.mesh_obj.get_node_to_xycoord(nodes[1]),
                                                     ptA=(xa, ya),
                                                     ptP=self.centroids[element])

                # Calcule les projections de vecteurs unitaires
                dETA = dA  # Équivalent, mais noté pour éviter la confusion
                PNKSI = np.dot(n, eKSI)  # Projection de n sur ξ
                PKSIETA = np.dot(eKSI, eETA)  # Projection de ξ sur η

                if bc_type == "DIRICHLET":
                    D = (1 / PNKSI) * gamma * (dA / dKSI)  # Direct gradient term

                    # Calcule le terme correction de cross-diffusion si activé
                    if self.cross_diffusion is True:
                        # Évaluation des phi aux noeuds de la face frontière
                        phi0 = bc_value(self.mesh_obj.get_node_to_xycoord(nodes[0])[0],
                                        self.mesh_obj.get_node_to_xycoord(nodes[0])[1])
                        phi1 = bc_value(self.mesh_obj.get_node_to_xycoord(nodes[1])[0],
                                        self.mesh_obj.get_node_to_xycoord(nodes[1])[1])
                        Sdc = -gamma * (PKSIETA / PNKSI) * (phi1 - phi0) / dETA * dA

                    A[element, element] += D
                    B[element] += D * bc_value(xa, ya) + Sdc
                elif bc_type == "NEUMANN":
                    B[element] += gamma * bc_value(xa, ya) * dA

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
                PNKSI = np.dot(n, eKSI)  # Projection de n sur ξ
                PKSIETA = np.dot(eKSI, eETA)  # Projection de ξ sur η

                D = (1 / PNKSI) * self.case.get_gamma() * (dA / dKSI)  # Direct gradient term

                ###### Calculer le cross-diffusion term ici pour dirichlet ici...
                # Calcule le terme correction de cross-diffusion si activé
                if self.cross_diffusion is True:
                    Sdc = -gamma * (PKSIETA / PNKSI) * np.dot((GRAD[elements[1]] + GRAD[elements[0]]) / 2, eETA) * dA

                # Remplissage de la matrice et du vecteur
                A[elements[0], elements[0]] += D
                A[elements[1], elements[1]] += D
                A[elements[0], elements[1]] -= D
                A[elements[1], elements[0]] -= D

                B[elements[0]] += Sdc
                B[elements[1]] -= Sdc

            # Ajout de la contribution du terme source sur les éléments et calcul de la solution analytique
            for i_elem in range(self.mesh_obj.get_number_of_elements()):
                B[i_elem] += self.volumes[i_elem] * self.case.source_term(self.centroids[i_elem][0],
                                                                          self.centroids[i_elem][1])
                PHI_EX[i_elem] = self.case.get_analytical_function()(self.centroids[i_elem][0],
                                                                     self.centroids[i_elem][1])

            # Résolution du problème
            PHI = linsolve.spsolve(sps.csr_matrix(A, dtype=np.float64), B)

            solver_moindrescarres.set_phi(PHI)
            GRAD = solver_moindrescarres.solve()

        self.case.set_solution(PHI, self.volumes, PHI_EX)


class GradientMoindresCarres:
    def __init__(self, case):
        self.case = case                 # Cas à résoudre
        self.mesh_obj = case.get_mesh()  # Maillage du cas
        self.bcdata = case.get_bc()      # Types de conditions frontière

    def set_centroids_and_volumes(self, centroids, volumes):
        self.centroids = centroids
        self.volumes = volumes

    def set_phi(self, phi):
        self.phi = phi

    # Calcule le gradient du cas étudié
    def solve(self):
        # Itinitialisation des matrices
        NTRI = self.mesh_obj.get_number_of_elements()
        ATA = np.zeros((NTRI, 2, 2))
        B = np.zeros((NTRI, 2))

        # Remplissage des matrices pour le cas d'une condition frontière (Dirichlet ou Neumann)
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces()):
            tag = self.mesh_obj.get_boundary_face_to_tag(i_face)  # Numéro de la frontière de la face
            bc_type, bc_value = self.bcdata[tag]  # Condition frontière (Dirichlet ou Neumann)
            element = self.mesh_obj.get_face_to_elements(i_face)[0]  # Élément de la face

            # Si Dirichlet ou Neumann
            # Détermination des positions des points et de la distance

            nodes = self.mesh_obj.get_face_to_nodes(i_face)
            xa = (self.mesh_obj.get_node_to_xycoord(nodes[0])[0] + self.mesh_obj.get_node_to_xycoord(nodes[1])[0]) / 2.
            ya = (self.mesh_obj.get_node_to_xycoord(nodes[0])[1] + self.mesh_obj.get_node_to_xycoord(nodes[1])[1]) / 2.
            xb, yb = self.centroids[element][0], self.centroids[element][1]
            dx, dy = xb - xa, yb - ya

            # Calcul du gradient
            if bc_type == 'DIRICHLET':
                dphi = bc_value(xa, ya) - self.phi[element]

            # Modification de la position du point sur l'arête si Neumann
            if bc_type == 'NEUMANN':
                nodes = self.mesh_obj.get_face_to_nodes(i_face)
                (xa, ya), (xb, yb) = self.mesh_obj.get_node_to_xycoord(nodes[0]), self.mesh_obj.get_node_to_xycoord(
                    nodes[1])
                dA = np.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)
                n = np.array([(yb - ya) / dA, -(xb - xa) / dA])

                dx, dy = np.dot([dx, dy], n) * n
                dphi = np.dot([dx, dy], n) * bc_value(xa, ya)

            ALS = np.array([[dx * dx, dx * dy], [dy * dx, dy * dy]])
            ATA[element] += ALS

            # Remplisage du membre de droite
            B[element] += (np.array([dx, dy]) * dphi)

        # Internal faces
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces(), self.mesh_obj.get_number_of_faces()):
            elements = self.mesh_obj.get_face_to_elements(i_face)
            dx, dy = self.centroids[elements[1]] - self.centroids[elements[0]]

            # Remplissage de la matrice ATA pour l'arête interne
            ALS = np.array([[dx * dx, dx * dy], [dy * dx, dy * dy]])
            ATA[elements[0]] += ALS
            ATA[elements[1]] += ALS

            # Remplisage du membre de droite
            dphi = self.phi[elements[0]] - self.phi[elements[1]]
            B[elements[0]] += (np.array([dx, dy]) * dphi)
            B[elements[1]] += (np.array([dx, dy]) * dphi)

        ATAI = np.array([np.linalg.inv(ATA[i_tri]) for i_tri in range(NTRI)])
        GRAD = np.array([np.dot(ATAI[i_tri], B[i_tri]) for i_tri in range(NTRI)])

        return GRAD
