import numpy as np


# Solveur utilisant la méthode des volumes finis avec diffusion seulement
class MethodeVolumesFinisDiffusion:
    def __init__(self, case, cross_diffusion):
        self.case = case  # Cas à résoudre
        self.mesh_obj = case.get_mesh()         # Maillage du cas
        self.bcdata = case.get_bc()             # Conditions frontières
        self.cross_diffusion = cross_diffusion  # Si le terme de cross diffusion est activé

        self.centroids = np.zeros([self.mesh_obj.get_number_of_elements(), 2])
        self.areas = np.zeros([self.mesh_obj.get_number_of_elements(), 1])

        self.dA = np.zeros([self.mesh_obj.get_number_of_faces(), 1])
        self.dKSI = np.zeros([self.mesh_obj.get_number_of_faces(), 1])

        # Pas necessaire à mon avis
        """self.normals = np.zeros([self.mesh_obj.get_number_of_faces(), 1])
        self.eKSI = np.zeros([self.mesh_obj.get_number_of_faces(), 1, 2])
        self.eETA = np.zeros([self.mesh_obj.get_number_of_faces(), 1, 2])"""

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

            self.areas[i_elem] = np.sum([np.linalg.det(area_matrices[i]) for i in range(len(nodes))]) / 2.

        # Calcul la longueur de la face et sa normal (pas trop sure que c'est nécessaire de faire une fonction poru ca)
        """def compute_lengths_and_unit_vectors(i_face):
            nodes = self.mesh_obj.get_face_to_nodes(i_face)
            elements = self.mesh_obj.get_face_to_elements(i_face)

            (xa, ya), (xb, yb) = self.mesh_obj.get_node_to_xycoord(nodes[0]), \
                                 self.mesh_obj.get_node_to_xycoord(nodes[1])
            (xA, yA), (xP, yP) = self.centroids[elements[1]], self.centroids[elements[0]]

            self.dA[i_face] = np.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)
            self.dKSI[i_face] = np.sqrt((xA - xP) ** 2 + (yA - yP) ** 2)

            self.normals[i_face] = np.array([(yb - ya) / self.dA[i_face], -(xb - xa) / self.dA[i_face]])
            self.eKSI[i_face] = np.array([(xA - xP) / self.dKSI[i_face], (yA - yP) / self.dKSI[i_face]])
            self.eETA[i_face] = np.array([(xb - xa) / self.dA[i_face], (yb - ya) / self.dKSI[i_face]])"""



        # Calculs pour les éléments (centres d'élément et aires)
        for i_elem in range(self.mesh_obj.get_number_of_elements()):
            find_centroid(i_elem)
            compute_area(i_elem)

        # Calculs pour les faces (longueurs des face)
        """for i_face in range(self.mesh_obj.get_number_of_faces()):
            compute_lengths_and_unit_vectors(i_face)"""

    def solve(self):
        # Itinitialisation des matrices
        NELEM = self.mesh_obj.get_number_of_elements()
        A = np.zeros((NELEM, NELEM))
        B = np.zeros(NELEM)
        Sdc = 0  # Cross-diffusion term reste nul si False

        for i_face in range(self.mesh_obj.get_number_of_boundary_faces()):
            tag = self.mesh_obj.get_boundary_face_to_tag(i_face)  # Numéro de la frontière de la face
            bc_type = self.bcdata[tag][0]  # Type de condition frontière (Dirichlet ou Neumann)
            bc_value = self.bcdata[tag][1]
            nodes = self.mesh_obj.get_face_to_nodes(i_face)
            elements = self.mesh_obj.get_face_to_elements(i_face)  # Élément de la face
            (xa, ya), (xb, yb) = self.mesh_obj.get_node_to_xycoord(nodes[0]), \
                                 self.mesh_obj.get_node_to_xycoord(nodes[1])
            (xA, yA), (xP, yP) = ((xa + xb)/2, (ya + yb)/2), self.centroids[elements[0]]
            dx, dy = (xb - xa), (yb - ya)
            dA = np.sqrt(dx ** 2 + dy ** 2)
            dKSI = np.sqrt((xA - xP) ** 2 + (yA - yP) ** 2)

            n = np.array([(yb - ya) / dA, -(xb - xa) / dA])
            eKSI = np.array([(xA - xP) / dKSI, (yA - yP) / dKSI])
            eETA = np.array([(xb - xa) / dA, (yb - ya) / dA])

            PNKSI = np.dot(n, eKSI)  # Projection de n sur ξ
            PKSIETA = np.dot(eKSI, eETA)  # Projection de ξ sur η


            if bc_type == "DIRICHLET":
                D = (1/PNKSI)*self.case.get_gamma()*(dA/dKSI)  # Direct gradient term

                # Calculer le cross-diffusion term ici pour dirichlet ici...
                if self.cross_diffusion is True:
                    Sdc = 0


                A[elements[0], elements[0]] += D
                B[elements[0]] += D*bc_value + Sdc

            elif bc_type == "NEUMANN":
                B[elements[0]] += self.case.get_gamma()*bc_value*dA


        # Remplissage de la matrice et du vecteur pour les faces internes
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces(), self.mesh_obj.get_number_of_faces()):
            nodes = self.mesh_obj.get_face_to_nodes(i_face)
            elements = self.mesh_obj.get_face_to_elements(i_face)

            (xa, ya), (xb, yb) = self.mesh_obj.get_node_to_xycoord(nodes[0]), \
                                 self.mesh_obj.get_node_to_xycoord(nodes[1])
            (xA, yA), (xP, yP) = self.centroids[elements[1]], self.centroids[elements[0]]

            dx, dy = (xb - xa), (yb - ya)
            dA = np.sqrt(dx ** 2 + dy ** 2)
            dKSI = np.sqrt((xA - xP) ** 2 + (yA - yP) ** 2)

            n = np.array([(yb - ya) / dA, -(xb - xa) / dA])
            eKSI = np.array([(xA - xP) / dKSI, (yA - yP) / dKSI])
            eETA = np.array([(xb - xa) / dA, (yb - ya) / dA])

            PNKSI = np.dot(n, eKSI)  # Projection de n sur ξ
            PKSIETA = np.dot(eKSI, eETA)  # Projection de ξ sur η

            D = (1/PNKSI)*self.case.get_gamma()*(dA/dKSI)  # Direct gradient term

            # Calculer le cross-diffusion term ici pour les faces internes...
            if self.cross_diffusion is True:
                Sdc = 0


            # Remplissage de la matrice et du vecteur
            A[elements[0], elements[0]] += D
            A[elements[1], elements[1]] += D
            A[elements[0], elements[1]] -= D
            A[elements[1], elements[0]] -= D

            B[elements[0]] += Sdc
            B[elements[1]] -= Sdc

        # Ajout de la contribution du terme source
        for i_elem in range(self.mesh_obj.get_number_of_elements()):
            B[i_elem] += self.case.source_term*self.areas[i_elem]

        PHI = np.linalg.solve(A, B)

        self.case.set_solution(PHI, self.areas)