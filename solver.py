import numpy as np


# Solveur utilisant la méthode des volumes finis avec diffusion seulement
class MethodeVolumesFinisDiffusion:
    def __init__(self, case):
        self.case = case  # Cas à résoudre
        self.mesh_obj = case.get_mesh()  # Maillage du cas
        self.bcdata = case.get_bc()  # Conditions frontières

        self.centroids = np.zeros([self.mesh_obj.get_number_of_elements(), 1, 2])
        self.find_centroids()



        # Appeler les fonctions permettant de calculer les informations suivantes :
        #   - la longueur des faces
        #   - Calculer l'aire des éléments
        #   - Déterminer le centre des faces ...

    # Détermine les positions centrales des éléments touchant l'arête
    def find_centroids(self):
        for i_elem in range(self.mesh_obj.get_number_of_elements()):
            # Détermine les numéros des noeuds pour l'élément
            start, end = self.mesh_obj.element_to_nodes_start[i_elem], self.mesh_obj.element_to_nodes_start[i_elem + 1]
            nodes = self.mesh_obj.element_to_nodes[start:end]

            # Détermine les centroides des éléments par la moyenne des coordonnées
            for i_node in nodes:
                self.centroids[i_elem] += [self.mesh_obj.get_node_to_xycoord(i_node)[0],
                                self.mesh_obj.get_node_to_xycoord(i_node)[1]]

            self.centroids[i_elem] /= len(nodes)


    """# Calcule les distances dx et dy entre les centres des éléments d'une arête ou entre
    def compute_distances(self, i_face):
        elements = self.mesh_obj.get_face_to_elements(i_face)  # Liste des éléments touchant l'arête

        # Dans le cas d'une arête interne (2 éléments)
        if elements[1] != -1:
            # Détermine les numéros des noeuds pour les éléments de gauche et de droite sous forme de 2 listes
            start = [self.mesh_obj.element_to_nodes_start[elements[0]],
                     self.mesh_obj.element_to_nodes_start[elements[1]]]
            end = [self.mesh_obj.element_to_nodes_start[elements[0] + 1],
                   self.mesh_obj.element_to_nodes_start[elements[1] + 1]]
            nodes = [self.mesh_obj.element_to_nodes[start[0]:end[0]], self.mesh_obj.element_to_nodes[start[1]:end[1]]]
        # Dans le cas d'une arête touchant une frontière (1 élément)
        else:
            # Détermine les numéros des noeuds pour l'élément de gauche et de l'arête sous forme de 2 listes
            start = self.mesh_obj.element_to_nodes_start[elements[0]]
            end = self.mesh_obj.element_to_nodes_start[elements[0] + 1]
            nodes = [self.mesh_obj.element_to_nodes[start:end], self.mesh_obj.get_face_to_nodes(i_face)]

        # Détermine les positions centrales des éléments et la distance entre celles-ci
        # center positions = (sum(xi)/nb_nodes et sum(yi)/nb_nodes) pour les 2 points à déterminer
        # Calcule l'aire de l'élément de gauche
        center_positions = np.zeros([2, 2])
        area_matrix = [np.zeros([2, 2]) for i in range(len(nodes[0]))]

        for i in range(len(nodes[0])):
            i_nodes_left = nodes[0][i]
            x, y = self.mesh_obj.get_node_to_xycoord(i_nodes_left)[0], self.mesh_obj.get_node_to_xycoord(i_nodes_left)[
                1]
            center_positions[0] += [x, y]

            # Construit les matrices pour calculer l'aire
            area_matrix[i][:][0] = [x, y]
            area_matrix[i - 1][:][1] = [x, y]

        for i_nodes_right in nodes[1]:
            center_positions[1] += [self.mesh_obj.get_node_to_xycoord(i_nodes_right)[0],
                                    self.mesh_obj.get_node_to_xycoord(i_nodes_right)[1]]

        nb_nodes = np.array([len(nodes[0]), len(nodes[1])])
        center_positions[0], center_positions[1] = center_positions[0] / nb_nodes[0], center_positions[1] / nb_nodes[1]

        dx = center_positions[1][0] - center_positions[0][0]  # Xd - Xg
        dy = center_positions[1][1] - center_positions[0][1]  # Yd - Yg

        # Aire pour TRI ou QUAD de l'élément de gauche
        self.areas[elements[0]] = np.sum([np.linalg.det(area_matrix[i]) for i in range(len(nodes[0]))]) / 2

        return dx, dy, center_positions

    # Calcule le vecteur normal au centre d'une arête
    def compute_normal(self, i_face):
        nodes = self.mesh_obj.get_face_to_nodes(i_face)
        (xa, ya), (xb, yb) = self.mesh_obj.get_node_to_xycoord(nodes[0]), self.mesh_obj.get_node_to_xycoord(nodes[1])
        dA = np.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)
        normal = np.array([(yb - ya) / dA, -(xb - xa) / dA])

        return normal

    # Calcule le gradient du cas étudié
    def solve(self):
        # Itinitialisation des matrices
        NTRI = self.mesh_obj.get_number_of_elements()
        ATA = np.zeros((NTRI, 2, 2))
        B = np.zeros((NTRI, 2))
        dphi_exact = np.zeros((NTRI, 2))
        self.areas = np.zeros(NTRI)

        # Remplissage des matrices pour le cas d'une condition frontière (Dirichlet ou Neumann)
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces()):
            tag = self.mesh_obj.get_boundary_face_to_tag(i_face)  # Numéro de la frontière de la face
            bc_type = self.bcdata[tag][0]  # Type de condition frontière (Dirichlet ou Neumann)
            elements = self.mesh_obj.get_face_to_elements(i_face)  # Élément de la face

            # Si Dirichlet ou Neumann
            if bc_type != 'LIBRE':
                # Détermination des positions des points et de la distance
                dx, dy, center_positions = self.compute_distances(i_face)
                Xtg, Ytg, Xta, Yta = center_positions[0][0], center_positions[0][1], center_positions[1][0], \
                                     center_positions[1][1]
                # Calcul du gradient
                dphi = self.phi(Xta, Yta) - self.phi(Xtg, Ytg)

                # Modification de la position du point sur l'arête si Neumann
                if bc_type == 'NEUMANN':
                    n = self.compute_normal(i_face)
                    dx, dy = np.dot([dx, dy], n) * n
                    Xan, Yan = np.array([Xtg + dx, Ytg + dy])
                    dphi = self.phi(Xan, Yan) - self.phi(Xtg, Ytg)

                ALS = np.array([[dx * dx, dx * dy], [dy * dx, dy * dy]])
                ATA[elements[0]] += ALS

                # Remplisage du membre de droite
                B[elements[0]] += (np.array([dx, dy]) * dphi)

                # Calcule la fonction analytique pour l'élément de gauche
                dphi_exact[elements[0]] = [self.dphi_function[0](Xtg), self.dphi_function[1](Ytg)]

        # Internal faces
        for i_face in range(self.mesh_obj.get_number_of_boundary_faces(), self.mesh_obj.get_number_of_faces()):
            elements = self.mesh_obj.get_face_to_elements(i_face)
            dx, dy, center_positions = self.compute_distances(i_face)

            # Remplissage de la matrice ATA pour l'arête interne
            ALS = np.array([[dx * dx, dx * dy], [dy * dx, dy * dy]])
            ATA[elements[0]] += ALS
            ATA[elements[1]] += ALS

            # Remplisage du membre de droite
            Xtg, Ytg, Xtd, Ytd = center_positions[0][0], center_positions[0][1], center_positions[1][0], \
                                 center_positions[1][1]
            dphi = self.phi(Xtd, Ytd) - self.phi(Xtg, Ytg)
            B[elements[0]] += (np.array([dx, dy]) * dphi)
            B[elements[1]] += (np.array([dx, dy]) * dphi)

            # Calcule la fonction analytique pour l'élément de gauche
            dphi_exact[elements[0]] = [self.dphi_function[0](Xtg), self.dphi_function[1](Ytg)]

        ATAI = np.array([np.linalg.inv(ATA[i_tri]) for i_tri in range(NTRI)])
        GRAD = np.array([np.dot(ATAI[i_tri], B[i_tri]) for i_tri in range(NTRI)])

        self.case.set_solutions(GRAD, dphi_exact, self.areas)
"""