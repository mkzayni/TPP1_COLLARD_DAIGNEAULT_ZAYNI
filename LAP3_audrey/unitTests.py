import numpy as np

class UnitTest:
    def __init__(self, case):
        self.case = case

    def test_relation_euler(self):
        nb_face = self.case.get_mesh().get_number_of_elements()
        nb_arete = self.case.get_mesh().get_number_of_faces()
        nb_sommet = self.case.get_mesh().get_number_of_nodes()
        nb_trou = self.case.get_nb_trou()

        if (nb_face - nb_arete + nb_sommet) == (1 - nb_trou):
            print("Vérification de la relation d'Euler : PASS")
        else:
            print("Vérification de la relation d'Euler : FAIL")

    def test_divergence(self):
        mesh_obj = self.case.get_mesh()
        F = np.array([5, 3])

        number_of_elements = mesh_obj.get_number_of_elements()
        div = np.zeros(number_of_elements)

        # Boundary faces
        for i_face in range(mesh_obj.get_number_of_faces()):
            element = mesh_obj.get_face_to_elements(i_face)
            nodes = mesh_obj.get_face_to_nodes(i_face)

            positions = [[mesh_obj.get_node_to_xcoord(nodes[0]), mesh_obj.get_node_to_ycoord(nodes[0])],
                         [mesh_obj.get_node_to_xcoord(nodes[1]), mesh_obj.get_node_to_ycoord(nodes[1])]]

            ds = np.sqrt((positions[0][0] - positions[1][0])**2 + (positions[0][1] - positions[1][1])**2)

            A = np.array([[positions[0][0] - positions[1][0],   positions[0][1] - positions[1][1]],
                          [positions[0][1] - positions[1][1], -(positions[0][0] - positions[1][0])]])
            b = np.array([0, ds])
            n = np.linalg.solve(A, b)

            # Ajout du flux à l'élément de gauche
            flux = np.dot(F, n) * ds
            div[element[0]] += flux

            # Si élément de droite existe soustraire le flux
            if element[1] != -1:
                flux = np.dot(F, n) * ds
                div[element[1]] -= flux


        if np.sum(div) < 1e-6:
            print("Vérification de la divergence d'un champ constant : PASS")
        else:
            print("Vérification de la divergence d'un champ constant : FAIL")

