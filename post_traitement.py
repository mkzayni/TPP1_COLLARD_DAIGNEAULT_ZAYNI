"""
Date :    8 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Mohamad Karim Zayni (2167132)
Utilité : TPP1 - Méthode des volumes finis avec diffusion
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

class PostTraitement:
    def __init__(self, exemple):
        self.exemple = exemple  # Nom de l'exemple post-traité
        self.data = []  # Initialisation du dictionnaire de données

    # Ajoute les données selon un nombre d'éléments'
    def set_data(self, case):
        self.data.append({'n': case.get_mesh().get_number_of_elements(),
                          'phi_num': case.get_solutions()[0],
                          'phi_exact': case.get_solutions()[1],
                          'area': case.get_areas_and_centroids()[0],
                          'position': case.get_areas_and_centroids()[1],
                          'time': case.get_times()})

    # Génère les graphiques des solutions numérique et analytique
    def show_solutions(self, mesh, title, save_path):
        Figure1, (NUM, EX) = plt.subplots(1, 2, figsize=(20, 8))

        Figure1.suptitle(title)

        # Set levels of color for the colorbar
        levels = np.linspace(np.min([self.data[mesh]['phi_num'], self.data[mesh]['phi_exact']]),
                             np.max([self.data[mesh]['phi_num'], self.data[mesh]['phi_exact']]), num=25)

        c = NUM.tricontourf(self.data[mesh]['position'][:, 0],
                            self.data[mesh]['position'][:, 1],
                            self.data[mesh]['phi_num'], levels=levels)
        plt.colorbar(c, ax=NUM)
        NUM.set_xlabel("L (m)")
        NUM.set_ylabel("H (m)")
        NUM.set_title("Solution numérique")

        c = EX.tricontourf(self.data[mesh]['position'][:, 0],
                           self.data[mesh]['position'][:, 1],
                           self.data[mesh]['phi_exact'], levels=levels)
        plt.colorbar(c, ax=EX)
        EX.set_xlabel("L (m)")
        EX.set_ylabel("H (m)")
        EX.set_title("Solution analytique MMS")

        plt.savefig(save_path, dpi=200)

    def show_plan_solutions(self, mesh, title, save_path, X_Coupe, Y_Coupe):
        # Chercher l'indice des éléments à un X donné
        def Coupe_X(Coordonnees, X, Solution, Analytique):
            Elements_ds_coupe = []
            Solution_coupe = []
            Analytique_coupe = []
            eps = 1e-6  # Précision
            for i in range(len(Coordonnees)):
                if np.abs(Coordonnees[i, 0] - X) < eps:
                    Elements_ds_coupe.append(Coordonnees[i, :])
                    Solution_coupe.append(Solution[i])
                    Analytique_coupe.append(Analytique[i])
            Elements_ds_coupe = np.array(Elements_ds_coupe)
            Solution_coupe = np.array(Solution_coupe)
            return Elements_ds_coupe, Solution_coupe, Analytique_coupe

        def Coupe_Y(Coordonnees, Y, Solution, Analytique):
            Elements_ds_coupe = []
            Solution_coupe = []
            Analytique_coupe = []
            eps = 1e-6  # Précision
            for i in range(len(Coordonnees)):
                if (np.abs(Coordonnees[i, 1] - Y) < eps):
                    Elements_ds_coupe.append(Coordonnees[i, :])
                    Solution_coupe.append(Solution[i])
                    Analytique_coupe.append(Analytique[i])
            Elements_ds_coupe = np.array(Elements_ds_coupe)
            Solution_coupe = np.array(Solution_coupe)
            return Elements_ds_coupe, Solution_coupe, Analytique_coupe

        Figure1, (COUPEX, COUPEY) = plt.subplots(1, 2, figsize=(20, 8))

        Figure1.suptitle(title)

        Centres = self.data[mesh]['position']

        Elem_ds_coupeX, Solution_coupeX, SolutionEX_coupeX = \
            Coupe_X(Centres, X_Coupe, self.data[mesh]['phi_num'], self.data[mesh]['phi_exact'])
        Elem_ds_coupeY, Solution_coupeY, SolutionEX_coupeY = \
            Coupe_Y(Centres, Y_Coupe, self.data[mesh]['phi_num'], self.data[mesh]['phi_exact'])

        COUPEX.plot(Solution_coupeX, Elem_ds_coupeX[:, 1], label="Solution Numérique")
        COUPEX.plot(SolutionEX_coupeX, Elem_ds_coupeX[:, 1], label="Solution MMS")
        COUPEX.set_xlabel("Température")
        COUPEX.set_ylabel("Y (m)")
        COUPEX.set_title(f"Solution dans une coupe à X = {X_Coupe}")
        COUPEX.legend()

        COUPEY.plot(Elem_ds_coupeY[:, 0], Solution_coupeY, label="Solution Numérique")
        COUPEY.plot(Elem_ds_coupeY[:, 0], SolutionEX_coupeY, label="Solution MMS")
        COUPEY.set_xlabel("X (m)")
        COUPEY.set_ylabel("Température")
        COUPEY.set_title(f"Solution dans une coupe à Y = {Y_Coupe}")
        COUPEY.legend()

        # Enregistrer
        plt.savefig(save_path, dpi=200)

    def show_mesh_differences(self, mesh1, mesh2, title, save_path):
        figure, (plot1, plot2) = plt.subplots(1, 2, figsize=(20, 8))

        figure.suptitle(title)

        # Set levels of color for the colorbar
        levels = np.linspace(np.min(np.append(self.data[mesh1[0]]['phi_num'], self.data[mesh2[0]]['phi_num'])),
                             np.max(np.append(self.data[mesh1[0]]['phi_num'], self.data[mesh2[0]]['phi_num'])), num=25)

        center1 = self.data[mesh1[0]]['position']
        c = plot1.tricontourf(center1[:, 0], center1[:, 1], self.data[mesh1[0]]['phi_num'], levels=levels)
        plot1.set_xlabel("L (m)")
        plot1.set_ylabel("H (m)")
        plot1.set_title(f"{mesh1[1]} à {self.data[mesh1[0]]['n']} éléments")
        plt.colorbar(c, ax=plot1)

        center2 = self.data[mesh2[0]]['position']
        c = plot2.tricontourf(center2[:, 0], center2[:, 1], self.data[mesh2[0]]['phi_num'], levels=levels)
        plot2.set_xlabel("L (m)")
        plot2.set_ylabel("H (m)")
        plot2.set_title(f"{mesh2[1]} à {self.data[mesh2[0]]['n']} éléments")
        plt.colorbar(c, ax=plot2)


        # Enregistrer
        plt.savefig(save_path, dpi=200)

    def show_time(self, save_path):
        Figure1, (Comp) = plt.subplots(1, 1, figsize=(15, 10))

        Figure1.suptitle("Comparaison Temps de calculs")

        taille = np.zeros(len(self.data))
        Time_Dense = np.zeros(len(self.data))
        Time_Sparse = np.zeros(len(self.data))
        for i in range(len(self.data)):
            taille[i] = self.data[i]['n']**2
            Time_Dense[i] = self.data[i]['time']['DENSE']
            Time_Sparse[i] = self.data[i]['time']['SPARSE']

        Comp.plot(taille,Time_Sparse,'.-',label="Méthode Sparse")
        Comp.plot(taille,Time_Dense,'.-',label="Méthode Dense")
        Comp.set_xlabel("Taille de la matrice (Nombre d'éléments)²")
        Comp.set_ylabel("Temps de résolution (s)")
        Comp.grid()
        Comp.legend()


        plt.savefig(save_path, dpi=200)

    def show_error(self):
        # Calcul l'erreur (en x), l'ajoute aux données et détermine l'ordre de convergence
        for i in range(len(self.data)):
            total_area = np.sum(self.data[i]['area'])
            area, n = self.data[i]['area'], self.data[i]['n']
            phi_num, phi_exact = self.data[i]['phi_num'], self.data[i]['phi_exact']
            E_L2 = np.sqrt(np.sum(area*(phi_num - phi_exact)**2)/total_area)
            self.data[i]['err_L2'] = E_L2
            self.data[i]['h'] = np.sqrt(total_area/n)

        p = np.polyfit(np.log([self.data[i]['h'] for i in range(len(self.data))]),
                  np.log([self.data[i]['err_L2'] for i in range(len(self.data))]), 1)

        # Graphique de l'erreur
        fig_E, ax_E = plt.subplots(figsize=(15, 10))
        fig_E.suptitle("Normes de l'erreur L² des solutions numériques sur une échelle de logarithmique "
                       "pour " + self.exemple, y=0.925)
        text = AnchoredText('Ordre de convergence: ' + str(round(p[0], 2)), loc='upper left')

        ax_E.loglog([self.data[i]['h'] for i in range(len(self.data))],
                  [self.data[i]['err_L2'] for i in range(len(self.data))], '.-')
        ax_E.minorticks_on()
        ax_E.grid(True, which="both", axis="both", ls="-")
        ax_E.set_xlabel('Grandeur (h)')
        ax_E.set_ylabel('Erreur (E)')
        ax_E.add_artist(text)