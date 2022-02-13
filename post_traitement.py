"""
Date :    2 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Zacharie Prigent (2096785)
Utilité : Post-traitement pour afficher les résultats numériques et analytiques
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

class PostTraitement:
    def __init__(self, exemple):
        self.exemple = exemple  # Nom de l'exemple post-traité
        self.data = []  # Initialisation du dictionnaire de données

    # Ajoute les données selon un nombre de noeuds
    def set_data(self, case):
        self.data.append({'n': case.get_mesh().get_number_of_nodes(),
                          'phi_num': case.get_solutions()[0],
                          'phi_exact': case.get_solutions()[2],
                          'area': case.get_solutions()[1]})

    # Génère les graphiques demandés
    def genere_graphiques(self):
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
        ax_E.grid(True,which="both", axis="both",ls="-")
        ax_E.set_xlabel('Grandeur (h)')
        ax_E.set_ylabel('Erreur (E)')
        ax_E.add_artist(text)