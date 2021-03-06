"""
Date :    8 février 2022
Auteurs : Audrey Collard-Daigneault (1920374) & Mohamad Karim Zayni (2167132)
Utilité : Effectuer le post-traitement des données
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

class PostTraitement:
    """
    Effectuer le post_traitement après la résolution du problème.

    Parmeters
    ---------
    exempmle: case
        L'exemple en cours de traitement


    Attributes
    ----------
    data: dictionniare
        Dictionnaire qui comporte toutes les données de nécessaire
        pour effectuer le post-traitement.

    """
    def __init__(self, exemple):
        self.exemple = exemple  # Nom de l'exemple post-traité
        self.data = []          # Initialisation du dictionnaire de données

    # Ajoute les données selon un nombre d'éléments'
    def set_data(self, case):
        """
        Modifier les données de l'exemple traité.

        Parameters
        ----------
        case: case
            Exemple Traité. 

        Returns
        -------
        None

        """
        self.data.append({'n': case.get_mesh().get_number_of_elements(),
                          'phi_num': case.get_solutions()[0],
                          'phi_exact': case.get_solutions()[1],
                          'area': case.get_areas_and_centroids()[0],
                          'position': case.get_areas_and_centroids()[1],
                          'time': case.get_times()})

    # Génère les graphiques des solutions numérique et analytique
    def show_solutions(self, mesh, title, save_path):
        """
        Affichage des graphiques qui montrent la différence entre la solution
        numérique et la solution analytique

        Parameters
        ----------
        mesh: mesh
            Maillage de l'exemple traité. 
        
        title: str
            Nom du document pour le sauvegarder
            
        save_path: str
            Nom du fichier de sauvegarde

        Returns
        -------
        None

        """
        Figure1, (NUM, EX) = plt.subplots(1, 2, figsize=(20, 8))

        Figure1.suptitle(title + f" avec {self.data[mesh]['n']} éléments")

        # Set levels of color for the colorbar
        levels = np.linspace(np.min([self.data[mesh]['phi_num'], self.data[mesh]['phi_exact']]),
                             np.max([self.data[mesh]['phi_num'], self.data[mesh]['phi_exact']]), num=30)

        # Solution numérique
        c = NUM.tricontourf(self.data[mesh]['position'][:, 0],
                            self.data[mesh]['position'][:, 1],
                            self.data[mesh]['phi_num'], levels=levels)
        plt.colorbar(c, ax=NUM)
        NUM.set_xlabel("L (m)")
        NUM.set_ylabel("H (m)")
        NUM.set_title("Solution numérique")

        # Solution analytique/MMS
        c = EX.tricontourf(self.data[mesh]['position'][:, 0],
                           self.data[mesh]['position'][:, 1],
                           self.data[mesh]['phi_exact'], levels=levels)
        plt.colorbar(c, ax=EX)
        EX.set_xlabel("L (m)")
        EX.set_ylabel("H (m)")
        EX.set_title("Solution analytique MMS/analytique")

        plt.savefig(save_path, dpi=200)

    def show_plan_solutions(self, mesh, title, save_path, X_Coupe, Y_Coupe):
        
        """
        Affichage des graphiques qui montrent les résultats dans des coupes 
        en X ou en Y

        Parameters
        ----------
        mesh: mesh
            Maillage de l'exemple traité. 
        
        title: str
            Nom du document pour le sauvegarder
            
        save_path: str
            Nom du fichier de sauvegarde
            
        X_coupe: float
            L'endroit de la coupe suivant la droite X=X_coupe

        Y_coupe: float
            L'endroit de la coupe suivant la droite Y=Y_coupe
        Returns
        -------
        None

        """
        # Chercher l'indice des éléments à un X ou Y donné
        def Coupe_X(Coordonnees, X, Solution, Analytique, Plan):
            Elements_ds_coupe = []
            Solution_coupe = []
            Analytique_coupe = []
            eps = 1e-6  # Précision
            for i in range(len(Coordonnees)):
                if np.abs(Coordonnees[i, Plan] - X) < eps:
                    Elements_ds_coupe.append(Coordonnees[i, :])
                    Solution_coupe.append(Solution[i])
                    Analytique_coupe.append(Analytique[i])
            Elements_ds_coupe = np.array(Elements_ds_coupe)
            Solution_coupe = np.array(Solution_coupe)
            return Elements_ds_coupe, Solution_coupe, Analytique_coupe

        Figure1, (COUPEX, COUPEY) = plt.subplots(1, 2, figsize=(20, 6))

        Figure1.suptitle(title)

        Centres = self.data[mesh]['position']

        Elem_ds_coupeX, Solution_coupeX, SolutionEX_coupeX = \
            Coupe_X(Centres, X_Coupe, self.data[mesh]['phi_num'], self.data[mesh]['phi_exact'], 0)
        Elem_ds_coupeY, Solution_coupeY, SolutionEX_coupeY = \
            Coupe_X(Centres, Y_Coupe, self.data[mesh]['phi_num'], self.data[mesh]['phi_exact'], 1)

        COUPEX.plot(Solution_coupeX, Elem_ds_coupeX[:, 1], label="Solution Numérique")
        COUPEX.plot(SolutionEX_coupeX, Elem_ds_coupeX[:, 1], '--', label="Solution MMS")
        COUPEX.set_xlabel("Température")
        COUPEX.set_ylabel("Y (m)")
        COUPEX.set_title(f"Solution dans une coupe à X = {X_Coupe}")
        COUPEX.legend()

        COUPEY.plot(Elem_ds_coupeY[:, 0], Solution_coupeY, label="Solution Numérique")
        COUPEY.plot(Elem_ds_coupeY[:, 0], SolutionEX_coupeY, '--', label="Solution MMS/analytique")
        COUPEY.set_xlabel("X (m)")
        COUPEY.set_ylabel("Température")
        COUPEY.set_title(f"Solution dans une coupe à Y = {Y_Coupe}")
        COUPEY.legend()

        # Enregistrer
        plt.savefig(save_path, dpi=200)

    def show_mesh_differences(self, mesh1, mesh2, title, save_path, diff=False):
        """
        Affichage des graphiques qui montrent les résultats entre deux types
        de maillage

        Parameters
        ----------
        mesh1: mesh
            Maillage 1 de l'exemple traité. 
        
        mesh2: mesh
            Maillage 2 de l'exemple traité. 
        
        title: str
            Nom du document pour le sauvegarder
            
        save_path: str
            Nom du fichier de sauvegarde
            
        diff: Bool
            Pour décider si on trace la différence (Erreur)  entre les deux maillages. 
        Returns
        -------
        None

        """
        if diff is True:
            figure, (plot1, plot2, plot3) = plt.subplots(1, 3, figsize=(28, 6))
        else:
            figure, (plot1, plot2) = plt.subplots(1, 2, figsize=(20, 6))

        figure.suptitle(title)
        # Set levels of color for the colorbar
        levels = np.linspace(np.min(np.append(self.data[mesh1[0]]['phi_num'], self.data[mesh2[0]]['phi_num'])),
                             np.max(np.append(self.data[mesh1[0]]['phi_num'], self.data[mesh2[0]]['phi_num'])), num=40)

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

        if diff is True:
            err = np.abs(self.data[mesh1[0]]['phi_num'] - self.data[mesh2[0]]['phi_num'])

            levels = np.linspace(np.min(err), np.max(err), num=40)

            c = plot3.tricontourf(center1[:, 0], center1[:, 1], err, levels=levels)
            plot3.set_xlabel("L (m)")
            plot3.set_ylabel("H (m)")
            plot3.set_title(f"Erreur absolue entre les maillages à {self.data[mesh1[0]]['n']} éléments")
            plt.colorbar(c, ax=plot3)


        # Enregistrer
        plt.savefig(save_path, dpi=200)

    def show_time(self, title, save_path):
        """
        Affichage des graphiques qui montre le temps de calcul entre les
        deux méthodes Sparse et Dense

        Parameters
        ----------
        
        title: str
            Nom du document pour le sauvegarder
            
        save_path: str
            Nom du fichier de sauvegarde
            
        Returns
        -------
        None

        """
        Figure1, (Comp) = plt.subplots(figsize=(15, 10))

        Figure1.suptitle(title)

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
        """
        Affichage des graphiques d'ordre de convergence et calcul de l'erreur
        par rapport au solution exacte. 

        Parameters
        ----------
        None
            
        Returns
        -------
        None

        """
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