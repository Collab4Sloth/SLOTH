import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import csv
from argparse import ArgumentParser, RawTextHelpFormatter


def find_and_read_csv(base_dir):
    # Créer un dictionnaire pour stocker les données organisées par ordre
    orders_dict = {}

    # Parcourir tous les dossiers dans le répertoire de base
    for root, dirs, files in os.walk(base_dir):
        # Vérifier si le dossier contient un fichier 'time_specialized.csv'
        if 'Problem1' in root and 'time_specialized.csv' in files:
            # Extraire la partie order et Nx du chemin du dossier
            base_folder = os.path.basename(os.path.dirname(root))
            parts = base_folder.split('_')

            # Assumer que le dossier est de la forme 'Saves_order_<order>_Nx<value>'
            if len(parts) > 2 and parts[0] == 'Saves' and parts[1] == 'order':
                try:
                    order = int(parts[2])  # Extraire l'ordre
                except ValueError:
                    continue  # Ignorer si la partie n'est pas un entier valide

                # Extraire la valeur de Nx
                nx_part = parts[-1]  # 'Nx...'
                try:
                    nx_value = int(nx_part.replace('Nx', ''))
                except ValueError:
                    continue  # Ignorer si la partie n'est pas un entier valide

                # Définir le chemin vers 'time_specialized.csv' dans 'Problem1'
                csv_path = os.path.join(root, 'time_specialized.csv')
                if os.path.isfile(csv_path):
                    # Lire le fichier CSV
                    df = pd.read_csv(csv_path)
                    if df.shape[1] >= 4:  # Vérifier s'il y a au moins 4 colonnes
                        # Extraire la 4ème colonne
                        fourth_column = df.iloc[:, 3].tolist()

                        # Calculer pi/n
                        pi_over_n = LL / nx_value

                        # Préparer les données sous forme de tableau
                        data = []
                        # TODO ici on trace la dernière ligne du csv, faire en sorte que l'utilistaeur puisse choisir
                        data.append([pi_over_n, fourth_column[-1]])

                        # Ajouter les données au dictionnaire sous la clé 'order <ordre>'
                        key = f'order {order}'
                        if key not in orders_dict:
                            orders_dict[key] = []

                        # Ajouter les données au tableau
                        orders_dict[key].extend(data)

    # Trier les données pour chaque ordre par la première colonne (pi/n)
    sorted_orders_dict = {
        key: sorted(orders_dict[key], key=lambda x: x[0])
        for key in sorted(orders_dict, reverse=True)
    }

    return sorted_orders_dict


def plot_data(orders_dict):
    color = ["red", "blue", "black", "brown", "green", "orange"]
    nom_fichier = "convergence_output.csv"
    with open(nom_fichier, mode='w', newline='', encoding='utf-8') as fichier_csv:
        writer = csv.writer(fichier_csv)
        i = 0
        for order, data in orders_dict.items():
            x_values = [row[0] for row in data]
            y_values = [row[1] for row in data]
            for j, x in enumerate(x_values):
                writer.writerow([order, x, y_values[j]])
            print(color[i])
            plt.loglog(x_values, y_values, marker='.', markersize=12,
                       linestyle="none", color=color[i], label=f'{order}')
            # Ajustement du modèle
            try:
                x_fit = np.linspace(min(x_values), max(x_values), 100)
                slope, intercept, r_value, p_value, std_err = linregress(
                    np.log(x_values), np.log(y_values))
                print("Order " + str(order) + ", exponant coefficient = " +
                      str(slope) + ", proportional coefficient = " + str(intercept))
                plt.plot(x_fit, np.exp(intercept)*x_fit**(slope),
                         color=color[i], linestyle='--', label=f'{order} fit: $\\propto h^{{{slope:.4f}}}$')
            except Exception as e:
                print(f"Error fitting curve for {order}: {e}")
            i += 1
        plt.title(r'Convergence study')
        plt.xlabel(r'Element size h [m]')
        plt.ylabel('$L_2$ error [-]')
        plt.legend()
        plt.grid(True, which="both", ls="--")
        plt.savefig('convergence_output.png')
        # plt.show()
        plt.close()


if __name__ == "__main__":
    parser = ArgumentParser(description='Export LogLog plot of L2error = f(element size)',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-l", "--domainSize", type=float,
                         default=1., help="Domain Size")
    args = parser.parse_args()
    LL = args.domainSize
    # Remplacer './' par le chemin vers votre répertoire de base
    base_directory = './'
    result = find_and_read_csv(base_directory)

    # Tracer les données
    plot_data(result)
