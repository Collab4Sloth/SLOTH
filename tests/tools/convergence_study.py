import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import csv
from argparse import ArgumentParser, RawTextHelpFormatter


def find_and_read_csv(base_dir, val):
    orders_dict = {}

    for root, dirs, files in os.walk(base_dir):                                 # Loop through folders
        if 'Problem1' in root and 'time_specialized.csv' in files:              
            base_folder = os.path.basename(os.path.dirname(root))               # Extract informations : FE order and number of elements
            parts = base_folder.split('_')
            if len(parts) > 2 and parts[0] == 'Saves' and parts[1] == 'order':  # Check if folder name is Saves_order_XX_Nx_YY_*
                try:
                    order = int(parts[2])                                       # Extract order
                except ValueError:
                    continue  
                nx_part = parts[3]                                              # Extract Nx and convert it to int
                try:
                    nx_value = int(nx_part.replace('Nx', ''))
                except ValueError:
                    continue 

                csv_path = os.path.join(root, 'time_specialized.csv')           # Check if 'time_specialized.csv' exist
                if os.path.isfile(csv_path):
                    df = pd.read_csv(csv_path)                                  # Read the CSV
                    if df.shape[1] >= 4:                                        # Check the validity of the file
                        fourth_column = df.iloc[:, 3].tolist()                  # Extract 4th col : L2 error value
                        third_column = df.iloc[:, 2].tolist()                   # Extract 3th col : time value
                        try:
                            index_val = third_column.index(val)                 # Find index of time value in file
                        except ValueError:  
                            index_val = -1                                      # Last time step if value don't exist                                
                            print("Convergence done for the last timestep")

                        pi_over_n = LL / nx_value

                        data = []
                        data.append([pi_over_n, fourth_column[index_val]])      # create list of list : [elements size, L2_error at the wanted time]

                        key = f'order {order}'
                        if key not in orders_dict:
                            orders_dict[key] = []

                        orders_dict[key].extend(data)                           # add list to a dictionnary for the specified order

    sorted_orders_dict = {                                                      # sort data
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
            try:
                x_fit = np.linspace(min(x_values), max(x_values), 100)
                slope, intercept, r_value, p_value, std_err = linregress(
                    np.log(x_values), np.log(y_values))                         # fit curve to find convergence order
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
        plt.close()


if __name__ == "__main__":
    parser = ArgumentParser(description='Export LogLog plot of L2error = f(element size)',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-l", "--domainSize", type=float,
                         default=1., help="Domain Size")
    parser.add_argument("-t", "--timeStep", type=float,
                         default=-1, help="Time Step")

    args = parser.parse_args()
    LL = args.domainSize
    val = args.timeStep
    base_directory = './'
    result = find_and_read_csv(base_directory,val)

    plot_data(result)

