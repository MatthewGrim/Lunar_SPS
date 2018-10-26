"""
Author: Rohan Ramasamy
Date: 26/10/2018

Script to compare results of orbital elements for J4 and HPOP
"""

import os
import numpy as np
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import convert_string_to_datetime


def read_data(fname):
    f = open(fname)
    lines = f.readlines()
    data = np.zeros((8, len(lines) - 2))
    for i, line in enumerate(lines[:-1]):
        if i == 0:
            continue

        components = line.split(",")

        # Time
        data[0, i - 1] = i
        
        # Orbital elements - sma, e, i, raan, peri, true anomaly, mean anomaly
        data[1, i - 1] = float(components[1])
        data[2, i - 1] = float(components[2])
        data[3, i - 1] = float(components[3])
        data[4, i - 1] = float(components[4])
        data[5, i - 1] = float(components[5])
        data[6, i - 1] = float(components[6])
        data[7, i - 1] = float(components[7])

    return data
        

def compare_propagators():
    main_dir = "Validation"
    orbit_dirs = ["Equatorial", "Polar"]
    prop_dirs = ["HPOP", "J4"]
    peri_apo_pairs = [(800, 800), (5000, 5000), (800, 5000)]

    for k, pair in enumerate(peri_apo_pairs):
        for i, orbit in enumerate(orbit_dirs):
            fig, ax = plt.subplots(2, 3, figsize=(10, 10))
            for j, prop in enumerate(prop_dirs):
                    name = "element_variation_{}_{}_coarse.csv".format(pair[1], pair[0])
                    path = os.path.join(main_dir, orbit, prop, name)
                    assert os.path.exists(path), path

                    data = read_data(path)
                    ax[0, 0].plot(data[0, :], data[1, :], label=prop)

                    ax[0, 1].plot(data[0, :], data[2, :], label=prop)
                    
                    ax[0, 2].plot(data[0, :], data[3, :], label=prop)
                    
                    ax[1, 0].plot(data[0, :], data[4, :], label=prop)
                    
                    ax[1, 1].plot(data[0, :], data[5, :], label=prop)
                    
                    ax[1, 2].plot(data[0, :], data[6, :], label=prop)
                    ax[1, 2].plot(data[0, :], data[7, :], label=prop)

            ax[0, 0].set_title("Semi Major Axis")
            ax[0, 1].set_title("Eccentricity")
            ax[0, 2].set_title("Inclination")
            ax[1, 0].set_title("RAAN")
            ax[1, 1].set_title("Argument of Perigee")
            ax[1, 2].set_title("True/Mean Anomaly")
            name = "{}: {}km by {}km".format(orbit, pair[0], pair[1])
            fig.suptitle(name)       
            plt.savefig(name)
            plt.close()


if __name__ == '__main__':
    compare_propagators()

