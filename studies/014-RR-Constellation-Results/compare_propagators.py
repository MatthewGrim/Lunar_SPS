"""
Author: Rohan Ramasamy
Date: 26/10/2018

Script to compare results of orbital elements for J4 and HPOP
"""

import os
import numpy as np
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import convert_string_to_datetime


def read_data(fname, start):
    f = open(fname)
    lines = f.readlines()
    data = np.zeros((8, len(lines) - 2))
    for i, line in enumerate(lines[:-1]):
        if i == 0:
            continue

        components = line.split(",")

        # Time
        data[0, i - 1] = (convert_string_to_datetime(components[0].split(':')) - start).total_seconds()
        
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
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    main_dir = "Validation"
    orbit_dirs = ["Equatorial", "Polar"]
    prop_dirs = ["HPOP", "J4"]
    peri_apo_pairs = [(800, 800), (5000, 5000), (800, 5000)]

    for k, pair in enumerate(peri_apo_pairs):
        for i, orbit in enumerate(orbit_dirs):
            fig, ax = plt.subplots(3, figsize=(15, 8))
            xlim = None
            for j, prop in enumerate(prop_dirs):
                    name = "element_variation_{}_{}_coarse.csv".format(pair[1], pair[0])
                    path = os.path.join(main_dir, orbit, prop, name)
                    assert os.path.exists(path), path

                    data = read_data(path, start)
                    data[0, :] /= 86400
                    ax[0].plot(data[0, :], data[1, :], label=prop)

                    ax[1].plot(data[0, :], data[2, :], label=prop)
                    
                    ax[2].plot(data[0, :], data[3, :], label=prop)
                    xlim = data[0, -1] if j == 0 else min(data[0, -1], xlim)


                    # ax[1, 0].plot(data[0, :], data[4, :], label=prop)
                    
                    # ax[1, 1].plot(data[0, :], data[5, :], label=prop)
                    
                    # ax[1, 2].plot(data[0, :], data[6, :], label=prop)
                    # ax[1, 2].plot(data[0, :], data[7, :], label=prop)

            ax[0].set_xlim([0.0, xlim])
            ax[1].set_xlim([0.0, xlim])
            ax[2].set_xlim([0.0, xlim])
            
            ax[0].set_ylabel("Semi Major Axis [$km$]")
            ax[1].set_ylabel("Eccentricity")
            ax[2].set_ylabel("Inclination [$^\circ$]")
            name = "{}: {}km by {}km".format(orbit, pair[0], pair[1])    
            plt.tight_layout()
            plt.savefig(os.path.join(main_dir, name))
            plt.close()


if __name__ == '__main__':
    compare_propagators()

