"""
Author: Rohan Ramasamy
Date: 12/11/2018

This script assesses the variation in azimuth, elevation and range for the optimised orbit from the IAC.
"""

import os
import numpy as np
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import convert_string_to_datetime


def read_data(fname, start, is_sat=True):
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
        if is_sat:
	        data[3, i - 1] = float(components[3])
	        data[4, i - 1] = float(components[4])
	        data[5, i - 1] = float(components[5])
	        data[6, i - 1] = float(components[6])

    return data
        

def plot_aer_sat_data():
    file = "aer_report.csv"

    fig, ax = plt.subplots(2, 3, figsize=(15, 6), sharex='col')
    path = os.path.join(main_dir, file)
    assert os.path.exists(path), path

    data = read_data(path, start)
    data[0, :] /= 86400
    ax[0, 0].scatter(data[0, :], data[1, :], label='azimuth')
    ax[0, 1].scatter(data[0, :], data[2, :], label='elevation')
    ax[0, 2].scatter(data[0, :], data[3, :], label='range')
    ax[1, 0].scatter(data[0, :], data[4, :], label='azimuth rate')
    ax[1, 1].scatter(data[0, :], data[5, :], label='elevation rate')
    ax[1, 2].scatter(data[0, :], data[6, :], label='range rate')

    ax[0, 0].set_xlim([0.0, 1.5])
    ax[0, 0].set_xlim([0.0, 1.5])
    ax[0, 1].set_xlim([0.0, 1.5])
    ax[0, 2].set_xlim([0.0, 1.5])
    ax[1, 0].set_xlim([0.0, 1.5])
    ax[1, 1].set_xlim([0.0, 1.5])
    ax[1, 2].set_xlim([0.0, 1.5])
    ax[0, 0].set_ylim([np.min(data[1, :]), np.max(data[1, :])])
    ax[0, 1].set_ylim([np.min(data[2, :]), np.max(data[2, :])])
    ax[0, 2].set_ylim([np.min(data[3, :]), np.max(data[3, :])])
    ax[1, 0].set_ylim([np.min(data[4, :]), np.max(data[4, :])])
    ax[1, 1].set_ylim([np.min(data[5, :]), np.max(data[5, :])])
    ax[1, 2].set_ylim([np.min(data[6, :]), np.max(data[6, :])])
    
    ax[0, 0].set_ylabel("Azimuth [$^\circ$]")
    ax[0, 1].set_ylabel("Elevation [$^\circ$]")
    ax[0, 2].set_ylabel("Range [$km$]")
    ax[1, 0].set_ylabel("Azimuth [$^\circ s^{-1}$]")
    ax[1, 1].set_ylabel("Elevation [$^\circ s^{-1}$]")
    ax[1, 2].set_ylabel("Range [$kms^{-1}$]")
    ax[1, 0].set_xlabel("Days]")
    ax[1, 1].set_xlabel("Days]")
    ax[1, 2].set_xlabel("Days]")

    plt.tight_layout()
    plt.savefig(os.path.join(main_dir, "AER_satellite_results"))
    plt.show()
    plt.close()

def plot_aer_target():
    file = "target_ae_report.csv"

    fig, ax = plt.subplots(2, figsize=(15, 6), sharex='col')
    path = os.path.join(main_dir, file)
    assert os.path.exists(path), path

    data = read_data(path, start, is_sat=False)
    data[0, :] /= 86400
    ax[0].scatter(data[0, :], data[1, :], label='azimuth')
    ax[1].scatter(data[0, :], data[2, :], label='elevation')

    ax[0].set_xlim([0.0, 1.5])
    ax[1].set_xlim([0.0, 1.5])
    ax[0].set_ylim([np.min(data[1, :]), np.max(data[1, :])])
    ax[1].set_ylim([np.min(data[2, :]), np.max(data[2, :])])
    
    ax[0].set_ylabel("Azimuth [$^\circ$]")
    ax[1].set_ylabel("Elevation [$^\circ$]")

    plt.tight_layout()
    plt.savefig(os.path.join(main_dir, "AER_target_results"))
    plt.show()
    plt.close()


if __name__ == '__main__':
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    main_dir = "AER_Validation"
    plot_aer_sat_data()
    plot_aer_target()