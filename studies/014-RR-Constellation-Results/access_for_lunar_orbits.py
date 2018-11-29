"""
Author: Rohan Ramasamy
Date: 26/11/2018

This script considers access times from a frozen orbit at a 1000km periselene and 10000km aposelene.
The analysis is similar to study 002/Access_and_Blackout_Durations.py
"""

import numpy as np
import os
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def read_element_data(fname, start):
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

def plot_classical_elements():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    main_dir = "LunarComms Orbit"
    sub_dir = "10000"

    name = "classical_elements_lunarComms.csv"
    path = os.path.join(main_dir, sub_dir, name)
    assert os.path.exists(path), path
    data = read_element_data(path, start)
    
    fig, ax = plt.subplots(3, figsize=(15, 8))
    data[0, :] /= 86400
    ax[0].plot(data[0, :], data[1, :])

    ax[1].plot(data[0, :], data[2, :])
    
    ax[2].plot(data[0, :], data[3, :])
    xlim = data[0, -1]

    ax[0].set_xlim([0.0, xlim])
    ax[1].set_xlim([0.0, xlim])
    ax[2].set_xlim([0.0, xlim])
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    
    ax[0].set_ylabel("Semi Major Axis [$km$]")
    ax[1].set_ylabel("Eccentricity")
    ax[2].set_ylabel("Inclination [$^\circ$]")
    name = "optimum_equatorial_orbits"    
    plt.tight_layout()
    plt.savefig(os.path.join(main_dir, sub_dir, "frozen_orbit_elements"))
    plt.close()


def get_frozen_orbit_access():
    main_dir = "LunarComms Orbit"
    sub_dir = "10000"
    sps_lighting_file = "Lighting_Times_lunarComms.csv"
    sps_access_file = "Access_Modified_LunarComms.csv"
    target_lighting = "Target_Lighting.csv"
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    sps_access = parse_csv_to_array(os.path.join(main_dir, sub_dir, sps_access_file), start)
    sps_lighting = parse_csv_to_array(os.path.join(main_dir, sub_dir, sps_lighting_file), start)
    target_lighting = parse_csv_to_array(os.path.join(main_dir, target_lighting), start)
    target_eclipse = invert_events_list(target_lighting, total_duration)
    
    sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
    sps_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)

    sps_active_days = [i / 86400.0 for i in sps_active[0]]
    sps_active_dur = [i / 3600.0 for i in sps_active[2]]   
    sps_blackout_days = [i / 86400.0 for i in sps_blackout[0]]
    sps_blackout_dur = [i / 3600.0 for i in sps_blackout[2]]   
    
    fig, ax = plt.subplots(2, figsize=(12, 7))

    ax[0].bar(sps_active_days, sps_active_dur)
    ax[0].set_ylabel("Active duration [hr]")
    ax[0].set_xlabel("Time [days]")

    ax[1].bar(sps_blackout_days, sps_blackout_dur)
    ax[1].set_ylabel("Blackout duration [hr]")
    ax[1].set_xlabel("Time [days]")

    plt.savefig(os.path.join(main_dir, sub_dir, "frozen_orbit_access"))
    plt.close()

if __name__ == '__main__':
    get_frozen_orbit_access()
    plot_classical_elements()

