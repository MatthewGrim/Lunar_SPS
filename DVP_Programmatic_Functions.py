""""
11/06/2018
Author: Darian van Paridon

This script contains function used in the generation of data via parametric scanning in STK, and for the post processing
and visualization of the data.

"""


import math
import os
import numpy as np
import time
import matplotlib.pyplot as plt


def vary_orbital_elements(resolution, min_perigee, max_perigee, max_apogee):

    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK

    radius_moon = 1737.0
    orbit_data = [0.0, 0.0]

    print('Generating orbital elements...')

    start_time = time.time()
    for j in range(1, int((max_perigee - min_perigee) / resolution) + 1):
        perigee = min_perigee + (j * resolution)
        apogee = perigee
        while apogee <= max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            apogee += resolution

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])
    end_time = time.time()

    print('Time required to generate {} sets of orbital elements: {} seconds'.format(len(orbit_data), round(end_time - start_time, 5)))

    return semi_maj_axis, eccentricity, orbit_data


def sort_data_list_into_array(orbit_data, resolution, data_list):

    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    # Initialize arrays for plotting
    data_array = np.zeros((len(unique_perigees), len(unique_apogees)))

    # Sort data into a arrays for 2D plot
    max_apogee = max(unique_apogees)
    num_apogees = list()
    for k in range(len(unique_perigees)):
        # Determine how many apogee steps required to reach the max value.
        num_apogees.append(int(math.floor((max_apogee - unique_perigees[k]) / resolution)))
        # For perigee value (index k), collect data points for corresponding apogees (index j)
        # where j starts at the first unique apogee equal to the current perigee
        j = 0
        while unique_apogees[j] < unique_perigees[k]:
            j += 1
        start = j
        for j in range(0, num_apogees[k] + 1):
            if k == 0:
                idx = j + start
            else:
                idx = int(np.sum(num_apogees) - num_apogees[k] + j + start)
            data_array[k, j + start] = data_list[idx]

    return data_array


def vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee):

    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK

    radius_moon = 1737.0

    print('Generating orbital elements...')

    start_time = time.time()
    resolution = np.array((10.0, 25.0, 50.0, 100.0))
    thresholds = np.array((250.0, 500.0, 1000.0))
    orbit_data = [0, 0]
    perigee = 0.0
    while perigee <= max_perigee:
        if 0.0 <= perigee < thresholds[0]:
            peri_step = resolution[0]
        elif thresholds[0] <= perigee < thresholds[1]:
            peri_step = resolution[1]
        elif thresholds[1] <= perigee < thresholds[2]:
            peri_step = resolution[2]
        elif thresholds[2] <= perigee:
            peri_step = resolution[3]
        perigee = perigee + peri_step
        apogee = perigee
        while apogee <= max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            if 0.0 <= apogee < thresholds[0]:
                apo_step = resolution[0]
            elif thresholds[0] <= apogee < thresholds[1]:
                apo_step = resolution[1]
            elif thresholds[1] <= apogee < thresholds[2]:
                apo_step = resolution[2]
            elif thresholds[2] <= apogee:
                apo_step = resolution[3]
            apogee += apo_step

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])
    end_time = time.time()

    print('Time required to generate {} sets of orbital elements: {} seconds'.format(len(orbit_data), round(end_time - start_time, 5)))

    return semi_maj_axis, eccentricity, orbit_data


def sort_data_list_with_incremented_resolution_into_array(orbit_data, data_list):

    r_moon = 1737.0
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])
    max_apogee = max(unique_apogees)

    # Initialize arrays for plotting
    data_array = np.zeros((len(unique_perigees), len(unique_apogees)))

    resolution = np.array((10.0, 25.0, 50.0, 100.0))
    thresholds = np.array((250.0, 500.0, 1000.0))
    steps_per_interval = np.array((thresholds[0] / resolution[0],
                                   (thresholds[1] - thresholds[0]) / resolution[1],
                                   (thresholds[2] - thresholds[1]) / resolution[2],
                                   (max_apogee - r_moon - thresholds[2]) / resolution[3]))

    thresholds = np.array((250.0 + r_moon, 500.0 + r_moon, 1000.0 + r_moon))

    # Sort data into a arrays for 2D plot
    num_apogees = []
    for k in range(len(unique_perigees)):
        # Determine how many apogee steps required to reach the max value from current perigee
        if r_moon <= unique_perigees[k] < thresholds[0]:
            num_steps = int((thresholds[0] - unique_perigees[k]) / resolution[0]) + int(np.sum(steps_per_interval[1:]))
        elif thresholds[0] <= unique_perigees[k] < thresholds[1]:
            num_steps = int((thresholds[1] - unique_perigees[k]) / resolution[1]) + int(np.sum(steps_per_interval[2:]))
        elif thresholds[1] <= unique_perigees[k] < thresholds[2]:
            num_steps = int((thresholds[2] - unique_perigees[k]) / resolution[2]) + int(np.sum(steps_per_interval[3:]))
        elif thresholds[2] <= unique_perigees[k]:
            num_steps = int((max_apogee - unique_perigees[k]) / resolution[3])
        # For perigee value (index k), collect data points for corresponding apogees (index j)
        num_apogees.append(num_steps + 1)
        # where j starts at the first unique apogee equal to the current perigee (circular orbit)
        j = 0
        while unique_apogees[j] < unique_perigees[k]:
            j += 1
        start = j
        for j in range(0, num_apogees[k]):
            idx = j + start
            data_array[k, j + start] = data_list[idx]

    return data_array


def write_data_to_file(stk_data_path, study_name, data, data_name):

    try:
        os.remove('{}/{}_{}.txt'.format(stk_data_path, data_name, study_name))
    except OSError:
        pass
    data_array = np.asarray(data)
    np.savetxt('{}/{}_{}.txt'.format(stk_data_path, data_name, study_name), data_array)


def read_data_from_file(stk_data_path, study_name, data_name):

    # Read in saved active time data, if the data is already processed and saved into a txt file
    # Comment this section out if the data is being processed, and the results are not yet written to a txt file
    data = []
    fh = open('{}/{}_{}.txt'.format(stk_data_path, data_name, study_name))
    for i, line in enumerate(fh):
        temp = line.split("\n")
        data.append(float(temp[0]))

    return data


def make_contour_plot(X, Y, data, title, fig_num):
    plt.figure(fig_num)
    plt.contourf(Y, X, data, 500)
    plt.colorbar()
    plt.title('{}'.format(title))
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.show()
