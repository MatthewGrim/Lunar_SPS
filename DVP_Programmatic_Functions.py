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

    for j in range(0, int((max_perigee - min_perigee) / resolution) + 1):
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
        num_apogees.append(int((max_apogee - unique_perigees[k]) / resolution))
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

    data_array[data_array == 0.0] = np.nan

    return data_array


def vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee, min_perigee=None,
                                                  resolutions=np.array((10.0, 25.0, 50.0, 100.0)), thresholds=np.array((100.0, 250.0, 1000.0))):
    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK
    if min_perigee is not None:
        assert min_perigee < max_perigee
    assert resolutions.shape == (4,)
    assert thresholds.shape == (3,)

    radius_moon = 1737.0
    orbit_data = [0, 0]
    if min_perigee is not None:
        perigee = min_perigee
        apogee = perigee    
        orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
    else:
        perigee = 0.0
    
    while perigee <= max_perigee:
        # Set perigee step
        if 0.0 <= perigee < thresholds[0]:
            peri_step = resolutions[0]
        elif thresholds[0] <= perigee < thresholds[1]:
            peri_step = resolutions[1]
        elif thresholds[1] <= perigee < thresholds[2]:
            peri_step = resolutions[2]
        elif thresholds[2] <= perigee:
            peri_step = resolutions[3]
                
        # Add apogees for current perigee
        perigee = perigee + peri_step
        apogee = perigee
        while apogee <= max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            if 0.0 <= apogee < thresholds[0]:
                apo_step = resolutions[0]
            elif thresholds[0] <= apogee < thresholds[1]:
                apo_step = resolutions[1]
            elif thresholds[1] <= apogee < thresholds[2]:
                apo_step = resolutions[2]
            elif thresholds[2] <= apogee:
                apo_step = resolutions[3]
            apogee += apo_step

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])

    return semi_maj_axis, eccentricity, orbit_data


def sort_incremented_resolution_data(orbit_data, data_list):

    # Get list of unique perigee and apogee altitudes
    r_moon = 1737.0
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])
    max_apogee = max(unique_apogees)

    data_array = np.zeros([len(unique_perigees), len(unique_apogees)])

    # Define resolution and bounds for incremented resolution approach
    resolution = [10.0, 25.0, 50.0, 100.0]
    thresholds = [100.0 + r_moon, 250.0 + r_moon, 1000.0 + r_moon]
    num_apogees = np.zeros([len(unique_perigees)])

    # Determine number of unique apogee steps between a given perigee value and max apogee value
    for j in range(len(unique_perigees)):
        if r_moon < unique_perigees[j] < thresholds[0]:
            num_apogees[j] = 1 + int((thresholds[0] - unique_perigees[j]) / resolution[0]) + int((thresholds[1] - thresholds[0]) / resolution[1]) + int((thresholds[2] - thresholds[1]) / resolution[2]) + int((max_apogee - thresholds[2]) / resolution[3])
        elif thresholds[0] <= unique_perigees[j] < thresholds[1]:
            num_apogees[j] = 1 + int((thresholds[1] - unique_perigees[j]) / resolution[1]) + int((thresholds[2] - thresholds[1]) / resolution[2]) + int((max_apogee - thresholds[2]) / resolution[3])
        elif thresholds[1] <= unique_perigees[j] < thresholds[2]:
            num_apogees[j] = 1 + int((thresholds[2] - unique_perigees[j]) / resolution[2]) + int((max_apogee - thresholds[2]) / resolution[3])
        elif thresholds[2] <= unique_perigees[j]:
            num_apogees[j] = 1 + int((max_apogee - unique_perigees[j]) / resolution[3])

    # Sort data
    start = np.zeros([len(unique_perigees)])
    for j in range(len(unique_perigees)):
        start[j] = np.sum(num_apogees[0:j])
        for k in range(int(num_apogees[j])):
            shift = int(len(unique_apogees) - num_apogees[j])
            data_array[j, shift + k] = data_list[k + int(start[j])]

    # Remove zeros for impossible orbit combos (with perigee altitude > apogee altitude)
    for i in range(len(data_array[0])):
        for j in range(len(data_array[1])):
            if j < i:
                data_array[i][j] = np.nan

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
    plt.contourf(X, Y, data, 500)
    plt.colorbar()
    plt.title('{}'.format(title))
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.show()


def get_event_statistics(X, Y, data, rover, title):
    fig, ax = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(15, 8))

    num_contours = 100
    im = ax[0, 0].contourf(X, Y, data["min_active_time"], num_contours)
    fig.colorbar(im, ax=ax[0, 0])
    im = ax[0, 1].contourf(X, Y, data["max_active_time"], num_contours)
    fig.colorbar(im, ax=ax[0, 1])
    im = ax[0, 2].contourf(X, Y, data["mean_active_time"], num_contours)
    fig.colorbar(im, ax=ax[0, 2])
    im = ax[0, 3].contourf(X, Y, data["std_active_time"], num_contours)
    fig.colorbar(im, ax=ax[0, 3])

    im = ax[1, 0].contourf(X, Y, data["min_blackout_time"], num_contours)
    fig.colorbar(im, ax=ax[1, 0])
    ax[1, 0].set_ylabel("Apogee Altitude [km]")
    im = ax[1, 1].contourf(X, Y, data["max_blackout_time"], num_contours)
    fig.colorbar(im, ax=ax[1, 1])
    ax[1, 1].set_ylabel("Apogee Altitude [km]")
    im = ax[1, 2].contourf(X, Y, data["mean_blackout_time"], num_contours)
    fig.colorbar(im, ax=ax[1, 2])
    ax[1, 2].set_ylabel("Apogee Altitude [km]")
    im = ax[1, 3].contourf(X, Y, data["std_blackout_time"], num_contours)
    fig.colorbar(im, ax=ax[1, 3])
    ax[1, 3].set_ylabel("Apogee Altitude [km]")

    ax[0, 0].set_title("Min Duration [hrs]")
    ax[0, 1].set_title("Max Duration [hrs]")
    ax[0, 2].set_title("Mean Duration [hrs]")
    ax[0, 3].set_title("Standard Deviation in Duration [hrs]")

    ax[0, 0].set_ylabel("Active Events\n Perigee Altitude [km]")
    ax[1, 0].set_ylabel("Blackout Events\n Perigee Altitude [km]")
    fig.suptitle(title)
    fig.tight_layout()
    plt.show()

    plt.figure()

    power_balance = data["mean_active_time"] / data["mean_blackout_time"] * (rover["operation_pwr"] - rover["hibernation_pwr"]) / rover["hibernation_pwr"]
    plt.contourf(X, Y, power_balance, num_contours)
    plt.colorbar()

    plt.xlabel("Apogee Altitude [km]")
    plt.ylabel("Perigee Altitude [km]")
    plt.title("Ratio of Power Available to Power Consumed")
    plt.show()


def calculate_orbital_perturbations(semi_maj_axis, eccentricity, inclination_ep, arg_perigee_ep):
    import sympy
    from sympy import cos, sin

    semi_maj_axis = [i * 1000.0 for i in semi_maj_axis]
    inclination_ep = inclination_ep * (np.pi / 180.0)
    arg_perigee_ep = arg_perigee_ep * (np.pi / 180.0)

    # Initialize solar system data
    mass_moon = 7.34767309e22
    r_moon = 1737.63e3

    # Zonal harmonics of moon from GLGM-1 gravitational model
    J2 = 0.00020374485
    J3 = 0.000008476
    mass_earth = 5.972e24
    mass_fraction = mass_earth / (mass_earth + mass_moon)
    seperation_earth_moon = 385000e3
    G = 6.67384e-11

    # Intermediate variables to simplify equations
    Ra = [r_moon / i for i in semi_maj_axis]
    e = [1.0 - j ** 2 for j in eccentricity]

    # Apparent mean motion of Earth around moon
    mean_motion_earth = np.sqrt(G * (mass_earth + mass_moon) / seperation_earth_moon ** 3)

    # Transform into frame used in Ely paper (apparent orbital plane of Earth around moon)
    relative_inclination = 6.8 * np.pi / 180.0
    RAAN_op = 0.0
    arg_perigee_op = arg_perigee_ep
    i = sympy.Symbol('i')
    inclination_op = sympy.solve(cos(relative_inclination) * cos(i) - sin(relative_inclination) * cos(RAAN_op) * sin(i), i)
    # Select positive solution
    inclination_op = np.asarray(inclination_op)
    inclination_op = inclination_op[inclination_op > 0]

    # initialize arrays
    dwdt_earth = np.zeros(len(semi_maj_axis))
    dedt_earth = np.zeros(len(semi_maj_axis))
    dedt_oblate = np.zeros(len(semi_maj_axis))
    dwdt_oblate = np.zeros(len(semi_maj_axis))
    didt_oblate = np.zeros(len(semi_maj_axis))

    # Calculate perturbations on initial orbits
    for j in range(len(eccentricity)):
        # Due to Earth's gravity
        dwdt_earth[j] = 3.0 * (5.0 * math.cos(inclination_op) ** 2 - 1.0 + eccentricity[j] ** 2 + 5.0 * (1 - eccentricity[j] ** 2 - math.cos(inclination_op) ** 2) * math.cos(2 * arg_perigee_op)) * mass_fraction * mean_motion_earth ** 2 / (8 * math.sqrt(1 - eccentricity[j] ** 2) * math.sqrt(G * mass_moon / semi_maj_axis[j] ** 3))
        dedt_earth[j] = 15.0 * mass_fraction * mean_motion_earth ** 2 * eccentricity[j] * np.sqrt(1 - eccentricity[j] ** 2) * math.sin(inclination_op) ** 2 * math.sin(2 * arg_perigee_op) / (8 * math.sqrt(G * mass_moon / semi_maj_axis[j] ** 3))

        # Due to oblateness of moon
        mean_motion_sat = np.sqrt(G * mass_moon / semi_maj_axis[j] ** 3)

        didt_oblate[j] = - 1.5 * Ra[j] ** 3 * J3 * eccentricity[j] * np.cos(inclination_ep) * np.cos(arg_perigee_ep) * (1.25 * np.sin(inclination_ep) ** 2 - 1) / (e[j] ** 2)
        dedt_oblate[j] = 1.5 * mean_motion_sat * Ra[j] ** 3 * J3 * np.sin(inclination_ep) * np.cos(arg_perigee_ep) * (np.sin(inclination_ep) ** 2 - 1) / (e[j] ** 2)

        # Manage singularity for circular orbits
        if eccentricity[j] == 0.0:
            # Non divergent term. Can be used, although for circular orbits it is not practically useful.
            # dwdt_oblate[j] = -0.75 * mean_motion_sat * Ra[j] ** 2 * J2 * (1 - 5 * np.cos(inclination_ep) ** 2)
            dwdt_oblate[j] = 0.0
        elif inclination_ep == 0.0:
            dwdt_oblate[j] = 0.0
        else:
            dwdt_oblate[j] = - 0.75 * mean_motion_sat * Ra[j] ** 2 * J2 * (1 - 5.0 * np.cos(inclination_ep) ** 2) / e[j] ** 2 \
            - 1.5 * mean_motion_sat * Ra[j] ** 3 * (J3 / (eccentricity[j] * e[j] ** 3)) * (np.sin(arg_perigee_ep) / np.sin(inclination_ep)) * ((1.25 * np.sin(inclination_ep) ** 2 - 1.0) * np.sin(inclination_ep) ** 2 + eccentricity[j] ** 2 * (1 - (35.0 / 4.0) * np.sin(inclination_ep) ** 2 * np.cos(inclination_ep) ** 2))

    # Combine effects of Earth and lunar oblateness
    dwdt_total = [i + j for i, j in zip(dwdt_oblate, dwdt_earth)]
    dedt_total = [i + j for i, j in zip(dedt_oblate, dedt_earth)]

    perturbations = [dwdt_total, dedt_total, didt_oblate]

    return perturbations


def determine_constellation_size(eccentricity, max_constellation_size, study_name):

    from DVP_general_SPS_functions import convert_string_to_datetime, parse_csv_to_array, invert_events_list
    import sympy
    from sympy import cos
    import numpy as np

    # Initialization
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Read in total active times
    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive")
    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)
    total_target_eclipse = np.sum(target_eclipse[2])

    # Calculate SPS constellation size required
    number_of_sps = []
    for i in range(len(total_active_time)):
        if math.isnan(total_active_time[i]):
            number_of_sps.append(0)
        else:
            sps_num = int(total_target_eclipse / total_active_time[i])
            if sps_num > max_constellation_size:
                number_of_sps.append(max_constellation_size)
            else:
                number_of_sps.append(int(total_target_eclipse / total_active_time[i]))

    plt.plot(number_of_sps)
    plt.show()

    # Calculate distribution of SPS in true anomaly
    num_sps = range(2, max(number_of_sps) + 1)

    return set_constellation_size(eccentricity, num_sps, study_name)


def set_constellation_size(eccentricity, num_sps, study_name):
    # Calculate distribution of SPS in true anomaly
    assert num_sps > 0
    unique_num_sps = range(num_sps + 1)

    if 'SouthPole' in study_name:
        angular_distribution = {}
        for i in unique_num_sps:
            angular_distribution['{}sps'.format(i)] = []
            for j in range(1, i):
                true_anomaly = j * (2 * np.pi / i)
                E = sympy.Symbol('E')
                if j <= (i / 2):
                    eccentric_anomaly = float(min(sympy.solve(cos(true_anomaly) - (cos(E) - eccentricity[i])/(1 - eccentricity[i] * cos(E)), E)))
                elif j > (i / 2):
                    eccentric_anomaly = float(max(sympy.solve(cos(true_anomaly) - (cos(E) - eccentricity[i])/(1 - eccentricity[i] * cos(E)), E)))
                angular_distribution['{}sps'.format(i)].append(round((180.0 / np.pi) * eccentric_anomaly - eccentricity[i] * np.sin(eccentric_anomaly), 4))

    elif 'Equatorial' in study_name:
        angular_distribution = {}
        for i in unique_num_sps:
            angular_distribution['{}sps'.format(i)] = []
            for j in range(0, i):
                arg_perigee = j * (2 * np.pi / i)
                angular_distribution['{}sps'.format(i)].append(round(arg_perigee, 4))


    return num_sps, angular_distribution
