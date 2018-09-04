""""
21/06/2018
Author: Darian van Paridon

This script contains functions used for the SPS constrained design tool.

"""


def rover_metrics(rover_name):

    rover = {}

    if rover_name == "amalia":
        # Team ITALIA AMALIA (intermediate)
        rover['rec_radius'] = 0.5
        rover['rec_efficiency'] = 0.40
        rover['operation_pwr'] = 100.0
        rover['hibernation_pwr'] = 7.0
        rover['battery_capacity'] = 100.0

    elif rover_name == "sorato":
        # ispace Sorato (miniature)
        rover['rec_radius'] = 0.1
        rover['rec_efficiency'] = 0.40
        rover['operation_pwr'] = 21.5
        rover['hibernation_pwr'] = 4.5
        rover['battery_capacity'] = 38.0

    elif rover_name == "curiosity":
        # NASA Curiosity (large)
        rover['rec_radius'] = 1.0
        rover['rec_efficiency'] = 0.40
        rover['operation_pwr'] = 270.0
        rover['hibernation_pwr'] = 23.5
        rover['battery_capacity'] = 1600.0
    else:
        print('Invalid rover name. Valid names: amalia, sorato, curiosity')

    return rover


def trans_metrics(selection):

    transmitter = {}

    if selection == '100kW':
        # IPG YLS10000
        transmitter['wavelength'] = 1070e-9
        transmitter['power'] = 100e3
        transmitter['mass'] = 3600.0
        transmitter['efficiency'] = 0.35

    elif selection == '15kW':
        # IPG YLS-CUT
        transmitter['wavelength'] = 1070e-9
        transmitter['power'] = 15e3
        transmitter['mass'] = 440.0
        transmitter['efficiency'] = 0.35

    elif selection == '4kW':
        # Fujikura
        # Mass is estimated based on specific power of IPG lasers
        transmitter['wavelength'] = 1080e-9
        transmitter['power'] = 4e3
        transmitter['mass'] = 150.0
        transmitter['efficiency'] = 0.26

    else:
        print('Select either \high power/ or \low power/')

    return transmitter


def enforce_constraints(data_set, data_type, constraints, constraint_name, constraint_type):

    import numpy as np

    if constraint_type == 'max':
        for i in range(len(data_set[data_type])):
            if data_set[data_type][i] > constraints[constraint_name]:
                for j in data_set:
                    data_set[j][i] = np.nan

    elif constraint_type == 'min':
        for i in range(len(data_set[data_type])):
            if data_set[data_type][i] < constraints[constraint_name]:
                for j in data_set:
                    data_set[j][i] = np.nan

    return data_set


def sort_data_lists(data_set, orbit_data):

    from Lunar_SPS.DVP_Programmatic_Functions import sort_incremented_resolution_data

    data_set_sorted = {}

    for j in data_set:
        data_set_sorted[j] = sort_incremented_resolution_data(orbit_data, data_set[j])

    return data_set_sorted
