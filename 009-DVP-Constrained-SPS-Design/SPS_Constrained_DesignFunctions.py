""""
21/06/2018
Author: Darian van Paridon

This script contains functions used for the SPS constrained design tool.

"""


def sps_constraints(rover):

    constraints = {}
    active_constraints = {}

    # Minimum pointing error of SPS system in radians
    constraints['point_error'] = 1e-6
    # Minimum reduction in overall blackout time in percent
    constraints['min_active_time'] = 33.0
    # Minimum power requirement at target in Watts
    constraints['min_power'] = rover['operation_pwr']
    # Maximum time rover can survive without recharging in hours
    constraints['max_blackout'] = rover['battery_capacity'] / rover['hibernation_pwr']
    # Maximum allowable skew in argument of perigee, in degrees per year
    constraints['max_arg_perigee_skew'] = 113.0

    # Specify which constraints are active
    active_constraints['point_error'] = 1
    active_constraints['min_active_time'] = 1
    active_constraints['min_power'] = 1
    active_constraints['max_blackout'] = 1
    active_constraints['max_arg_perigee_skew'] = 0

    return constraints, active_constraints


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

    if selection == 'high power':
        # IPG YLS10000
        transmitter['wavelength'] = 1070e-9
        transmitter['power'] = 100e3
        transmitter['mass'] = 3600.0
        transmitter['efficiency'] = 0.35

    elif selection == 'low power':
        # IPG YLS-CUT
        transmitter['wavelength'] = 1070e-9
        transmitter['power'] = 15e3
        transmitter['mass'] = 440.0
        transmitter['efficiency'] = 0.35
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

    from DVP_Programmatic_Functions import sort_incremented_resolution_data

    data_set_sorted = {}

    for j in data_set:
        data_set_sorted[j] = sort_incremented_resolution_data(orbit_data, data_set[j])

    return data_set_sorted
