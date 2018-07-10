""""
21/06/2018
Author: Darian van Paridon

This script contains functions used for the SPS constrained design tool.

"""


def study_initialization(study_name):

    from DVP_general_SPS_functions import convert_string_to_datetime
    from DVP_Programmatic_Functions import vary_orbital_elements, vary_orbital_elements_incrementing_resolution

    study = {}
    if study_name == 'Brandhorst_1000.0kmRes':
        study['start'] = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
        study['end'] = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
        study['duration'] = (study['end'] - study['start']).total_seconds()

        # Difference in subsequent apogee/perigee radii
        resolution = 1000.0
        max_perigee = 10000.0
        min_perigee = 1000.0
        max_apogee = 54000.0

        # Get orbital data
        semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements(resolution, min_perigee, max_perigee, max_apogee)
        study['semi-maj-axis'] = semi_maj_axis
        study['eccentricity'] = eccentricity
        study['orbits'] = orbit_data
        study['inclination'] = 0.0
        study['arg_perigee'] = 0.0

    elif study_name == 'SouthPole_IncrementedRes_Inertial':
        study['start'] = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
        study['end'] = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
        study['duration'] = (study['end'] - study['start']).total_seconds()

        # Set bounds on parametric scan
        max_perigee = 5000.0
        max_apogee = 5000.0

        # Get orbit data set
        semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)
        study['semi-maj-axis'] = semi_maj_axis
        study['eccentricity'] = eccentricity
        study['orbits'] = orbit_data
        study['inclination'] = 90.0
        study['arg_perigee'] = 90.0

    elif study_name == 'Equatorial_IncrementedRes':
        study['start'] = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
        study['end'] = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
        study['duration'] = (study['end'] - study['start']).total_seconds()

        # Set bounds on parametric scan
        max_perigee = 5000.0
        max_apogee = 5000.0

        # Get orbit data set
        semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)
        study['semi-maj-axis'] = semi_maj_axis
        study['eccentricity'] = eccentricity
        study['orbits'] = orbit_data
        study['inclination'] = 0.0
        study['arg_perigee'] = 0.0

    else:
        print('Invalid study name')

    return study


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


def sort_data_lists(data_set, orbit_data, study_name):

    from DVP_Programmatic_Functions import sort_incremented_resolution_data, sort_data_list_into_array

    data_set_sorted = {}

    if study_name == 'Brandhorst_1000.0kmRes':
        for j in data_set:
            data_set_sorted[j] = sort_data_list_into_array(orbit_data, 1000.0, data_set[j])
    elif study_name == 'SouthPole_IncrementedRes_Inertial' or 'Equatorial_IncrementedRes':
        for j in data_set:
            data_set_sorted[j] = sort_incremented_resolution_data(orbit_data, data_set[j])

    return data_set_sorted
