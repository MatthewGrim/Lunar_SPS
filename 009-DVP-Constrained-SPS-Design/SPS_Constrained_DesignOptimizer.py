""""
01/06/2018
Author: Darian van Paridon

This script takes the processed data for the Lunar SPS for south pole target orbit parametric scan, and applies
constraints relevant to the SPS and target systems to remove data points which do no the meet the requirements. Then
the transmitter aperture size is varied in order to optimize the link efficiency within the remaining design space.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *
from SPS_Constrained_DesignFunctions import *


def calculate_link_eff(trans_radius):

    # INITIALIZATION
    ####################################################################################################################
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Name study
    study_name = 'SouthPole_IncrementedRes_Inertial'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Get orbit data
    max_perigee = 5000.0
    max_apogee = 5000.0
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)
    ####################################################################################################################

    # SELECT TRANSMITTER
    ####################################################################################################################
    # IPG YLS10000 (high power)
    wavelength = 1070e-9
    trans_power = 100e3

    # IPG YLS-CUT (low power)
    # wavelength = 1070e-9
    # trans_power = 15e3
    ####################################################################################################################

    # SELECT RECEIVER
    ####################################################################################################################
    rover = rover_metrics('sorato')
    ####################################################################################################################

    # SPECIFY DESIGN CONSTRAINTS
    ####################################################################################################################
    constraints, active_constraints = sps_constraints(rover)
    ####################################################################################################################

    # READ IN DATA FILES
    ####################################################################################################################
    data_set = {"total_active_time": [], 'total_blackout_time': [], 'max_active_time': [],
                'max_blackout_time': [], 'mean_range': [], 'mean_active_time': [], 'mean_blackout_time': []}

    data_set['total_active_time'] = read_data_from_file(stk_data_path, study_name, "TotalActive_Inertial_Extended")
    data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackout_Inertial_Extended")
    data_set['max_active_time'] = read_data_from_file(stk_data_path, study_name, "MaxActive_Inertial_Extended")
    data_set['max_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MaxBlackout_Inertial_Extended")
    data_set['mean_range'] = read_data_from_file(stk_data_path, study_name, "MeanRange_Inertial_Extended")
    data_set['mean_active_time'] = read_data_from_file(stk_data_path, study_name, "MeanActive_Inertial_Extended")
    data_set['mean_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MeanBlackout_Inertial_Extended")
    ####################################################################################################################

    # ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS
    ####################################################################################################################
    # Orbital perturbations on argument of perigee [0], eccentricity [1], and inclination [2]
    perturbations = calculate_orbital_perturbations(semi_maj_axis, eccentricity)
    # Calculate skew in argument of perigee in degrees per year
    data_set['arg_perigee_skew'] = [i * (365.0 * 24.0 * 3600.0) * 180.0 / np.pi for i in perturbations[0]]
    ####################################################################################################################

    # CALCULATE LINK EFFICIENCY, MEAN POWER, and TOTAL ENERGY DELIVERED
    # POINTING ERROR CONSTRAINT ALSO APPLIED
    ####################################################################################################################
    data_set['mean_link_efficiency'] = []
    data_set['mean_power_received'] = []
    for i in range(len(data_set['mean_range'])):
        # Minimum beam radius as defined by pointing error
        min_beam_radius = rover['rec_radius'] + (constraints['point_error'] * data_set['mean_range'][i] * 1000.0)
        # Actual beam radius as defined by Gaussian beam divergence
        surf_beam_radius = trans_radius * np.sqrt(
            1 + (wavelength * (data_set['mean_range'][i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
        # Calculate link efficiency
        if surf_beam_radius < min_beam_radius:
            if active_constraints['point_error'] == 1:
                data_set['mean_link_efficiency'].append(np.nan)
                data_set['mean_power_received'].append(
                    data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * trans_power)
                # apply constraint all other data lists
                for j in data_set:
                    data_set[j][i] = np.nan
            else:
                data_set['mean_link_efficiency'].append(1.0)
                data_set['mean_power_received'].append(
                    data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * trans_power)
        else:
            data_set['mean_link_efficiency'].append((rover['rec_radius'] / surf_beam_radius) ** 2)
            data_set['mean_power_received'].append(
                data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * trans_power)
    # Calculate total energy delivered to receiver
    data_set['total_energy'] = []
    data_set['total_energy'] = [i * j for i, j in zip(data_set['mean_power_received'], data_set['total_active_time'])]
    ####################################################################################################################

    # ENFORCE CONSTRAINTS
    ####################################################################################################################
    # Convert data to appropriate units for applying constraints
    data_set['max_blackout_time'] = [i / 3600.0 for i in data_set['max_blackout_time']]
    data_set['total_active_time'] = [100.0 * i / total_duration for i in data_set['total_active_time']]
    # Remove data points for which not enough power is delivered on average
    if active_constraints['min_power'] == 1:
        data_set = enforce_constraints(data_set, 'mean_power_received', constraints, 'min_power', 'min')
    else:
        pass
    # Remove data points for which blackout durations exceed the limit
    if active_constraints['max_blackout'] == 1:
        data_set = enforce_constraints(data_set, 'max_blackout_time', constraints, 'max_blackout', 'max')
    else:
        pass
    # Remove data points for which the overall blackout time is not sufficiently reduced
    if active_constraints['min_active_time'] == 1:
        data_set = enforce_constraints(data_set, 'total_active_time', constraints, 'min_active_time', 'min')
    else:
        pass
    # Remove data points for which the overall blackout time is not sufficiently reduced
    if active_constraints['max_arg_perigee_skew'] == 1:
        data_set = enforce_constraints(data_set, 'arg_perigee_skew', constraints, 'max_perigee_skew', 'max')
    else:
        pass
    ####################################################################################################################

    # Find best link efficiency
    best_orbit_idx = np.nanargmax(data_set['mean_link_efficiency'])

    return 1 - data_set['mean_link_efficiency'][best_orbit_idx]


def optimize_link_efficiency():

    from scipy.optimize import minimize_scalar

    optimum = minimize_scalar(calculate_link_eff, bounds=(0, 0.3345), method='bounded')
    print('Transmitter aperture radius: {} cm'.format(round(optimum.x * 100.0, 2)))

    return optimum
