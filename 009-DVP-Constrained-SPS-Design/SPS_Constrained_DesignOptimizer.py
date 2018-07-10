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


def calculate_link_eff(trans_radius, args):

    # INITIALIZATION
    ####################################################################################################################
    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Retrieve transmitter metrics
    transmitter = trans_metrics(args[0])
    # Retrieve rover/target metrics
    rover = rover_metrics(args[1])
    # Retrieve constraints
    constraints = args[2]
    active_constraints = args[3]
    # Retrieve study name
    study_name = args[4]
    study = study_initialization(study_name)

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    ####################################################################################################################

    # READ IN DATA FILES
    ####################################################################################################################
    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive")
    max_blackout_time = read_data_from_file(stk_data_path, study_name, "MaxBlackout")
    mean_range = read_data_from_file(stk_data_path, study_name, "MeanRange")
    ####################################################################################################################

    # CALCULATE LINK EFFICIENCY, MEAN POWER, and TOTAL ENERGY DELIVERED
    ####################################################################################################################
    mean_link_efficiency = []
    mean_power_received = []
    for i in range(len(mean_range)):
        if math.isnan(mean_range[i]):
            mean_link_efficiency.append(0.0)
        else:
            # Actual beam radius as defined by Gaussian beam divergence
            surf_beam_radius = trans_radius * np.sqrt(
                1 + (transmitter['wavelength'] * (mean_range[i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
            # Calculate link efficiency
            if surf_beam_radius <= rover['rec_radius']:
                mean_link_efficiency.append(1.0)
            else:
                mean_link_efficiency.append((rover['rec_radius'] / surf_beam_radius) ** 2)
            # Calculate mean power received at target
            mean_power_received.append(mean_link_efficiency[i] * rover['rec_efficiency'] * transmitter['power'])
    ####################################################################################################################

    # ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS
    ####################################################################################################################
    # Orbital perturbations on argument of perigee [0], eccentricity [1], and inclination [2]
    perturbations = calculate_orbital_perturbations(study['semi-maj-axis'], study['eccentricity'], study_name)
    arg_perigee_skew = [abs(i * study['duration'] * 180.0 / np.pi) for i in perturbations[0]]
    ####################################################################################################################

    # ENFORCE POINTING CONSTRAINTS
    ####################################################################################################################
    if active_constraints['point_error'] == 1:
        for i in range(len(mean_link_efficiency)):
            # Minimum beam radius as defined by pointing error
            min_beam_radius = rover['rec_radius'] + (constraints['point_error'] * mean_range[i] * 1000.0)
            # Actual beam radius as defined by Gaussian beam divergence
            surf_beam_radius = trans_radius * np.sqrt(
                1 + (transmitter['wavelength'] * (mean_range[i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
            # Check pointing error constraint
            if surf_beam_radius < min_beam_radius:
                mean_link_efficiency[i] = 0.0
    else:
        pass
    ####################################################################################################################

    # ENFORCE POWER CONSTRAINT
    ####################################################################################################################
    if active_constraints['min_power'] == 1:
        for i in range(len(mean_power_received)):
            if mean_power_received[i] < constraints['min_power']:
                mean_link_efficiency[i] = 0.0
    else:
        pass
    ####################################################################################################################

    # ENFORCING BLACKOUT DURATION CONSTRAINTS
    ####################################################################################################################
    if active_constraints['max_blackout'] == 1:
        # Remove data points for which blackout durations exceed the limit
        for i in range(len(max_blackout_time)):
            if max_blackout_time[i] > constraints['max_blackout'] * 3600.0:
                mean_link_efficiency[i] = 0.0
    else:
        pass
    ####################################################################################################################

    # ENFORCE ACTIVE DURATION CONSTRAINT
    ####################################################################################################################
    if active_constraints['min_active_time'] == 1:
        # Remove data points for which the overall blackout time is not sufficiently reduced
        for i in range(len(total_active_time)):
            if (100.0 * total_active_time[i] / study['duration']) < constraints['min_active_time']:
                mean_link_efficiency[i] = 0.0
    else:
        pass
    ####################################################################################################################

    # ENFORCE ORBITAL STABILITY CONSTRAINT
    ####################################################################################################################
    if active_constraints['max_arg_perigee_skew'] == 1:
        # Remove data points for which the overall blackout time is not sufficiently reduced
        for i in range(len(arg_perigee_skew)):
            if abs(arg_perigee_skew[i]) / 2.0 > constraints['max_arg_perigee_skew']:
                mean_link_efficiency[i] = 0.0
    else:
        pass
    ####################################################################################################################

    # Find best link efficiency
    best_orbit_idx = np.nanargmax(mean_link_efficiency)

    return 1.0 - mean_link_efficiency[best_orbit_idx]


def optimize_link_efficiency(trans_selection, rover_selection, constraints, active_constraints, study_name):

    from scipy.optimize import minimize_scalar

    args = [trans_selection, rover_selection, constraints, active_constraints, study_name]
    optimum = minimize_scalar(calculate_link_eff, bounds=(1e-3, 0.75), method='bounded', args=args)

    return optimum
