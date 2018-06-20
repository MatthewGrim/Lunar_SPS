""""
01/06/2018
Author: Darian van Paridon

This script takes the processed data for the Lunar SPS for south pole target orbit parametric scan, and applies
constraints relevant to the SPS and target systems to remove data points which do no the meet the requirements. Then
the transmitter aperture size is varied in order to optimize the link efficiency within the remaining design space.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *


def calculate_link_eff(trans_radius):

    # INITIALIZATION
    ####################################################################################################################
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name study
    study_name = 'SouthPole_IncrementedRes_Inertial'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
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
    # Team ITALIA AMALIA (intermediate)
    rec_radius = 0.5
    rec_efficiency = 0.40
    operation_pwr = 100.0
    hibernation_pwr = 7.0
    battery_capacity = 100.0

    # ispace Sorato (miniature)
    # rec_radius = 0.1
    # rec_efficiency = 0.40
    # operation_pwr = 21.5
    # hibernation_pwr = 4.5
    # battery_capacity = 38

    # NASA Curiosity (large)
    # rec_radius = 1
    # rec_efficiency = 0.40
    # operation_pwr = 270
    # hibernation_pwr = 23.5
    # battery_capacity = 1600.0
    ####################################################################################################################

    # SPECIFY DESIGN CONSTRAINTS
    ####################################################################################################################
    # Minimum pointing error of SPS system in radians
    point_error = 1e-6
    # Minimum reduction in overall blackout time in percent
    min_active_time = 18.0
    # Minimum power requirement at target in Watts
    min_power = operation_pwr
    # Maximum time rover can survive without recharging in hours
    rover_blackout_limit = battery_capacity / hibernation_pwr
    ####################################################################################################################

    # READ IN DATA FILES
    ####################################################################################################################
    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive_Inertial_Extended")
    total_blackout_time = read_data_from_file(stk_data_path, study_name, "TotalBlackout_Inertial_Extended")
    max_active_time = read_data_from_file(stk_data_path, study_name, "MaxActive_Inertial_Extended")
    max_blackout_time = read_data_from_file(stk_data_path, study_name, "MaxBlackout_Inertial_Extended")
    mean_range = read_data_from_file(stk_data_path, study_name, "MeanRange_Inertial_Extended")
    mean_active_time = read_data_from_file(stk_data_path, study_name, "MeanActive_Inertial_Extended")
    mean_blackout_time = read_data_from_file(stk_data_path, study_name, "MeanBlackout_Inertial_Extended")
    ####################################################################################################################

    # CALCULATE LINK EFFICIENCY, MEAN POWER, and TOTAL ENERGY DELIVERED
    ####################################################################################################################
    mean_link_efficiency = []
    mean_power_received = []
    for i in range(len(mean_range)):
        # Minimum beam radius as defined by pointing error
        min_beam_radius = rec_radius + (point_error * mean_range[i] * 1000.0)
        # Actual beam radius as defined by Gaussian beam divergence
        surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (mean_range[i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
        # Calculate link efficiency
        if surf_beam_radius <= min_beam_radius:
            mean_link_efficiency.append(1.0)
        else:
            mean_link_efficiency.append((rec_radius / surf_beam_radius) ** 2)
        # Calculate mean power received at target
        mean_power_received.append(mean_link_efficiency[i] * rec_efficiency * trans_power)
    ####################################################################################################################

    # ENFORCE POINTING CONSTRAINTS
    ####################################################################################################################
    for i in range(len(mean_link_efficiency)):
        # Minimum beam radius as defined by pointing error
        min_beam_radius = rec_radius + (point_error * mean_range[i] * 1000.0)
        # Actual beam radius as defined by Gaussian beam divergence
        surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (mean_range[i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
        # Check pointing error constraint
        if surf_beam_radius < min_beam_radius:
            mean_link_efficiency[i] = 0.0

            # Apply effect of constraint to other data sets
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            max_active_time[i] = np.nan
            max_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_power_received[i] = np.nan
    ####################################################################################################################

    # ENFORCE POWER CONSTRAINT
    ####################################################################################################################
    for i in range(len(mean_power_received)):
        if mean_power_received[i] < min_power:
            mean_power_received[i] = np.nan

            # Apply effect of constraint to other data sets
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            max_active_time[i] = np.nan
            max_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_link_efficiency[i] = 0.0
    ####################################################################################################################

    # ENFORCING BLACKOUT DURATION CONSTRAINTS
    ####################################################################################################################
    # Remove data points for which blackout durations exceed the limit
    for i in range(len(max_blackout_time)):
        if max_blackout_time[i] > rover_blackout_limit * 3600.0:
            max_blackout_time[i] = np.nan

            # Apply constraint to all other data sets
            max_active_time[i] = np.nan
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_link_efficiency[i] = 0.0
            mean_power_received[i] = np.nan
    ####################################################################################################################

    # ENFORCE ACTIVE DURATION CONSTRAINT
    ####################################################################################################################
    # Remove data points for which the overall blackout time is not sufficiently reduced
    for i in range(len(total_active_time)):
        if (100.0 * total_active_time[i] / total_duration) < min_active_time:
            total_active_time[i] = np.nan

            # Apply constraint to all other data sets
            max_active_time[i] = np.nan
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_link_efficiency[i] = 0.0
            mean_power_received[i] = np.nan
    ####################################################################################################################

    # Find best link efficiency
    best_orbit_idx = np.nanargmax(mean_link_efficiency)

    return 1 - mean_link_efficiency[best_orbit_idx]


def optimize_link_efficiency():

    from scipy.optimize import minimize_scalar

    optimum = minimize_scalar(calculate_link_eff, bounds=(0, 0.31), method='bounded')
    print(optimum)
    print("Best mean link efficiency: {} %".format(round(100.0 * (1 - optimum.fun), 2)))

    return optimum