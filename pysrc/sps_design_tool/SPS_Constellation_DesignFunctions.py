""""
21/06/2018
Author: Darian van Paridon

This script contains functions used for the SPS constrained design tool.
"""

import re
import numpy as np
import os

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import read_data_from_file
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import convert_string_to_datetime
from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import vary_orbital_elements_incrementing_resolution
from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import sort_incremented_resolution_data
from Lunar_SPS.pysrc.utils.physical_constants import PhysicalConstants


def study_initialization(study_name, **kwargs):
    study = dict()
    study['start'] = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    study['end'] = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    study['duration'] = (study['end'] - study['start']).total_seconds()

    # Get orbit data set
    max_perigee = kwargs.get('max_perigee', 5000.0)
    max_apogee = kwargs.get('max_apogee', 5000.0)
    min_perigee = kwargs.get('min_perigee', 800.0)
    resolutions = kwargs.get('resolutions', None)
    thresholds = kwargs.get('thresholds', None)
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee,
                                                                                            min_perigee=min_perigee,
                                                                                            resolutions=resolutions,
                                                                                            thresholds=thresholds)

    if 'NorthPole' in study_name:
        study['semi-maj-axis'] = semi_maj_axis
        study['eccentricity'] = eccentricity
        study['orbits'] = orbit_data
        study['inclination'] = 90.0
        study['arg_perigee'] = 90.0
        if 'ManyToOne' in study_name:
            study['constellation_variable'] = 'meananom'

    elif 'Equatorial' in study_name:
        study['semi-maj-axis'] = semi_maj_axis
        study['eccentricity'] = eccentricity
        study['orbits'] = orbit_data
        study['inclination'] = 0.0
        study['arg_perigee'] = 0.0
        if 'ManyToOne' in study_name:
            study['constellation_variable'] = 'argperi'

    else:
        print('Invalid study name')

    return study


def rover_metrics(rover_name):
    def approximate_rec_radius(rover, receiver_solar_efficiency=0.2):
        solar_intensity = PhysicalConstants.solar_irradiance
        receiver_area = rover['daytime_power'] / (solar_intensity * receiver_solar_efficiency)
        return np.sqrt(receiver_area / np.pi)

    rover = {}
    if "amalia" in rover_name:
        # Team ITALIA AMALIA (intermediate)
        rover['operation_pwr'] = 100.0
        rover['daytime_power'] = 100.0
        rover['rec_efficiency'] = 0.5
        rover['hibernation_pwr'] = 7.0
        rover['battery_capacity'] = 100.0
        rover['rec_radius'] = approximate_rec_radius(rover,  receiver_solar_efficiency=0.2)
    elif "sorato" in rover_name:
        # ispace Sorato (miniature)
        rover['operation_pwr'] = 21.5
        rover['daytime_power'] = 21.5
        rover['rec_efficiency'] = 0.5
        rover['hibernation_pwr'] = 4.5
        rover['battery_capacity'] = 38.0
        rover['rec_radius'] = approximate_rec_radius(rover)
    elif "rover_concept_five" in rover_name:
        # Designed for 1 microradian
        rover['operation_pwr'] = 240.0
        rover['daytime_power'] = 150.0
        rover['rec_efficiency'] = 0.5
        rover['hibernation_pwr'] = 48.0
        rover['battery_capacity'] = 300.0
        rover['rec_radius'] = approximate_rec_radius(rover, receiver_solar_efficiency=0.06)
    elif "rover_concept_six" in rover_name:
        # Designed for 1 microradian
        rover['operation_pwr'] = 280.0
        rover['daytime_power'] = 150.0
        rover['rec_efficiency'] = 0.5
        rover['hibernation_pwr'] = 28.0
        rover['battery_capacity'] = 300.0
        rover['rec_radius'] = approximate_rec_radius(rover, receiver_solar_efficiency=0.1)
    elif "rover_concept_seven" in rover_name:
        # Designed for 1 microradian
        rover['operation_pwr'] = 200.0
        rover['daytime_power'] = 150.0
        rover['rec_efficiency'] = 0.5
        rover['hibernation_pwr'] = 20.0
        rover['battery_capacity'] = 300.0
        rover['rec_radius'] = approximate_rec_radius(rover, receiver_solar_efficiency=0.1)
    else:
        print('Invalid rover name: {}. Valid names: amalia, sorato, curiosity'.format(rover_name))

    return rover


def trans_metrics(selection):

    transmitter = {}

    if '100kW' in selection:
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
        transmitter['wavelength'] = 1070e-9
        transmitter['power'] = 4e3
        transmitter['mass'] = 150.0
        transmitter['efficiency'] = 0.26

    else:
        print('Select either 100kW, 15kW, or 4kW')

    return transmitter


def enforce_constraints(data_set, data_type, constraints, constraint_name, constraint_type):

    import numpy as np
    import sys

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

    if np.isnan(data_set[data_type]).all():
        print('No feasible designs remaining, {} constraint too restrictive.'.format(constraint_name))
        sys.exit()

    return data_set


def determine_number_of_sps_for_active_time(stk_data_path, study_name, study, constraints, active_constraints, max_num_sps):

    # Determine what size of SPS constellation is necessary in order to achieve the total active time constraint
    num_sps = 0
    all_nans = True
    while all_nans:
        num_sps += 1
        if num_sps == max_num_sps:
            break
        else:
            data_set = {}
            data_set['total_active_time'] = read_data_from_file(stk_data_path, study_name, "TotalActive_{}SPS".format(num_sps))

            # Apply active constraints
            if active_constraints['min_active_time'] == 1:
                # Remove data points for which the overall blackout time is not sufficiently reduced
                for i in range(len(data_set['total_active_time'])):
                    if (100.0 * data_set['total_active_time'][i] / study['duration']) < constraints['min_active_time']:
                        data_set['total_active_time'][i] = np.nan
            else:
                pass

        all_nans = np.isnan(data_set['total_active_time']).all()

    return num_sps


def read_in_processed_data_reports(stk_data_path, study_name, num_sps, use_storage):
    data_set = {}
    if use_storage:
        data_set['total_active_time'] = read_data_from_file(stk_data_path, study_name, "TotalActiveWithStorage_{}SPS".format(num_sps))
        data_set['max_active_time'] = read_data_from_file(stk_data_path, study_name, "MaxActiveWithStorage_{}SPS".format(num_sps))
        data_set['mean_active_time'] = read_data_from_file(stk_data_path, study_name, "MeanActiveWithStorage_{}SPS".format(num_sps))
        data_set['min_active_duration'] = read_data_from_file(stk_data_path, study_name, 'MinActiveWithStorage_{}SPS'.format(num_sps))

        data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackoutWithStorage_{}SPS".format(num_sps))
        data_set['max_blackout_duration'] = read_data_from_file(stk_data_path, study_name, "MaxBlackoutWithStorage_{}SPS".format(num_sps))
        data_set['mean_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MeanBlackoutWithStorage_{}SPS".format(num_sps))
    else:
        data_set['total_active_time'] = read_data_from_file(stk_data_path, study_name, "TotalActive_{}SPS".format(num_sps))
        data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackout_{}SPS".format(num_sps))
        data_set['max_active_time'] = read_data_from_file(stk_data_path, study_name, "MaxActive_{}SPS".format(num_sps))
        data_set['mean_active_time'] = read_data_from_file(stk_data_path, study_name, "MeanActive_{}SPS".format(num_sps))
        data_set['min_active_duration'] = read_data_from_file(stk_data_path, study_name, 'MinActive_{}SPS'.format(num_sps))

        data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackout_{}SPS".format(num_sps))
        data_set['max_blackout_duration'] = read_data_from_file(stk_data_path, study_name, "MaxBlackout_{}SPS".format(num_sps))
        data_set['mean_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MeanBlackout_{}SPS".format(num_sps))

    data_set['mean_range'] = np.loadtxt(os.path.join(stk_data_path, "MeanRange_{}SPS_{}.txt".format(num_sps, study_name)))
    data_set['max_range'] = np.loadtxt(os.path.join(stk_data_path, "MeanMaxRange_{}SPS_{}.txt".format(num_sps, study_name)))
    data_set['min_range'] = np.loadtxt(os.path.join(stk_data_path, "MeanMinRange_{}SPS_{}.txt".format(num_sps, study_name)))

    # STORED POWER
    data_set['total_stored_power_time'] = np.loadtxt(os.path.join(stk_data_path, "TotalStoredPower_{}SPS_{}.txt".format(num_sps, study_name)))
    data_set['max_stored_power_time'] = np.loadtxt(os.path.join(stk_data_path, "MaxStoredPower_{}SPS_{}.txt".format(num_sps, study_name)))
    data_set['mean_stored_power_time'] = np.loadtxt(os.path.join(stk_data_path, "MeanStoredPower_{}SPS_{}.txt".format(num_sps, study_name)))

    data_set['mean_range'] = np.loadtxt(os.path.join(stk_data_path, "MeanRange_{}SPS_{}.txt".format(num_sps, study_name)))
    data_set['max_range'] = np.loadtxt(os.path.join(stk_data_path, "MeanMaxRange_{}SPS_{}.txt".format(num_sps, study_name)))
    data_set['min_range'] = np.loadtxt(os.path.join(stk_data_path, "MeanMinRange_{}SPS_{}.txt".format(num_sps, study_name)))

    return data_set


def sort_data_lists(data_set, transmitter, orbit_data, study_name, **kwargs):
    data_set_sorted = {}

    if 'IncrementedRes' in study_name:
        resolutions = kwargs.get('resolutions', None)
        thresholds = kwargs.get('thresholds', None)
        for j in data_set:
            data_set_sorted[j] = sort_incremented_resolution_data(orbit_data, data_set[j],
                                                                  resolution=resolutions, thresholds=thresholds)
        transmitter['radius'] = sort_incremented_resolution_data(orbit_data, transmitter['radius'],
                                                                 resolution=resolutions, thresholds=thresholds)
    else:
        raise RuntimeError('This code does not handle data sets without incremented resolutions')

    # Extract unique perigees and apogees tested for plotting
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    r_moon = PhysicalConstants.r_moon
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    # Reduce perigee and apogee to altitudes instead of radii
    perigee_altitudes = [i - r_moon for i in unique_perigees]
    apogee_altitudes = [i - r_moon for i in unique_apogees]

    return data_set_sorted, transmitter, perigee_altitudes, unique_perigees, apogee_altitudes, unique_apogees


def calculate_link_efficiency_and_power_delivered_for_single_rover(rover, data_set, transmitter, constraints, active_constraints, duration, include_misalignment):
    # Initialize lists
    data_set['min_link_efficiency'] = []
    data_set['min_power_received'] = []
    data_set['mean_link_efficiency'] = []
    data_set['mean_power_received'] = []
    # Cycle through all access events
    for i in range(0, len(data_set['max_range'])):
        # Calculate average misalignment loss
        if include_misalignment:
            theta_offset = np.asarray([np.arcsin(r_moon / (np.sqrt(2) * data_set['min_range'][i] * 1000.0)),
                                        np.arcsin(r_moon / (np.sqrt(2) * data_set['max_range'][i] * 1000.0)),
                                        np.arcsin(r_moon / (np.sqrt(2) * data_set['mean_range'][i] * 1000.0))])
            misalignment_loss = np.cos(theta_offset)
        else:
            misalignment_loss = np.ones((3,))

        # If range is "np.nan", no access periods exist for that orbit, therefore treat orbit as infeasible design
        if data_set['max_range'][i] == np.nan:
            data_set['min_link_efficiency'].append(np.nan)
            data_set['mean_link_efficiency'].append(np.nan)
            data_set['min_power_received'].append(np.nan)
            data_set['mean_power_received'].append(np.nan)
        # Otherwise calculate efficiency and power, apply pointing constraints
        else:
            # Minimum beam radius as defined by pointing error, at minimum [0] and maximum [1] and mean [2] ranges
            min_beam_radius = [rover['rec_radius'] + (constraints['point_error'] * data_set['min_range'][i] * 1000.0),
                               rover['rec_radius'] + (constraints['point_error'] * data_set['max_range'][i] * 1000.0),
                               rover['rec_radius'] + (constraints['point_error'] * data_set['mean_range'][i] * 1000.0)]

            # Actual beam radius as defined by Gaussian beam divergence. idx [0] = min, idx [1] = max, idx [2] = mean
            surf_beam_radius = [transmitter['radius'][i] * np.sqrt(1 + (transmitter['wavelength'] * (data_set['min_range'][i] * 1000.0) / (np.pi * transmitter['radius'][i] ** 2)) ** 2),
                                transmitter['radius'][i] * np.sqrt(1 + (transmitter['wavelength'] * (data_set['max_range'][i] * 1000.0) / (np.pi * transmitter['radius'][i] ** 2)) ** 2),
                                transmitter['radius'][i] * np.sqrt(1 + (transmitter['wavelength'] * (data_set['mean_range'][i] * 1000.0) / (np.pi * transmitter['radius'][i] ** 2)) ** 2)]

            # IF constraint is active, remove data points which violate pointing constraint
            # Checking pointing constraint at minimum mean and maximum ranges
            if active_constraints['point_error'] == 1:
                if surf_beam_radius[0] < min_beam_radius[0] or surf_beam_radius[1] < min_beam_radius[1] or surf_beam_radius[2] < min_beam_radius[2]:
                    data_set['min_link_efficiency'].append(np.nan)
                    data_set['mean_link_efficiency'].append(np.nan)
                    data_set['min_power_received'].append(np.nan)
                    data_set['mean_power_received'].append(np.nan)
                    # Apply constraint all other data lists
                    for j in data_set:
                        data_set[j][i] = np.nan
                else:
                    data_set['min_link_efficiency'].append((rover['rec_radius'] / surf_beam_radius[1]) ** 2 * misalignment_loss[1])
                    data_set['mean_link_efficiency'].append((rover['rec_radius'] / surf_beam_radius[2]) ** 2 * misalignment_loss[2])
                    data_set['min_power_received'].append(data_set['min_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])
                    data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])

            # ELSE calculate as normal
            else:
                # Calculate minimum link efficiency, as well as power delivered based on max surface beam size
                if surf_beam_radius[1] <= rover['rec_radius']:
                    data_set['min_link_efficiency'].append(misalignment_loss[1])
                else:
                    data_set['min_link_efficiency'].append((rover['rec_radius'] / surf_beam_radius[1]) ** 2 * misalignment_loss[1])
                data_set['min_power_received'].append(data_set['min_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])
                
                # Calculate mean link efficiency, as well as power delivered based on mean surface beam size
                if surf_beam_radius[2] <= rover['rec_radius']:
                    data_set['mean_link_efficiency'].append(misalignment_loss[2])
                else:
                    data_set['mean_link_efficiency'].append((rover['rec_radius'] / surf_beam_radius[2]) ** 2 * misalignment_loss[2])
                data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])

    # Calculate total energy delivered to receiver based on mean power
    data_set['total_energy'] = []
    data_set['total_energy'] = [i * j * duration / 100.0 for i, j in zip(data_set['mean_power_received'], data_set['total_active_time'])]
    
    return data_set

def optimize_transmitter_radius(data_set, constraints, active_constraints, transmitter, rover):
    """
    This function should determine the optimum transmitter size for every orbit, subject to pointing constraints.

    The function is inherantly conservative, in that it assumes we are beaming at the maximum range. In some applications,
    we may not need to do this, and can improve the link efficiency.
    """
    transmitter['radius'] = np.zeros((len(data_set['max_range'])))
    for i, z_max in enumerate(data_set['max_range']):
        # Optimum from gaussian optics at maximum range
        trans_radius = np.sqrt(z_max * 1000.0 * transmitter['wavelength'] / np.pi)

        # Reduce transmitter size to increase beam size to meet pointing errors
        if active_constraints['point_error']:
            trans_radius = 1.12
            max_iter = 30000
            iterations = 0
            while True:
                break_loop = True
                # Loop through range values to test pointing requirement
                for range_loc in ['min_range', 'max_range', 'mean_range']:
                    surf_beam_radius = trans_radius * np.sqrt(1 + (transmitter['wavelength'] * (data_set[range_loc][i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
                    w_pointing = rover['rec_radius'] + (constraints['point_error'] * data_set[range_loc][i] * 1000.0)
                    if surf_beam_radius < w_pointing:
                        trans_radius *= 0.9999
                        break_loop = False
                    else:
                        pass

                # Break loop if beam is bigger than the pointing error
                if break_loop:
                    break
                
                # If max iterations is exceeded raise a runtime error
                iterations += 1
                if iterations > max_iter:
                    print(max_iter)
                    print(iterations)
                    raise RuntimeError('Failed to find solution! Current radius: {}'.format(trans_radius))

            # Last check that beam pointing is satisfied
            for range_loc in ['min_range', 'max_range', 'mean_range']:
                surf_beam_radius = trans_radius * np.sqrt(1 + (transmitter['wavelength'] * (data_set[range_loc][i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
                w_pointing = rover['rec_radius'] + (constraints['point_error'] * data_set[range_loc][i] * 1000.0)
                assert surf_beam_radius > w_pointing

        transmitter['radius'][i] = trans_radius

    return transmitter

def get_transmitter_power_and_surf_flux(transmitter, rover, sorted_data_set, constraints, active_constraints):
    max_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['max_range'] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)
    min_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['min_range'] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)
    mean_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['mean_range'] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)

    if active_constraints['transmitter_pwr_optimization'] == 1:
        min_allowable_transmitter_pwr = constraints['min_power'] / (sorted_data_set['min_link_efficiency'] * rover['rec_efficiency'])

        # Adjust mean power delivered according to the minimum allowable transmitter power
        sorted_data_set['min_power_received'] = sorted_data_set['min_power_received'] * (min_allowable_transmitter_pwr / transmitter['power'])
        sorted_data_set['mean_power_received'] = sorted_data_set['mean_power_received'] * (min_allowable_transmitter_pwr / transmitter['power'])

        # Assume transmitter operates at minimum allowable power
        transmitter['power'] = min_allowable_transmitter_pwr
    else:
        pass

    # Calculate flux variations
    surf_flux = {}
    surf_flux ['mean'] = transmitter['power'] / (np.pi * mean_surf_beam_radius ** 2)
    surf_flux['min'] = transmitter['power'] / (np.pi * max_surf_beam_radius ** 2)
    surf_flux['max'] = transmitter['power'] / (np.pi * min_surf_beam_radius ** 2)

    return sorted_data_set, surf_flux, transmitter   


def optimize_transmitter_power(transmitter, rover, sorted_data_set, best_orbit_idx, constraints, active_constraints):

    max_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['max_range'][best_orbit_idx] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)
    min_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['min_range'][best_orbit_idx] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)
    mean_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['mean_range'][best_orbit_idx] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)

    if active_constraints['transmitter_pwr_optimization'] == 1:
        # Determine minimum allowable transmitter power based on rover requirements
        min_allowable_transmitter_pwr = constraints['min_power'] * max_surf_beam_radius ** 2 / (rover['rec_efficiency'] * rover['rec_radius'] ** 2)

        # Adjust mean power delivered according to the minimum allowable transmitter power
        sorted_data_set['mean_power_received'] = sorted_data_set['mean_power_received'] * (min_allowable_transmitter_pwr / transmitter['power'])

        # Switch to smallest transmitter which accommodates power
        if min_allowable_transmitter_pwr < 4000.0:
            transmitter = trans_metrics('4kW')
        elif 4000.0 < min_allowable_transmitter_pwr < 15000.0:
            transmitter = trans_metrics('15kW')
        elif min_allowable_transmitter_pwr > 15000.0:
            transmitter = trans_metrics('100kW')

        # Assume transmitter operates at minimum allowable power
        transmitter['power'] = min_allowable_transmitter_pwr
    else:
        pass

    # Calculate flux variations
    surf_flux = {}
    surf_flux ['mean'] = transmitter['power'] / (np.pi * mean_surf_beam_radius ** 2)
    surf_flux['min'] = transmitter['power'] / (np.pi * max_surf_beam_radius ** 2)
    surf_flux['max'] = transmitter['power'] / (np.pi * min_surf_beam_radius ** 2)

    return sorted_data_set, surf_flux, transmitter
