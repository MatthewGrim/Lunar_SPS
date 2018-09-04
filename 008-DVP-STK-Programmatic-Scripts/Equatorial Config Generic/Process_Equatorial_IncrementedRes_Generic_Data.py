"""
Author: Rohan Ramasamy
Date: 04/09/2018

This script processes the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the polar lunar SPS configuration based on feasible pointing and power
constraints.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import *
import sys


def process_data(max_constellation_size):
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # get orbital elements
    max_perigee = 5000.0
    max_apogee = 5000.0
    sma, ecc, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee, min_perigee=800.0,
                                                                         resolutions=np.array((50.0, 100.0, 100.0, 250.0)),
                                                                         thresholds=np.array((1000.0, 1500.0, 2500.0)))

    # Get angular distributions for each sps
    arg_perigee = set_constellation_size(max_constellation_size)
    arg_perigees = arg_perigee['{}sps'.format(max_constellation_size)]

    # Get pathway to SPS data directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)
    study_name = 'Equatorial_IncrementedRes_Generic'
    stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

    # Import target illumination events
    target_lighting_raw = os.path.join(stk_data_path, 'DVP_{}_Target_Lighting.csv'.format(study_name))
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Initialize lists
    data = dict()
    data['total_active_time'] = []
    data['total_blackout_time'] = []
    data['max_active_time'] = []
    data['mean_active_time'] = []
    data['min_active_time'] = []
    data['std_active_time'] = []
    data['max_blackout_time'] = []
    data['mean_blackout_time'] = []
    data['min_blackout_time'] = []
    data['std_blackout_time'] = []
    data['mean_range'] = np.zeros((len(arg_perigees), len(orbit_data) - 1))
    data['mean_max_range'] = np.zeros((len(arg_perigees), len(orbit_data) - 1))
    data['mean_min_range'] = np.zeros((len(arg_perigees), len(orbit_data) - 1))

    for i in range(1, len(orbit_data)):
        print('Progress: {}%'.format(round(100.0 * (i - 1) / (len(orbit_data) - 2), 2)))
        print("Perigee radius: {} km, Apogee radius: {} km".format(orbit_data[i][0], orbit_data[i][1]))
        sps_active_total = []

        for j, arg in enumerate(arg_perigees):
            sim_name = 'DVP_{}_{}perigee{}apogee_{}argperi'.format(study_name,
                                                                   orbit_data[i][0],
                                                                   orbit_data[i][1],
                                                                   arg)

            # Import access and lighting for SPS
            sps_lighting = parse_csv_to_array('{}/{}_lighting.csv'.format(stk_data_path, sim_name), start)
            sps_access = parse_csv_to_array('{}/{}_access.csv'.format(stk_data_path, sim_name), start)

            # Get active time of individual satellite and combine events with total
            sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
            sps_active_total = combine_events(sps_active_total, sps_active)

            # Get the mean, min and max range for the satellite
            sps_range = import_range_data_statistics('{}/{}_range'.format(stk_data_path, sim_name), stk_data_path)
            data['mean_range'][j, i] = np.sum([(range * duration) / np.sum(sps_active[2]) for range, duration in zip(sps_range[2], sps_active[2])])
            data['mean_min_range'][j, i] = np.sum([(range * duration) / np.sum(sps_active[2]) for range, duration in zip(sps_range[0], sps_active[2])])
            data['mean_max_range'][j, i] = np.sum([(range * duration) / np.sum(sps_active[2]) for range, duration in zip(sps_range[1], sps_active[2])])

        # Determine active time statistics
        data['total_active_time'].append(np.sum(sps_active[2]))
        data['max_active_time'].append(np.max(sps_active[2]))
        data['mean_active_time'].append(np.mean(sps_active[2]))
        data['min_active_time'].append(np.min(sps_active[2]))
        data['std_active_time'].append(np.std(sps_active[2]))

        # Determine the target blackout statistics
        target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        data['total_blackout_time'].append(np.sum(target_blackout[2]))
        data['max_blackout_time'].append(np.max(target_blackout[2]))
        data['mean_blackout_time'].append(np.mean(target_blackout[2]))
        data['min_blackout_time'].append(np.min(target_blackout[2]))
        data['std_blackout_time'].append(np.std(target_blackout[2]))

    # MEAN RANGE
    postfix_name = "{}SPS".format(max_constellation_size)
    write_data_to_file(stk_data_path, study_name, data['mean_range'], "MeanRange_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['mean_min_range'], "MeanMinRange_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['mean_max_range'], "MeanMaxRange_{}".format(postfix_name))

    # ACTIVE TIME DATA
    write_data_to_file(stk_data_path, study_name, data['total_active_time'], "TotalActive_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['max_active_time'], "MaxActive_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['mean_active_time'], "MeanActive_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['min_active_time'], "MinActive_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['std_active_time'], "StdActive_{}".format(postfix_name))

    # BLACKOUT TIME DATA
    write_data_to_file(stk_data_path, study_name, data['total_blackout_time'], "TotalBlackout_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['max_blackout_time'], "MaxBlackout_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['mean_blackout_time'], "MeanBlackout_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['min_blackout_time'], "MinBlackout_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['std_blackout_time'], "StdBlackout_{}".format(postfix_name))


if __name__ == '__main__':
    process_data(3)

