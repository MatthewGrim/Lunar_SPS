"""
Author: Rohan Ramasamy
Date: 04/09/2018

This function processes the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the polar lunar SPS configuration based on feasible pointing and power
constraints.
"""

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def process_stk_data(max_constellation_size, study_name, constellation_variable, start, end, **kwargs):
    """
    max_constellation_size: number of satellites in the constellation
    study_name: name of the study being processed
    oonstellation_variable: angle varied for particular constellation type. For equatorial this is the argument of
                            perigee (argperi). For south pole orbits this is the mean anomaly (meananom).
    start: list representing the start time of the simulation
    end: list representing the end time of the simulation
    """
    total_duration = (end - start).total_seconds()

    # get orbital elements
    max_perigee = kwargs.get('max_perigee', 5000.0)
    max_apogee = kwargs.get('max_apogee', 5000.0)
    min_perigee = kwargs.get('min_perigee', 800.0)
    resolutions = kwargs.get('resolutions', None)
    thresholds = kwargs.get('thresholds', None)
    sma, ecc, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee,
                                                                         min_perigee=min_perigee,
                                                                         resolutions=resolutions,
                                                                         thresholds=thresholds)

    # Get angular distributions for each sps
    arg_perigee = set_constellation_size(max_constellation_size)
    key = '{}sps'.format(max_constellation_size)
    arg_perigees = arg_perigee[key]

    # Get pathway to SPS data directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(os.path.dirname(issue_folder))
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
    data['total_station_keeping'] = np.zeros((len(orbit_data) - 1))
    data['total_stored_power_time'] = np.zeros((len(orbit_data) - 1))
    data['max_stored_power_time'] = np.zeros((len(orbit_data) - 1))
    data['mean_stored_power_time'] = np.zeros((len(orbit_data) - 1))
    data['mean_range'] = np.zeros((len(orbit_data) - 1))
    data['mean_max_range'] = np.zeros((len(orbit_data) - 1))
    data['mean_min_range'] = np.zeros((len(orbit_data) - 1))

    for i in range(1, len(orbit_data)):
        print('Progress: {}%'.format(round(100.0 * (i - 1) / (len(orbit_data) - 2), 2)))
        print("Perigee radius: {} km, Apogee radius: {} km".format(orbit_data[i][0], orbit_data[i][1]))
        sps_active_total = None

        total_time = 0.0
        for j, arg in enumerate(arg_perigees):
            sim_name = 'DVP_{}_{}perigee{}apogee_{}{}'.format(study_name,
                                                                   orbit_data[i][0],
                                                                   orbit_data[i][1],
                                                                   arg,
                                                                   constellation_variable)

            # Import access and lighting for SPS
            sps_lighting = parse_csv_to_array('{}/{}_lighting.csv'.format(stk_data_path, sim_name), start)
            sps_access = parse_csv_to_array('{}/{}_access.csv'.format(stk_data_path, sim_name), start)

            # Get active time of individual satellite and combine events with total
            sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
            if j == 0:
                sps_active_total = sps_active
            else:
                sps_active_total = combine_events(sps_active_total, sps_active)

            # Get the mean, min and max range for the satellite
            total_time += np.sum(sps_active[2])
            sps_range = import_range_data_statistics('{}_range'.format(sim_name), stk_data_path)
            data['mean_range'][i - 1] += np.sum([r * duration for r, duration in zip(sps_range[2], sps_active[2])])
            data['mean_min_range'][i - 1] += np.sum([r * duration for r, duration in zip(sps_range[0], sps_active[2])])
            data['mean_max_range'][i - 1] += np.sum([r * duration for r, duration in zip(sps_range[1], sps_active[2])])

            # Determine station-keeping/battery-charging events
            sps_station_keeping_events = determine_battery_chargeup_events(sps_lighting, sps_access, total_duration)
            data['total_station_keeping'][i - 1] += np.sum(sps_station_keeping_events[2])

            # Stored power events
            sps_eclipse = invert_events_list(sps_lighting, total_duration)
            sps_use_stored_power = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)
            # Check if stored power events exist, if not insert nan
            if np.sum(sps_use_stored_power) == 0.0:
                pass
            else:
                data['total_stored_power_time'][i - 1] += np.sum(sps_use_stored_power[2])
                data['max_stored_power_time'][i - 1] += max(sps_use_stored_power[2])
                data['mean_stored_power_time'][i - 1] += np.mean(sps_use_stored_power[2])

        # Divide by total time across all satellites to get average link requirements
        data['mean_range'][i - 1] /= total_time
        data['mean_min_range'][i - 1] /= total_time
        data['mean_max_range'][i - 1] /= total_time

        # Divide by number of satellites to get average satellite requirements
        data['total_station_keeping'][i - 1] /= max_constellation_size
        data['total_stored_power_time'][i - 1] /= max_constellation_size
        data['max_stored_power_time'][i - 1] /= max_constellation_size
        data['mean_stored_power_time'][i - 1] /= max_constellation_size

        # Determine active time statistics
        data['total_active_time'].append(np.sum(sps_active_total[2]))
        data['max_active_time'].append(np.max(sps_active_total[2]))
        data['mean_active_time'].append(np.mean(sps_active_total[2]))
        data['min_active_time'].append(np.min(sps_active_total[2]))
        data['std_active_time'].append(np.std(sps_active_total[2]))

        # Determine the target blackout statistics
        target_blackout = determine_blackout_data(sps_active_total, target_eclipse, total_duration)
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

    # STATION KEEPING
    write_data_to_file(stk_data_path, study_name, data['total_station_keeping'], "TotalStationKeeping_{}".format(postfix_name))

    # STORED POWER
    write_data_to_file(stk_data_path, study_name, data['total_stored_power_time'], "TotalStoredPower_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['max_stored_power_time'], "MaxStoredPower_{}".format(postfix_name))
    write_data_to_file(stk_data_path, study_name, data['mean_stored_power_time'], "MeanStoredPower_{}".format(postfix_name))

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

