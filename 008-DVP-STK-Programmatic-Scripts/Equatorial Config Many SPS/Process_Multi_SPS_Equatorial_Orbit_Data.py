""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the polar lunar SPS configuration based on feasible pointing and power
constraints.

There is a section for reading and processing STK data reports, and writing the processed
data to a file. The following section is for reading the processed data files if they have
already been written. Comment out one section or the other as necessary.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *
import sys


def main():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # INITIALIZATION
    ####################################################################################################################
    # Set bounds on parametric scan
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Get orbit data set
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name study
    study_name = 'Equatorial_IncrementedRes_ManytoOne'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    ####################################################################################################################

    # READ AND PROCESS STK DATA REPORTS
    ####################################################################################################################
    # Initialize lists
    data = {}
    data['total_active_time'] = []
    data['total_blackout_time'] = []
    data['max_active_time'] = []
    data['mean_active_time'] = []
    data['max_blackout_time'] = []
    data['mean_blackout_time'] = []
    data['mean_range'] = []
    data['mean_max_range'] = []
    data['mean_min_range'] = []
    data['total_stored_power_time'] = []
    data['mean_stored_power_time'] = []
    data['max_stored_power_time'] = []
    data['total_station_keeping'] = []
    data['min_active_time'] = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Cycle through available orbit configurations and calculate active/blackout durations
    # Comment out this loop if reading processed data in from a txt file
    for i in range(1, len(orbit_data)):
        print('Progress: {}%'.format(round(100.0 * (i - 1) / (len(orbit_data) - 2), 2)))
        print("Perigee radius: {} km, Apogee radius: {} km".format(orbit_data[i][0], orbit_data[i][1]))

        # Import SPS illumination and access events
        try:
            sps_lighting_one = parse_csv_to_array('{}\DVP_{}_{}perigee{}apogee_0.0argperi_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        except IOError:
            print('{}\DVP_{}_{}perigee{}apogee_00.0argperi_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]))
            write_data_files(stk_data_path, study_name, data)
            sys.exit()

        try:
            sps_access_one = parse_csv_to_array('{}\DVP_{}_{}perigee{}apogee_0.0argperi_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        except IOError:
            print('{}\DVP_{}_{}perigee{}apogee_0.0argperi_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]))
            write_data_files(stk_data_path, study_name, data)
            sys.exit()

        try:
            sps_lighting_two = parse_csv_to_array('{}\DVP_{}_{}perigee{}apogee_180.0argperi_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        except IOError:
            print('{}\DVP_{}_{}perigee{}apogee_180.0argeperi_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]))
            write_data_files(stk_data_path, study_name, data)
            sys.exit()

        try:
            sps_access_two = parse_csv_to_array('{}\DVP_{}_{}perigee{}apogee_180.0argperi_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        except IOError:
            print('{}\DVP_{}_{}perigee{}apogee_180.0argeperi_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]))
            write_data_files(stk_data_path, study_name, data)
            sys.exit()

        # If no access periods exist, insert nan into each list - infeasible orbit design
        if not hasattr(sps_access_one[0], "__len__") or not hasattr(sps_access_two[0], "__len__"):
            data['total_active_time'].append(np.nan)
            data['max_active_time'].append(np.nan)
            data['mean_active_time'].append(np.nan)

            # Blackout equivalent to total eclipse in this case
            data['total_blackout_time'].append(np.sum(target_eclipse[2]))
            data['max_blackout_time'].append(max(target_eclipse[2]))
            data['mean_blackout_time'].append(np.mean(target_eclipse[2]))

            data['mean_range'].append(np.nan)
            data['mean_min_range'].append(np.nan)
            data['mean_max_range'].append(np.nan)

            data['total_stored_power_time'].append(np.nan)
            data['max_stored_power_time'].append(np.nan)
            data['mean_stored_power_time'].append(np.nan)

            data['total_station_keeping'].append(np.nan)
            data['min_active_time'].append(np.nan)

        # Otherwise determine active periods
        else:
            # Calculate individual and combined active durations
            sps_active_one = determine_SPS_active_time(sps_lighting_one, target_eclipse, sps_access_one)
            sps_active_two = determine_SPS_active_time(sps_lighting_two, target_eclipse, sps_access_two)
            sps_active = combine_events(sps_active_one, sps_active_two)

            # IF no active period, put nan
            if np.sum(sps_active) == 0.0:
                data['total_active_time'].append(np.nan)
                data['max_active_time'].append(np.nan)
                data['mean_active_time'].append(np.nan)
                data['mean_range'].append(np.nan)
                data['mean_min_range'].append(np.nan)
                data['mean_max_range'].append(np.nan)
                # Blackout equivalent to total eclipse in this case
                data['total_blackout_time'].append(np.sum(target_eclipse[2]))
                data['max_blackout_time'].append(max(target_eclipse[2]))
                data['mean_blackout_time'].append(np.mean(target_eclipse[2]))
                data['total_station_keeping'].append(np.nan)
                data['min_active_time_2sps'].append(np.nan)
                data['min_active_time_1sps'].append(np.nan)

            # Otherwise calculate blackout period
            else:
                data['total_active_time'].append(np.sum(sps_active[2]))
                data['max_active_time'].append(max(sps_active[2]))
                data['mean_active_time'].append(np.mean(sps_active[2]))

                # Determine the total and maximum target blackout durations
                target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
                data['total_blackout_time'].append(np.sum(target_blackout[2]))
                data['max_blackout_time'].append(max(target_blackout[2]))
                data['mean_blackout_time'].append(np.mean(target_blackout[2]))
                data['min_active_time'].append(min(sps_active[2]))

                # Determine the mean range, averaging mean range per access event in time
                try:
                    sps_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_180.0argperi_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
                except IOError:
                    print('{}\DVP_{}_{}perigee{}apogee_180.0meananom_range.txt'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]))
                    write_data_files(stk_data_path, study_name, data)
                    sys.exit()
                data['mean_range'].append(np.sum([(i * j) / np.sum(sps_access_two[2]) for i, j in zip(sps_range[2], sps_access_two[2])]))
                data['mean_min_range'].append(np.sum([(i * j) / np.sum(sps_access_two[2]) for i, j in zip(sps_range[0], sps_access_two[2])]))
                data['mean_max_range'].append( np.sum([(i * j) / np.sum(sps_access_two[2]) for i, j in zip(sps_range[1], sps_access_two[2])]))

                # Determine station-keeping/battery-charging events
                sps_station_keeping_events_two = determine_battery_chargeup_events(sps_lighting_two, sps_access_two, total_duration)
                data['total_station_keeping'].append(np.sum(sps_station_keeping_events_two[2]))

            # Determine events when SPS is in range but eclipsed
            sps_eclipse = invert_events_list(sps_lighting_two, total_duration)
            sps_use_stored_power = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access_two)
            if np.sum(sps_use_stored_power) == 0.0:
                data['total_stored_power_time'].append(np.nan)
                data['max_stored_power_time'].append(np.nan)
                data['mean_stored_power_time'].append(np.nan)
            else:
                data['total_stored_power_time'].append(np.sum(sps_use_stored_power[2]))
                data['max_stored_power_time'].append(max(sps_use_stored_power[2]))
                data['mean_stored_power_time'].append(np.mean(sps_use_stored_power[2]))

    write_data_files(stk_data_path, study_name, data)
    ####################################################################################################################


def write_data_files(stk_data_path, study_name, data):
    # WRITE PROCESSED DATA TO FILE
    ####################################################################################################################
    # Write data to a file so that the file can be read as opposed to importing and processing
    # data every time. This can help speed up debugging/analysis. But first, remove the file
    # if it already exists and create a new one.

    # ACTIVE TIME
    write_data_to_file(stk_data_path, study_name, data['total_active_time'], 'TotalActive_2SPS')
    write_data_to_file(stk_data_path, study_name, data['max_active_time'], 'MaxActive_2SPS')
    write_data_to_file(stk_data_path, study_name, data['mean_active_time'], 'MeanActive_2SPS')
    write_data_to_file(stk_data_path, study_name, data['min_active_time'], 'MinActive_2SPS')

    # BLACKOUT TIME
    write_data_to_file(stk_data_path, study_name, data['total_blackout_time'], 'TotalBlackout_2SPS')
    write_data_to_file(stk_data_path, study_name, data['max_blackout_time'], 'MaxBlackout_2SPS')
    write_data_to_file(stk_data_path, study_name, data['mean_blackout_time'], 'MeanBlackout_2SPS')

    # RANGE
    write_data_to_file(stk_data_path, study_name, data['mean_range'], "MeanRange_180.0argperi")
    write_data_to_file(stk_data_path, study_name, data['mean_min_range'], "MeanMinRange_180.0argperi")
    write_data_to_file(stk_data_path, study_name, data['mean_max_range'], "MeanMaxRange_180.0argperi")

    # STORED POWER EVENTS
    write_data_to_file(stk_data_path, study_name, data['total_stored_power_time'], "TotalStoredPowerEvent_180.0argperi")
    write_data_to_file(stk_data_path, study_name, data['mean_stored_power_time'], "MeanStoredPowerEvent_180.0argperi")
    write_data_to_file(stk_data_path, study_name, data['max_stored_power_time'], "MaxStoredPowerEvent_180.0argperi")

    # STATION KEEPING EVENTS
    write_data_to_file(stk_data_path, study_name, data['total_station_keeping'], "TotalStationKeeping_180.0argperi")
    ####################################################################################################################


main()
