""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data reports obtained programmatically from STK for evaluating
the optimal SPS orbit for the equatorial lunar SPS configuration.
"""

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def main():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Difference in subsequent apogee/perigee radii
    r_moon = 1737.0
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name of study
    study_name = 'Equatorial_IncrementedRes'

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Get orbital data
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # READ AND PROCESS DATA REPORTS
    ####################################################################################################################
    # Initialize lists
    total_active_time = []
    total_blackout_time = []
    max_active_time = []
    mean_active_time = []
    std_active_time = []
    min_active_time = []
    max_blackout_time = []
    min_blackout_time = []
    mean_blackout_time = []
    std_blackout_time = []
    mean_range = []
    mean_max_range = []
    mean_min_range = []
    total_stored_power_time = []
    mean_stored_power_time = []
    max_stored_power_time = []
    total_station_keeping = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Cycle through available orbit configurations and calculate active/blackout durations
    for i in range(1, len(orbit_data)):
        print('Progress: {}%'.format(round(100.0 * (i-1) / (len(orbit_data) - 2), 2)))
        print("Perigee altitude: {} km, Apogee altitude: {} km".format(orbit_data[i][0] - r_moon, orbit_data[i][1] - r_moon))

        # Import access and lighting for SPS
        sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)

        # If no access periods exist, insert nan into each list - infeasible orbit design
        if not hasattr(sps_access[0], "__len__"):
            total_active_time.append(np.nan)
            max_active_time.append(np.nan)
            mean_active_time.append(np.nan)
            min_active_time.append(np.nan)
            std_active_time.append(np.nan)

            total_blackout_time.append(np.sum(target_eclipse[2]))
            max_blackout_time.append(max(target_eclipse[2]))
            mean_blackout_time.append(np.mean(target_eclipse[2]))
            min_blackout_time.append(min(target_eclipse[2]))
            std_blackout_time.append(np.std(target_eclipse[2]))

            mean_range.append(np.nan)
            mean_min_range.append(np.nan)
            mean_max_range.append(np.nan)

            total_stored_power_time.append(np.nan)
            max_stored_power_time.append(np.nan)
            mean_stored_power_time.append(np.nan)

            total_station_keeping.append(np.nan)

        # Otherwise insert results into list
        else:
            # Calculate active times
            sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
            # Check if there is active times, if not insert nan
            if np.sum(sps_active) == 0.0:
                # Active times
                total_active_time.append(np.nan)
                max_active_time.append(np.nan)
                mean_active_time.append(np.nan)
                min_active_time.append(np.nan)
                std_active_time.append(np.nan)

                # Range
                mean_range.append(np.nan)
                mean_min_range.append(np.nan)
                mean_max_range.append(np.nan)

                # Blackout times
                total_blackout_time.append(np.sum(target_eclipse[2]))
                max_blackout_time.append(max(target_eclipse[2]))
                mean_blackout_time.append(np.mean(target_eclipse[2]))
                min_blackout_time.append(min(target_eclipse[2]))
                std_blackout_time.append(np.std(target_eclipse[2]))

                # Station-keeping/battery-charge times
                total_station_keeping.append(np.nan)

            # If active times, calculate blackouts and range data the insert statistics
            else:
                # Active times
                total_active_time.append(np.sum(sps_active[2]))
                max_active_time.append(max(sps_active[2]))
                mean_active_time.append(np.mean(sps_active[2]))
                min_active_time.append(min(sps_active[2]))
                std_active_time.append(np.std(sps_active[2]))

                # Range
                sps_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
                assert (len(sps_range[0]) == len(sps_access[0]))
                mean_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[2], sps_access[2])]))
                mean_min_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[0], sps_access[2])]))
                mean_max_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[1], sps_access[2])]))

                # Blackout times
                target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
                total_blackout_time.append(np.sum(target_blackout[2]))
                max_blackout_time.append(max(target_blackout[2]))
                mean_blackout_time.append(np.mean(target_blackout[2]))
                min_blackout_time.append(min(target_blackout[2]))
                std_blackout_time.append(np.std(target_blackout[2]))

                # Determine station-keeping/battery-charging events
                sps_station_keeping_events = determine_battery_chargeup_events(sps_lighting, sps_access, total_duration)
                total_station_keeping.append(np.sum(sps_station_keeping_events[2]))

            # Stored power events
            sps_eclipse = invert_events_list(sps_lighting, total_duration)
            sps_use_stored_power = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)
            # Check if stored power events exist, if not insert nan
            if np.sum(sps_use_stored_power) == 0.0:
                total_stored_power_time.append(np.nan)
                max_stored_power_time.append(np.nan)
                mean_stored_power_time.append(np.nan)
            else:
                total_stored_power_time.append(np.sum(sps_use_stored_power[2]))
                max_stored_power_time.append(max(sps_use_stored_power[2]))
                mean_stored_power_time.append(np.mean(sps_use_stored_power[2]))
    ####################################################################################################################

    # WRITING DATA TO FILES
    ####################################################################################################################
    # Writing processed data to files to save time in future analysis. Function will overwrite old data files.

    # MEAN RANGE
    write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange")
    write_data_to_file(stk_data_path, study_name, mean_min_range, "MeanMinRange")
    write_data_to_file(stk_data_path, study_name, mean_max_range, "MeanMaxRange")

    # ACTIVE TIME DATA
    write_data_to_file(stk_data_path, study_name, total_active_time, "TotalActive")
    write_data_to_file(stk_data_path, study_name, max_active_time, "MaxActive")
    write_data_to_file(stk_data_path, study_name, mean_active_time, "MeanActive")
    write_data_to_file(stk_data_path, study_name, min_active_time, "MinActive")
    write_data_to_file(stk_data_path, study_name, std_active_time, "StdActive")

    # BLACKOUT TIME DATA
    write_data_to_file(stk_data_path, study_name, total_blackout_time, "TotalBlackout")
    write_data_to_file(stk_data_path, study_name, max_blackout_time, "MaxBlackout")
    write_data_to_file(stk_data_path, study_name, mean_blackout_time, "MeanBlackout")
    write_data_to_file(stk_data_path, study_name, min_blackout_time, "MinBlackout")
    write_data_to_file(stk_data_path, study_name, std_blackout_time, "StdBlackout")

    # STORED POWER EVENTS
    write_data_to_file(stk_data_path, study_name, total_stored_power_time, "TotalStoredPowerEvent")
    write_data_to_file(stk_data_path, study_name, mean_stored_power_time, "MeanStoredPowerEvent")
    write_data_to_file(stk_data_path, study_name, max_stored_power_time, "MaxStoredPowerEvent")

    # STATION KEEPING EVENTS
    write_data_to_file(stk_data_path, study_name, total_station_keeping, "TotalStationKeeping")
    ####################################################################################################################


main()
