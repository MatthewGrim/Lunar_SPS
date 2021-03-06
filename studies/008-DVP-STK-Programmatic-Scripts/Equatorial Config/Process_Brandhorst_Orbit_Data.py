""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the Brandhorst lunar SPS configuration.
"""

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def main():
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Difference in subsequent apogee/perigee radii
    r_moon = 1737.0
    resolution = 1000.0
    max_perigee = 10000.0
    min_perigee = 1000.0
    max_apogee = 54000.0

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name of study
    study_name = 'Brandhorst_{}kmRes'.format(resolution)

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Get orbital data
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements(resolution, min_perigee, max_perigee, max_apogee)

    # READ AND PROCESS DATA REPORTS
    ####################################################################################################################
    # Initialize lists
    total_active_time = []
    total_blackout_time = []
    max_active_time = []
    mean_active_time = []
    max_blackout_time = []
    mean_blackout_time = []
    mean_range = []
    mean_max_range = []
    mean_min_range = []
    total_stored_power_time = []
    mean_stored_power_time = []
    max_stored_power_time = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Cycle through available orbit configurations and calculate active/blackout durations
    for i in range(1, len(orbit_data)):
        print('Progress: {}%'.format(round(100.0 * (i-1) / (len(orbit_data) - 2), 2)))
        print("Perigee altitude: {} km, Apogee altitude: {} km".format(orbit_data[i][0] - r_moon, orbit_data[i][1] - r_moon))

        # Import access and lighting for SPS
        sps_lighting_raw = '{}/DVP_{}_{}perigee_{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
        sps_access_raw = '{}/DVP_{}_{}perigee_{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        sps_access = parse_csv_to_array(sps_access_raw, start)

        # Determine the total SPS active durations
        sps_active = determine_SPS_active_time(sps_lighting, sps_access, target_lighting)
        total_active_time.append(np.sum(sps_active[2]))
        max_active_time.append(max(sps_active[2]))
        mean_active_time.append(np.mean(sps_active[2]))

        # Determine the total target blackout durations
        target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        total_blackout_time.append(np.sum(target_blackout[2]))
        max_blackout_time.append(max(target_blackout[2]))
        mean_blackout_time.append(np.mean(target_blackout[2]))

        # Import range statistics and store mean range
        sps_range = import_range_data_statistics('{}/DVP_{}_{}perigee_{}apogee_range.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        mean_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[2], sps_access[2])]))
        mean_min_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[0], sps_access[2])]))
        mean_max_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[1], sps_access[2])]))

        # Determine events when SPS is in range but eclipsed
        sps_eclipse = invert_events_list(sps_lighting, total_duration)
        sps_use_stored_power = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)
        total_stored_power_time.append(np.sum(sps_use_stored_power[2]))
        max_stored_power_time.append(max(sps_use_stored_power[2]))
        mean_stored_power_time.append(np.mean(sps_use_stored_power[2]))
    ####################################################################################################################

    # WRITING DATA TO FILES
    ####################################################################################################################
    # Writing processed data to files to save time in future analysis. Comment this out if the processed data files
    # already exist and are being read out (see next section). Function will overwrite old data files.

    # MEAN RANGE
    write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange_Brandhorst")
    write_data_to_file(stk_data_path, study_name, mean_min_range, "MeanMinRange_Brandhorst")
    write_data_to_file(stk_data_path, study_name, mean_max_range, "MeanMaxRange_Brandhorst")

    # ACTIVE TIME DATA
    write_data_to_file(stk_data_path, study_name, total_active_time, "TotalActive_Brandhorst")
    write_data_to_file(stk_data_path, study_name, max_active_time, "MaxActive_Brandhorst")
    write_data_to_file(stk_data_path, study_name, mean_active_time, "MeanActive_Brandhorst")

    # BLACKOUT TIME DATA
    write_data_to_file(stk_data_path, study_name, total_blackout_time, "TotalBlackout_Brandhorst")
    write_data_to_file(stk_data_path, study_name, max_blackout_time, "MaxBlackout_Brandhorst")
    write_data_to_file(stk_data_path, study_name, mean_blackout_time, "MeanBlackout_Brandhorst")

    # STORED POWER EVENTS
    write_data_to_file(stk_data_path, study_name, total_stored_power_time, "TotalStoredPowerEvent_Brandhorst")
    write_data_to_file(stk_data_path, study_name, mean_stored_power_time, "MeanStoredPowerEvent_Brandhorst")
    write_data_to_file(stk_data_path, study_name, max_stored_power_time, "MaxStoredPowerEvent_Brandhorst")
    ####################################################################################################################


main()
