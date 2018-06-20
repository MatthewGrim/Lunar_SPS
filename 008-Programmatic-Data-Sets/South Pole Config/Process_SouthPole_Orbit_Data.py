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
    study_name = 'SouthPole_IncrementedRes_Inertial'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    ####################################################################################################################

    # READ AND PROCESS STK DATA REPORTS
    ####################################################################################################################
    # Initialize lists
    total_active_time = []
    total_blackout_time = []
    max_active_time = []
    mean_active_time = []
    max_blackout_time = []
    mean_blackout_time = []
    mean_range = []

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
        sps_lighting_raw = '{}\DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
        sps_access_raw = '{}\DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        sps_access = parse_csv_to_array(sps_access_raw, start)

        # Determine the total and maximum SPS active durations
        sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
        total_active_time.append(np.sum(sps_active[2]))
        max_active_time.append(max(sps_active[2]))
        mean_active_time.append(np.mean(sps_active[2]))

        # Determine the total and maximum target blackout durations
        target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        total_blackout_time.append(np.sum(target_blackout[2]))
        max_blackout_time.append(max(target_blackout[2]))
        mean_blackout_time.append(np.mean(target_blackout[2]))

        # Determine the mean range, averaging mean range per access event in time
        sps_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
        mean_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[2], sps_access[2])]))
    ####################################################################################################################

    # WRITE PROCESSED DATA TO FILE
    ####################################################################################################################
    # Write data to a file so that the file can be read as opposed to importing and processing
    # data every time. This can help speed up debugging/analysis. But first, remove the file
    # if it already exists and create a new one.

    # TOTAL ACTIVE TIME
    write_data_to_file(stk_data_path, study_name, total_active_time, 'TotalActive_Inertial_Extended')
    # TOTAL BLACKOUT TIME
    write_data_to_file(stk_data_path, study_name, total_blackout_time, 'TotalBlackout_Inertial_Extended')
    # MAX ACTIVE TIME
    write_data_to_file(stk_data_path, study_name, max_active_time, 'MaxActive_Inertial_Extended')
    # MAX BLACKOUT TIME
    write_data_to_file(stk_data_path, study_name, max_blackout_time, 'MaxBlackout_Inertial_Extended')
    # MEAN ACTIVE TIME
    write_data_to_file(stk_data_path, study_name, mean_active_time, 'MeanActive_Inertial_Extended')
    # MEAN BLACKOUT TIME
    write_data_to_file(stk_data_path, study_name, mean_blackout_time, 'MeanBlackout_Inertial_Extended')
    # MEAN RANGE
    write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange_Inertial_Extended")
    ####################################################################################################################


main()
