""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the polar lunar SPS configuration based on feasible pointing and power
constraints.

There is a section for reading and processing STK data reports, and wri ting the processed
data to a file. The following section is for reading the processed data files if they have
already been written. Comment out one section or the other as necessary.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import *



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
    min_active_time = []
    max_blackout_time = []
    mean_blackout_time = []
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
        min_active_time.append(min(sps_active[2]))
        
        # Determine station-keeping/battery-charging events
        sps_station_keeping_events = determine_battery_chargeup_events(sps_lighting, sps_access, total_duration)
        total_station_keeping.append(np.sum(sps_station_keeping_events[2]))

        # Determine the total and maximum target blackout durations
        target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        total_blackout_time.append(np.sum(target_blackout[2]))
        max_blackout_time.append(max(target_blackout[2]))
        mean_blackout_time.append(np.mean(target_blackout[2]))

        # Determine the mean range, averaging mean range per access event in time
        sps_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
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

    # WRITE PROCESSED DATA TO FILE
    ####################################################################################################################
    # Write data to a file so that the file can be read as opposed to importing and processing
    # data every time. This can help speed up debugging/analysis. But first, remove the file
    # if it already exists and create a new one.

    # ACTIVE TIME
    write_data_to_file(stk_data_path, study_name, total_active_time, 'TotalActive')
    write_data_to_file(stk_data_path, study_name, max_active_time, 'MaxActive')
    write_data_to_file(stk_data_path, study_name, mean_active_time, 'MeanActive')
    write_data_to_file(stk_data_path, study_name, min_active_time, 'MinActive')

    # BLACKOUT TIME
    write_data_to_file(stk_data_path, study_name, total_blackout_time, 'TotalBlackout')
    write_data_to_file(stk_data_path, study_name, max_blackout_time, 'MaxBlackout')
    write_data_to_file(stk_data_path, study_name, mean_blackout_time, 'MeanBlackout')

    # RANGE
    write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange")
    write_data_to_file(stk_data_path, study_name, mean_min_range, "MeanMinRange")
    write_data_to_file(stk_data_path, study_name, mean_max_range, "MeanMaxRange")

    # STORED POWER EVENTS
    write_data_to_file(stk_data_path, study_name, total_stored_power_time, "TotalStoredPowerEvent")
    write_data_to_file(stk_data_path, study_name, mean_stored_power_time, "MeanStoredPowerEvent")
    write_data_to_file(stk_data_path, study_name, max_stored_power_time, "MaxStoredPowerEvent")
    
    # STATION KEEPING EVENTS
    write_data_to_file(stk_data_path, study_name, total_station_keeping, "TotalStationKeeping")
    ####################################################################################################################


main()
