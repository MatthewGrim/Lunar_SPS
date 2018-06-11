""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the polar lunar SPS configuration.
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
    study_name = 'SouthPole_IncrementedRes'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    ####################################################################################################################

    # PARAMETRIC SCAN
    ###################################################################################################################
    # Import report data generated via STK Connect Commands
    total_active_time = []
    max_active_duration = []
    mean_link_efficiency = []
    mean_power_received = []
    mean_range = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)

    # Cycle through available orbit configurations and calculate active/blackout durations
    # Comment out if reading processed data in from a txt file
    # for i in range(1, len(orbit_data)):
    #     print('Progress: {}%'.format(round(100.0 * (i - 1) / (len(orbit_data) - 2), 2)))
    #     print("Perigee radius: {} km, Apogee radius: {} km".format(orbit_data[i][0], orbit_data[i][1]))
    #
    #     # CALCULATING ACTIVE TIMES
    #     ################################################################################################################
    #     # Import SPS illumination and access events
    #     sps_lighting_raw = '{}\DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
    #     sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    #     sps_access_raw = '{}\DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
    #     sps_access = parse_csv_to_array(sps_access_raw, start)
    #
    #     # Determine the total and maximum SPS active durations
    #     sps_active = determine_SPS_active_time(sps_lighting, sps_access, target_lighting)
    #     total_active_time.append(np.sum(sps_active[2]))
    #     max_active_duration.append(max(sps_active[2]))
    #
    #     # CALCULATING MEAN RANGE
    #     ################################################################################################################
    #     # Import SPS range data
    #     # First entry to sps_range is minimum range for each access event
    #     # Second entry is maximum range for each access event
    #     # Third entry is mean range for each access event
    #     sps_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
    #     mean_range.append(np.mean(sps_range[2]))
    #
    ####################################################################################################################

    # WRITE PROCESSED DATA TO FILE
    ################################################################################################################
    # Write data to a file so that the file can be read as opposed to importing and processing
    # data every time. This can help speed up debugging/analysis. But first, remove the file
    # if it already exists and create a new one.

    # # TOTAL ACTIVE TIME
    # write_data_to_file(stk_data_path, study_name, total_active_time, 'TotalActive')
    # # MEAN RANGE
    # write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange")

    ####################################################################################################################

    # READ IN DATA FILES IF THEY EXIST
    ####################################################################################################################
    # Read in saved active time data, if the data is already processed and saved into a txt file
    # Comment this section out if the data is being processed, and the results are not yet written to a txt file

    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive")

    mean_range = read_data_from_file(stk_data_path, study_name, "MeanRange")
    ####################################################################################################################

    # ENFORCING POINTING AND POWER CONSTRAINTS
    ################################################################################################################
    # Transmitter parameters
    wavelength = 1070e-9
    trans_radius = 0.3
    trans_power = 100e3
    # Receiver parameters for AMALIA rover receiver
    rec_radius = 0.5
    rec_efficiency = 0.35
    # Pointing error in radians
    point_error = 1e-6
    # Minimum power received by rover
    min_power = 300.0

    for j in mean_range:
        min_beam_radius = rec_radius + (point_error * j * 1000.0)
        surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (j * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)

        # Calculate mean link efficiency
        link_eff_temp = (rec_radius / surf_beam_radius) ** 2

        # Check pointing error constraint
        if surf_beam_radius < min_beam_radius:
            link_eff_temp = np.nan

        mean_pwr_temp = rec_efficiency * trans_power * link_eff_temp

        # Check power delivery constraint
        if mean_pwr_temp < min_power:
            mean_pwr_temp = np.nan

        # Get mean link efficiency
        mean_link_efficiency.append(link_eff_temp)
        # Get mean power received
        mean_power_received.append(mean_pwr_temp)

    # APPLY CONSTRAINTS TO REMAINING DATA SETS
    ####################################################################################################################

    for i in range(len(mean_power_received)):
        if math.isnan(mean_power_received[i]):
            total_active_time[i] = np.nan

    ####################################################################################################################

    # REORGANIZE DATA FOR PLOTTING
    ####################################################################################################################

    active_times_sorted = sort_data_list_with_incremented_resolution_into_array(orbit_data, total_active_time)

    link_efficiency_sorted = sort_data_list_with_incremented_resolution_into_array(orbit_data, mean_link_efficiency)

    power_received_sorted = sort_data_list_with_incremented_resolution_into_array(orbit_data, mean_power_received)

    # PLOTTING THE ACTIVE TIMES, LINK EFFICIENCIES, AND OBJECTIVES
    ####################################################################################################################
    # Find unique perigees and apogees tested for plotting
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    r_moon = 1737.0
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    # Find orbit which results in highest performance according to objective function
    from numpy import unravel_index
    best_orbit_idx = unravel_index(np.nanargmax(active_times_sorted), active_times_sorted.shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    print('Best orbit --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Total active time --> {} hours, or {}%'.format(round(active_times_sorted[best_orbit_idx] / 3600.0, 2), round(
        100.0 * active_times_sorted[best_orbit_idx]/ total_duration, 2)))
    print('Link efficiency --> {}%'.format(round(link_efficiency_sorted[best_orbit_idx] * 100.0, 5)))
    print('Power received --> {} W'.format(round(power_received_sorted[best_orbit_idx], 2)))

    # Generate figures
    perigee_altitudes = [i - r_moon for i in unique_perigees]
    apogee_altitudes = [i - r_moon for i in unique_apogees]
    make_contour_plot(perigee_altitudes, apogee_altitudes, power_received_sorted, "Mean Power Received [W]", 1)


main()
