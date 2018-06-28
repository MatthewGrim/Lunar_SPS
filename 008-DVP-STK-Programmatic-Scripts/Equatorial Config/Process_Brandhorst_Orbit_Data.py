""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the Brandhorst lunar SPS configuration.
"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *
import os


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
    max_active_time = []
    mean_active_time = []
    total_blackout_time = []
    max_blackout_time = []
    mean_blackout_time = []
    mean_range = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Cycle through available orbit configurations and calculate active/blackout durations
    # Comment out if reading processed data in from a txt file
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
        sps_range = import_range_data('{}/DVP_{}_{}perigee_{}apogee_range.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
        mean_range.append(np.mean(sps_range[1]))
    ####################################################################################################################

    # WRITING DATA TO FILES
    ####################################################################################################################
    # Writing processed data to files to save time in future analysis. Comment this out if the processed data files
    # already exist and are being read out (see next section). Function will overwrite old data files.

    # MEAN RANGE
    write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange_Equatorial")
    # ACTIVE TIME DATA
    write_data_to_file(stk_data_path, study_name, total_active_time, "TotalActive_Equatorial")
    write_data_to_file(stk_data_path, study_name, max_active_time, "MaxActive_Equatorial")
    write_data_to_file(stk_data_path, study_name, mean_active_time, "MeanActive_Equatorial")
    # BLACKOUT TIME DATA
    write_data_to_file(stk_data_path, study_name, total_blackout_time, "TotalBlackout_Equatorial")
    write_data_to_file(stk_data_path, study_name, max_blackout_time, "MaxBlackout_Equatorial")
    write_data_to_file(stk_data_path, study_name, mean_blackout_time, "MeanBlackout_Equatorial")
    ####################################################################################################################

    # READ IN DATA FILES IF THEY ARE ALREADY WRITTEN
    ####################################################################################################################
    # Read in saved active time data, if the data is already processed and saved into a txt file
    # Comment this section out if the data is being processed, and the results are not yet written to a txt file

    # MEAN RANGE
    # mean_range = read_data_from_file(stk_data_path, study_name, "MeanRange_Equatorial")
    # # ACTIVE TIME
    # total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive_Equatorial")
    # max_active_time = read_data_from_file(stk_data_path, study_name, "MaxActive_Equatorial")
    # mean_active_time = read_data_from_file(stk_data_path, study_name, "MeanActive_Equatorial")
    # # BLACKOUT TIME
    # total_blackout_time = read_data_from_file(stk_data_path, study_name, "TotalBlackout_Equatorial")
    # max_blackout_time = read_data_from_file(stk_data_path, study_name, "MaxBlackout_Equatorial")
    # mean_blackout_time = read_data_from_file(stk_data_path, study_name, "MeanBlackout_Equatorial")

    ####################################################################################################################

    # APPLY POINTING AND POWER CONSTRAINTS
    ####################################################################################################################
    # ENFORCING POINTING AND POWER CONSTRAINTS
    ####################################################################################################################
    # # Transmitter parameters
    # # IPG YLS10000 Industrial HP Laser
    # wavelength = 1070e-9
    # trans_radius = 0.3
    # trans_power = 100e3
    # trans_mass = 3600.0
    # trans_eff = 0.35
    # IPG YLS-CUT
    # wavelength = 1070e-9
    # trans_radius = 0.1
    # trans_power = 15e3
    # trans_mass = 440.0
    # trans_eff = 0.35
    # # Receiver parameters for AMALIA rover receiver
    # rec_radius = 0.5
    # rec_efficiency = 0.40
    # # Pointing error in radians
    # point_error = 1e-6
    # # Minimum power received at target
    # min_power = 100.0
    #
    # # Initialize lists
    # mean_link_efficiency = []
    # mean_power_received = []
    #
    # for j in mean_range:
    #
    #     # Minimum beam radius as defined by pointing error
    #     min_beam_radius = rec_radius + (point_error * j * 1000.0)
    #
    #     # Actual beam radius as defined by Gaussian beam divergence
    #     surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (j * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
    #
    #     # Calculate mean link efficiency
    #     link_eff_temp = (rec_radius / surf_beam_radius) ** 2
    #
    #     # Check pointing error constraint
    #     if surf_beam_radius < min_beam_radius:
    #         link_eff_temp = np.nan
    #
    #     mean_pwr_temp = rec_efficiency * trans_power * link_eff_temp
    #
    #     # Check power delivery constraint
    #     if mean_pwr_temp < min_power:
    #         mean_pwr_temp = np.nan
    #
    #     # Get mean link efficiency
    #     mean_link_efficiency.append(link_eff_temp)
    #     # Get mean power received
    #     mean_power_received.append(mean_pwr_temp)

    # Apply effect of constraints to remaining data set
    # for i in range(len(mean_power_received)):
    #     if math.isnan(mean_power_received[i]):
    #         total_active_time[i] = np.nan
    #         total_blackout_time[i] = np.nan
    #
    #         max_active_time[i] = np.nan
    #         max_blackout_time[i] = np.nan
    #
    #         mean_active_time[i] = np.nan
    #         mean_blackout_time[i] = np.nan

    # Calculate total energy delivered to receiver
    # total_energy = [i * j for i, j in zip(mean_power_received, total_active_time)]
    ####################################################################################################################

    # ENFORCING BLACKOUT DURATION CONSTRAINTS
    ####################################################################################################################
    # Determine the maximum duration for which a rover could hibernate without completely draining the battery
    # The following is for the AMALIA rover
    # Rover battery capacity in Watt-hours
    # battery_capacity = 100
    # Minimum power consumption, for hibernation mode
    # hibernation_pwr = 7
    # Maximum time rover can survive without recharging (in seconds)
    # rover_blackout_limit = (battery_capacity / hibernation_pwr) * 3600

    # # Remove data points for which blackout durations exceed the limit
    # for i in range(len(max_blackout_time)):
    #     if max_blackout_time[i] > rover_blackout_limit:
    #         max_blackout_time[i] = np.nan
    #
    #         # Apply constraint to all other data sets
    #         total_energy[i] = np.nan
    #         max_active_time[i] = np.nan
    #
    #         total_active_time[i] = np.nan
    #         total_blackout_time[i] = np.nan
    #
    #         mean_active_time[i] = np.nan
    #         mean_blackout_time[i] = np.nan
    #
    #         mean_link_efficiency[i] = np.nan
    #         mean_power_received[i] = np.nan

    # EVALUATE REQUIRED POWER GENERATOR METRICS
    ####################################################################################################################
    # Power generator parameter - Stretched lens array SquareRigger platform
    # generator_eff = 0.4
    # generator_spec_pwr = 300
    # # Power required from generator
    # generator_pwr = trans_power / (generator_eff * trans_eff)
    # # Mass of generator
    # generator_mass = generator_pwr / generator_spec_pwr
    ####################################################################################################################

    # REORGANIZE DATA FOR PLOTTING
    ####################################################################################################################
    total_active_times_sorted = sort_data_list_into_array(orbit_data, resolution, [i / 3600.0 for i in total_active_time])
    max_blackout_times_sorted = sort_data_list_into_array(orbit_data, resolution, [i / 3600.0 for i in max_blackout_time])
    total_blackout_times_sorted = sort_data_list_into_array(orbit_data, resolution, [i / 3600.0 for i in total_blackout_time])
    link_efficiency_sorted = sort_data_list_into_array(orbit_data, resolution, mean_link_efficiency)
    mean_power_sorted = sort_data_list_into_array(orbit_data, resolution, mean_power_received)
    total_energy_sorted = sort_data_list_into_array(orbit_data, resolution, total_energy)

    # Get list of unique perigees and apogees tested
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])
    apogee_altitudes = [i - r_moon for i in unique_apogees]
    perigee_altitudes = [i - r_moon for i in unique_perigees]

    # Find orbit which results in highest performance according to objective function
    # from numpy import unravel_index
    # best_orbit_idx = unravel_index(objective.argmax(), objective.shape)
    # best_perigee = unique_perigees[best_orbit_idx[0]]
    # best_apogee = unique_apogees[best_orbit_idx[1]]
    # print('Best orbit --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee, 2), round(best_apogee, 2)))
    # print('Total active time --> {} hours, or {}%'.format(round(active_times[best_orbit_idx], 2), round(100.0 * active_times[best_orbit_idx] * 3600.0 / total_duration, 2)))
    # print('Link efficiency --> {}%'.format(round(link_efficiency[best_orbit_idx] * 100.0, 5)))
    # print('Power received --> {} W'.format(round(power_received[best_orbit_idx], 2)))
    ####################################################################################################################

    # PLOTTING THE ACTIVE TIMES, LINK EFFICIENCIES, AND OBJECTIVES
    ####################################################################################################################
    # Generate figure showing access time as a function of perigee/apogee
    plt.figure(1)
    plt.contourf(apogee_altitudes, perigee_altitudes, total_active_times_sorted, 500)
    plt.colorbar()
    plt.title('Total Active Time [hrs]')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()

    # Generate figure showing mean link efficiency as a function of perigee/apogee
    plt.figure(2)
    plt.contourf(apogee_altitudes, perigee_altitudes, mean_power_sorted, 500)
    plt.colorbar()
    plt.title('Mean Power Received [W]')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()

    # Generate figure showing mean link efficiency as a function of perigee/apogee
    plt.figure(3)
    plt.contourf(apogee_altitudes, perigee_altitudes, link_efficiency_sorted * 100.0, 500)
    plt.colorbar()
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()

    # Generate figure showing mean link efficiency as a function of perigee/apogee
    plt.figure(3)
    plt.contourf(apogee_altitudes, perigee_altitudes, total_blackout_times_sorted, 500)
    plt.colorbar()
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()

    # Generate figure showing objetive function (combining total access and link efficiency)
    # Mark orbit which results in highest performance
    # plt.figure(4)
    # plt.contourf(unique_apogees, unique_perigees, objective, 1000)
    # plt.colorbar()
    # plt.plot(best_apogee, best_perigee, marker='x', markersize=5, color='green')
    # plt.title('Sum of Mean Link Efficiency and Total Active Time')
    # plt.xlabel('Apogee (Dist to Focus) [km]')
    # plt.ylabel('Perigee (Dist to Focus) [km]')
    # plt.show()


main()
