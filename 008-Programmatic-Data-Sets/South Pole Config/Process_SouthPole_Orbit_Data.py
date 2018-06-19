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
    # for i in range(1, len(orbit_data)):
    #     print('Progress: {}%'.format(round(100.0 * (i - 1) / (len(orbit_data) - 2), 2)))
        # print("Perigee radius: {} km, Apogee radius: {} km".format(orbit_data[i][0], orbit_data[i][1]))

        # Import SPS illumination and access events
        # sps_lighting_raw = '{}\DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        # sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
        # sps_access_raw = '{}\DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        # sps_access = parse_csv_to_array(sps_access_raw, start)

        # Determine the total and maximum SPS active durations
        # sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
        # total_active_time.append(np.sum(sps_active[2]))
        # max_active_time.append(max(sps_active[2]))
        # mean_active_time.append(np.mean(sps_active[2]))

        # Determine the total and maximum target blackout durations
        # target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        # total_blackout_time.append(np.sum(target_blackout[2]))
        # max_blackout_time.append(max(target_blackout[2]))
        # mean_blackout_time.append(np.mean(target_blackout[2]))

        # Determine the mean range, averaging mean range per access event in time
        # sps_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
        # mean_range.append(np.sum([(i * j) / np.sum(sps_access[2]) for i, j in zip(sps_range[2], sps_access[2])]))
    ####################################################################################################################

    # WRITE PROCESSED DATA TO FILE
    ####################################################################################################################
    # Write data to a file so that the file can be read as opposed to importing and processing
    # data every time. This can help speed up debugging/analysis. But first, remove the file
    # if it already exists and create a new one.

    # TOTAL ACTIVE TIME
    # write_data_to_file(stk_data_path, study_name, total_active_time, 'TotalActive_Inertial_Extended')
    # TOTAL BLACKOUT TIME
    # write_data_to_file(stk_data_path, study_name, total_blackout_time, 'TotalBlackout_Inertial_Extended')
    # MAX ACTIVE TIME
    # write_data_to_file(stk_data_path, study_name, max_active_time, 'MaxActive_Inertial_Extended')
    # MAX BLACKOUT TIME
    # write_data_to_file(stk_data_path, study_name, max_blackout_time, 'MaxBlackout_Inertial_Extended')
    # MEAN ACTIVE TIME
    # write_data_to_file(stk_data_path, study_name, mean_active_time, 'MeanActive_Inertial_Extended')
    # MEAN BLACKOUT TIME
    # write_data_to_file(stk_data_path, study_name, mean_blackout_time, 'MeanBlackout_Inertial_Extended')
    # MEAN RANGE
    # write_data_to_file(stk_data_path, study_name, mean_range, "MeanRange_Inertial_Extended")
    ####################################################################################################################

    # READ IN DATA FILES IF THEY EXIST
    ####################################################################################################################
    # Read in saved active time data, if the data is already processed and saved into a txt file
    # Comment this section out if the data is being processed, and the results are not yet written to a txt file

    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive_Inertial_Extended")

    total_blackout_time = read_data_from_file(stk_data_path, study_name, "TotalBlackout_Inertial_Extended")

    max_active_time = read_data_from_file(stk_data_path, study_name, "MaxActive_Inertial_Extended")

    max_blackout_time = read_data_from_file(stk_data_path, study_name, "MaxBlackout_Inertial_Extended")

    mean_range = read_data_from_file(stk_data_path, study_name, "MeanRange_Inertial_Extended")

    mean_active_time = read_data_from_file(stk_data_path, study_name, "MeanActive_Inertial_Extended")

    mean_blackout_time = read_data_from_file(stk_data_path, study_name, "MeanBlackout_Inertial_Extended")
    ####################################################################################################################

    # ENFORCING POINTING AND POWER CONSTRAINTS
    ####################################################################################################################
    # Transmitter parameters
    # IPG YLS10000 Industrial HP Laser
    wavelength = 1070e-9
    trans_radius = 0.25
    trans_power = 100e3
    trans_mass = 3600.0
    trans_eff = 0.35
    # IPG YLS-CUT
    # wavelength = 1070e-9
    # trans_radius = 0.1
    # trans_power = 15e3
    # trans_mass = 440.0
    # trans_eff = 0.35
    # Receiver parameters for AMALIA rover receiver
    rec_radius = 0.5
    rec_efficiency = 0.40
    # Pointing error in radians
    point_error = 1e-6
    # Minimum power received at target
    min_power = 300.0

    # Initialize lists
    mean_link_efficiency = []
    mean_power_received = []

    for j in mean_range:

        # Minimum beam radius as defined by pointing error
        min_beam_radius = rec_radius + (point_error * j * 1000.0)

        # Actual beam radius as defined by Gaussian beam divergence
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

    # Apply effect of constraints to remaining data set
    for i in range(len(mean_power_received)):
        if math.isnan(mean_power_received[i]):
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan

            max_active_time[i] = np.nan
            max_blackout_time[i] = np.nan

            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan

    # Calculate total energy delivered to receiver
    total_energy = [i * j for i, j in zip(mean_power_received, total_active_time)]
    ####################################################################################################################

    # ENFORCING BLACKOUT DURATION CONSTRAINTS
    ####################################################################################################################
    # Determine the maximum duration for which a rover could hibernate without completely draining the battery
    # The following is for the AMALIA rover
    # Rover battery capacity in Watt-hours
    battery_capacity = 100
    # Minimum power consumption, for hibernation mode
    hibernation_pwr = 7
    # Maximum time rover can survive without recharging (in seconds)
    rover_blackout_limit = battery_capacity / hibernation_pwr * 3600

    # Remove data points for which blackout durations exceed the limit
    for i in range(len(max_blackout_time)):
        if max_blackout_time[i] > rover_blackout_limit:
            max_blackout_time[i] = np.nan

            # Apply constraint to all other data sets
            total_energy[i] = np.nan
            max_active_time[i] = np.nan

            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan

            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan

            mean_link_efficiency[i] = np.nan
            mean_power_received[i] = np.nan

    # EVALUATE REQUIRED POWER GENERATOR METRICS
    ####################################################################################################################
    # Power generator parameter - Stretched lens array SquareRigger platform
    generator_eff = 0.4
    generator_spec_pwr = 300
    # Power required from generator
    generator_pwr = trans_power / (generator_eff * trans_eff)
    # Mass of generator
    generator_mass = generator_pwr / generator_spec_pwr
    ####################################################################################################################

    # PLOTTING THE ACTIVE TIMES, LINK EFFICIENCIES, AND OBJECTIVES
    ####################################################################################################################
    # Reorganize the data lists into 2D array
    total_active_times_sorted = sort_incremented_resolution_data(orbit_data, total_active_time)
    total_blackout_times_sorted = sort_incremented_resolution_data(orbit_data, total_blackout_time)

    max_active_times_sorted = sort_incremented_resolution_data(orbit_data, max_active_time)
    max_blackout_times_sorted = sort_incremented_resolution_data(orbit_data, max_blackout_time)

    mean_active_times_sorted = sort_incremented_resolution_data(orbit_data, mean_active_time)
    mean_blackout_times_sorted = sort_incremented_resolution_data(orbit_data, mean_blackout_time)

    mean_range_sorted = sort_incremented_resolution_data(orbit_data, mean_range)
    link_efficiency_sorted = sort_incremented_resolution_data(orbit_data, mean_link_efficiency)
    power_received_sorted = sort_incremented_resolution_data(orbit_data, mean_power_received)
    total_energy_sorted = sort_incremented_resolution_data(orbit_data, total_energy)

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
    # Define weighted objective function
    weight_active_time = 0.525
    weight_energy = 0.475
    weight_blackout_time = 0.0
    # Normalize design variables for objective function evaluation
    active_times_normalized = total_active_times_sorted / np.nanmax(total_active_times_sorted)
    total_energy_normalized = total_energy_sorted / np.nanmax(total_energy_sorted)
    blackout_times_normalized = mean_blackout_times_sorted / np.nanmax(mean_blackout_times_sorted)
    # Calculate objective function
    objective = np.array([weight_active_time * x + weight_energy * y + weight_blackout_time * z
                          for x, y, z in zip(active_times_normalized,
                                             total_energy_normalized,
                                             (1 - blackout_times_normalized))])

    # Find best orbit according to weighted objective function
    from numpy import unravel_index
    best_orbit_idx = unravel_index(np.nanargmax(objective), objective.shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]

    # Display results
    print('Generator and transmitter combined mass --> {} kg'.format(round(generator_mass + trans_mass, 2)))
    print('Best orbit --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Total active time --> {} hours, or {}%'.format(round(total_active_times_sorted[best_orbit_idx] / 3600.0, 2), round(
        100.0 * total_active_times_sorted[best_orbit_idx] / total_duration, 2)))
    print('Max active period duration --> {} hours'.format(round(max_active_times_sorted[best_orbit_idx] / 3600.0, 2)))
    print('Total blackout time --> {} hours, or {} %'.format(round(total_blackout_times_sorted[best_orbit_idx] / 3600.0, 2), round(
        100.0 * total_blackout_times_sorted[best_orbit_idx] / total_duration, 2)))
    print('Max blackout period duration --> {} hours'.format(round(max_blackout_times_sorted[best_orbit_idx] / 3600.0, 2)))
    print('Link efficiency --> {}%'.format(round(link_efficiency_sorted[best_orbit_idx] * 100.0, 5)))
    print('Power received per access period --> {} W'.format(round(power_received_sorted[best_orbit_idx], 2)))
    print('Average power received throughout mission --> {} MJ per year'.format(round(total_energy_sorted[best_orbit_idx] / 2e6, 2)))

    # Reduce perigee and apogee to altitudes instead of radii
    perigee_altitudes = [i - r_moon for i in unique_perigees]
    apogee_altitudes = [i - r_moon for i in unique_apogees]

    # make_contour_plot(perigee_altitudes, apogee_altitudes, [i / total_duration for i in total_energy_sorted], "Average Power [W]", 1)
    # make_contour_plot(perigee_altitudes, apogee_altitudes, total_energy_sorted, "Total Energy [J]", 2)
    # make_contour_plot(perigee_altitudes, apogee_altitudes, [i / 3600.0 for i in max_blackout_times_sorted], "Max Blackout Durations [h]", 3)

    # Weighted objective
    plt.figure(4)
    plt.contourf(apogee_altitudes, perigee_altitudes, objective, 500)
    plt.colorbar()
    plt.scatter(best_apogee - r_moon, best_perigee - r_moon, marker='x')
    textstr = 'Mean Power Received: {} W\nTotal Active Time Weighting: {}%\nMean Blackout Time Weighting {}%\nTotal Energy Weighting: {}%'\
        .format(round(power_received_sorted[best_orbit_idx], 2), weight_active_time * 100.0, weight_blackout_time * 100.0, weight_energy * 100.0)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(100, 4000, textstr, fontsize=7, verticalalignment='top', bbox=props)
    plt.title("Weighted Objective Function")
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.show()


main()
