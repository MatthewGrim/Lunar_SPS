""""
01/06/2018
Author: Darian van Paridon

This script imports and process the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the Brandhorst lunar SPS configuration.
"""

from DVP_general_SPS_functions import *
import os as os
import math


def vary_orbital_elements(resolution):

    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK

    # This is a copy-paste from Brandhorst_Programmatic_Script. Ensure parameters are the same as
    # they were when the reports were generated.

    radius_moon = 1737.0
    # Set resolution of data points in km
    max_perigee = 10000.0
    min_perigee = 5000.0
    min_apogee = 45000.0
    max_apogee = 55000.0
    orbit_data = [0.0, 0.0]

    for j in range(1, int((max_perigee - min_perigee) / resolution) + 1):
        perigee = min_perigee + (j * resolution)
        apogee = min_apogee
        while apogee <= max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            apogee += resolution

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])

    return semi_maj_axis, eccentricity, orbit_data


def main():
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Difference in subsequent apogee/perigee radii
    resolution = 500.0

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name of study
    study_name = 'Equatorial_{}kmRes_5min_Check'.format(resolution)

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Get orbital data
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements(resolution)

    # # Import report data generated via STK Connect Commands
    total_active_time = []
    max_active_duration = []
    mean_link_efficiency = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)

    # PARAMETRIC SCAN
    ####################################################################################################################
    # Cycle through available orbit configurations and calculate active/blackout durations
    # Comment out if reading processed data in from a txt file
    for i in range(1, len(orbit_data)):
        print('Progress: {}%'.format(round(100.0 * (i-1) / (len(orbit_data) - 2), 2)))
        print("Perigee radius: {} km, Apogee radius: {} km".format(orbit_data[i][0], orbit_data[i][1]))

        # CALCULATING ACCESS TIME
        ################################################################################################################
        # Import SPS illumination and access events
        sps_lighting_raw = '{}/DVP_{}_{}perigee_{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
        sps_access_raw = '{}/DVP_{}_{}perigee_{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1])
        sps_access = parse_csv_to_array(sps_access_raw, start)

        # Determine the total and maximum SPS active durations
        sps_active = determine_SPS_active_time(sps_lighting, sps_access, target_lighting)
        total_active_time.append(np.sum(sps_active[2]))
        max_active_duration.append(max(sps_active[2]))
        ################################################################################################################

        # CALCULATING LINK EFFICIENCY
        ################################################################################################################
        # Transmitter parameters
        wavelength = 1070.0e-9
        trans_radius = 0.005
        # Receiver size for AMALIA rover receiver
        rec_radius = 0.5
        # Ratio of transmitter to receiver radius
        trans_rec_radius_ratio = rec_radius / trans_radius

        # Import SPS range data from statistical report
        # First entry of sps_range is minimum range for each access period
        # Second entry is maximum range for each access period
        # Third entry is mean range for each access period
        link_eff_temp = []
        sps_range = import_range_data_statistics('DVP_{}_{}perigee_{}apogee_range'.format(study_name, orbit_data[i][0], orbit_data[i][1]), stk_data_path)
        for j in sps_range[2]:
            surf_beam_radius = trans_radius * np.sqrt(
                1 + (wavelength * (j * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
            if surf_beam_radius <= rec_radius:
                link_eff_temp.append(1.0)
            else:
                link_eff_temp.append((rec_radius / surf_beam_radius) ** 2)
        # Take mean of every mean link efficiency
        mean_link_efficiency.append(np.mean(link_eff_temp))
    ####################################################################################################################

    # WRITING DATA TO FILES
    ####################################################################################################################
    # Writing processed data to files to save time in future analysis. Comment this out if the processed data files
    # already exist and are being read out (see next section).

    # MEAN LINK EFFICIENCY
    try:
        os.remove('{}/{}_TransRecRatio_MeanLinkEfficiency_{}.txt'.format(stk_data_path, trans_rec_radius_ratio, study_name))
    except OSError:
        pass
    mean_link_efficiency = np.asarray(mean_link_efficiency)
    np.savetxt('{}/{}_RecTransRatio_MeanLinkEfficiency_{}.txt'.format(stk_data_path, trans_rec_radius_ratio, study_name), mean_link_efficiency)

    # TOTAL ACTIVE TIME
    try:
        os.remove('{}/TotalActive_{}.txt'.format(stk_data_path, study_name))
    except OSError:
        pass
    total_active_time = np.asarray(total_active_time)
    np.savetxt('{}/TotalActive_{}.txt'.format(stk_data_path, study_name), total_active_time)

    # MAX ACTIVE TIME
    try:
        os.remove('{}/MaxActive_{}.txt'.format(stk_data_path, study_name))
    except OSError:
        pass
    max_active_time = np.asarray(max_active_duration)
    np.savetxt('{}/MaxActive_{}.txt'.format(stk_data_path, study_name), max_active_time)
    ####################################################################################################################

    # READ IN DATA FILES IF THEY ARE ALREADY WRITTEN
    ####################################################################################################################
    # Read in saved active time data, if the data is already processed and saved into a txt file
    # Comment this section out if the data is being processed, and the results are not yet written to a txt file

    # # MEAN LINK EFFICIENCY
    # mean_link_efficiency = []
    # fh = open('{}/{}_ RecTransRatio_MeanLinkEfficiency_{}.txt'.format(stk_data_path, round(trans_rec_radius_ratio, 1), study_name))
    # for i, line in enumerate(fh):
    #     temp = line.split("\n")
    #     mean_link_efficiency.append(float(temp[0]))
    #
    # # TOTAL ACTIVE TIME
    # total_active_time = []
    # fh = open('{}/TotalActive_{}.txt'.format(stk_data_path, study_name))
    # for i, line in enumerate(fh):
    #     temp = line.split("\n")
    #     total_active_time.append(float(temp[0]))
    ####################################################################################################################

    # REORGANIZE DATA FOR PLOTTING
    ####################################################################################################################
    # Get list of unique perigees and apogees tested
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    # Initialize arrays for plotting
    link_efficiency = np.zeros((len(unique_perigees), len(unique_apogees)))
    active_times = np.zeros((len(unique_perigees), len(unique_apogees)))
    objective = np.zeros((len(unique_perigees), len(unique_apogees)))

    # Sort data into a arrays for 2D plot
    max_apogee = max(unique_apogees)
    num_apogees = list()
    for k in range(len(unique_perigees)):
        # Determine how many apogee steps required to reach the max value.
        num_apogees.append(int(math.floor((max_apogee - unique_perigees[k]) / resolution)))
        # For perigee value (index k), collect data points for corresponding apogees (index j)
        # where j starts at the first unique apogee equal to the current perigee
        j = 0
        while unique_apogees[j] < unique_perigees[k]:
            j += 1
        start = j
        for j in range(0, num_apogees[k] + 1):
            if k == 0:
                idx = j + start
            else:
                idx = int(np.sum(num_apogees) - num_apogees[k] + j + start)
            active_times[k, j + start] = total_active_time[idx] / 3600.0
            link_efficiency[k, j + start] = mean_link_efficiency[idx]
            # Compute objective function, which is equally weighted on total access time and mean link efficiency
            objective[k, j + start] = 0.5 * ((mean_link_efficiency[idx] / np.max(mean_link_efficiency)) + (total_active_time[idx] / np.max(total_active_time)))

    # # SORTING ALGORITHM VERSION II - FOR CHECKING OSCILLATION IN THE EQUATORIAL CASE
    # max_apogee = max(unique_apogees)
    # num_apogees = list()
    # resolution = 500.0
    # for k in range(len(unique_perigees)):
    #     for j in range(len(unique_apogees)):
    #         idx = int((np.sum(num_apogees) * k) + j)
    #         active_times[k, j] = total_active_time[idx] / 3600.0
    #         link_efficiency[k, j] = mean_link_efficiency[idx]
    #         objective[k, j] = 0.5 * ((mean_link_efficiency[idx] / np.max(mean_link_efficiency)) + (total_active_time[idx] / np.max(total_active_time)))

    # Find orbit which results in highest performance according to objective function
    from numpy import unravel_index
    best_orbit_idx = unravel_index(objective.argmax(), objective.shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    print('Best orbit --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee, 2), round(best_apogee, 2)))
    print('Total active time --> {} hours, or {}%'.format(active_times[best_orbit_idx], round(active_times[best_orbit_idx] * 3600.0 / total_duration, 2)))
    print('Link efficiency --> {}%'.format(round(link_efficiency[best_orbit_idx] * 100.0)))
    ####################################################################################################################

    # PLOTTING THE ACTIVE TIMES, LINK EFFICIENCIES, AND OBJECTIVES
    ####################################################################################################################
    # Generate figure showing access time as a function of perigee/apogee
    plt.figure(1)
    plt.contourf(unique_apogees, unique_perigees, active_times, 1000)
    plt.colorbar()
    plt.title('Total Active Time [hrs]')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()

    # Generate figure showing mean link efficiency as a function of perigee/apogee
    plt.figure(2)
    plt.contourf(unique_apogees, unique_perigees, link_efficiency * 100.0, 1000)
    plt.colorbar()
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()

    # Generate figure showing objetive function (combining total access and link efficiency)
    # Mark orbit which results in highest performance
    plt.figure(3)
    plt.contourf(unique_apogees, unique_perigees, objective, 1000)
    plt.colorbar()
    plt.plot(best_apogee, best_perigee, marker='x', markersize=5, color='green')
    plt.title('Sum of Mean Link Efficiency and Total Active Time')
    plt.xlabel('Apogee (Dist to Focus) [km]')
    plt.ylabel('Perigee (Dist to Focus) [km]')
    plt.show()


main()
