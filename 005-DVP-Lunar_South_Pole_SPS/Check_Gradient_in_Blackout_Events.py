""""
21/05/2018
Author: Darian van Paridon

This script is for analyzing a single strip through the lunar south pole data set to compare the maximum blackout durations
to the mean blackout durations as a function of the perigee altitude of an orbit with fixed apogee altitude.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import vary_orbital_elements_incrementing_resolution, read_data_from_file, sort_incremented_resolution_data
import os


def main():
    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)
    # Name study
    study_name = 'SouthPole_IncrementedRes_Inertial'
    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    # Initialize simulation times
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Select perigee altitude along which to check gradients in blackout event durations
    perigee = 40.0

    # Get orbit data set
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(5000.0, 5000.0)

    # Find unique perigees and apogees tested for plotting
    r_moon = 1737.0
    perigee_altitudes = [orbit_data[1][0] - r_moon]
    apogee_altitudes = [orbit_data[1][1] - r_moon]
    for i in range(1, len(orbit_data)):
        if (orbit_data[i][0] - r_moon) > max(perigee_altitudes):
            perigee_altitudes.append(orbit_data[i][0] - r_moon)
        if (orbit_data[i][1] - r_moon) > max(apogee_altitudes):
            apogee_altitudes.append(orbit_data[i][1] - r_moon)

    # Import blackout event data for this study
    data_set = {}
    data_set['total_blackout_time'] = sort_incremented_resolution_data(orbit_data, read_data_from_file(stk_data_path, study_name, "TotalBlackout"))
    data_set['mean_blackout_time'] = sort_incremented_resolution_data(orbit_data, read_data_from_file(stk_data_path, study_name, "MeanBlackout"))
    data_set['max_blackout_time'] = sort_incremented_resolution_data(orbit_data, read_data_from_file(stk_data_path, study_name, "MaxBlackout"))

    # Get max blackout data and mean blackout data strips nearest selected perigee altitude
    g = min(enumerate(perigee_altitudes), key=lambda x: abs(x[1] - perigee))
    max_blackouts = data_set['max_blackout_time'][:][g[0]]
    mean_blackouts = data_set['mean_blackout_time'][:][g[0]]
    total_blackouts = data_set['total_blackout_time'][:][g[0]]

    plt.figure(1)
    plt.subplot(311)
    plt.plot(apogee_altitudes, [i / 3600.0 for i in max_blackouts])
    plt.title('Max Blackout Durations')
    plt.ylabel('Hours')
    plt.subplot(312)
    plt.plot(apogee_altitudes, [i / 3600.0 for i in mean_blackouts])
    plt.title('Mean Blackout Durations')
    plt.ylabel('Hours')
    plt.subplot(313)
    plt.plot(apogee_altitudes, [100.0 * i / total_duration for i in total_blackouts])
    plt.title('Total Blackout Time')
    plt.ylabel('Percent')
    plt.xlabel('Perigee Altitude')
    plt.show()


main()
