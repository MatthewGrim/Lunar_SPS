""""
21/05/2018
Author: Darian van Paridon

This script is for analyzing a single strip through the lunar south pole data set to compare the maximum blackout durations
to the mean blackout durations as a function of the perigee altitude of an orbit with fixed apogee altitude.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import vary_orbital_elements_incrementing_resolution
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

    # Set limit for orbit size to select out of data set which is contains orbits up to 5000 km altitude.
    max_perigee = 1700.0
    max_apogee = 1700.0

    # Get orbit data set
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # Find unique perigees and apogees tested for plotting
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    r_moon = 1737.0
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    # Initialize data sets
    max_blackouts = []
    mean_blackouts = []

    # Import lunar target data
    target_lighting = parse_csv_to_array("{}/DVP_{}_Target_Lighting.csv".format(stk_data_path, study_name), start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Cycle through unique perigees within data set, with apogee fixed at max value as defined above
    for i in unique_perigees:
        # Import sps data for orbit
        sps_access = parse_csv_to_array("{}/DVP_{}_{}perigee{}apogee_access.csv".format(stk_data_path, study_name, i, max_apogee + r_moon), start)
        sps_lighting = parse_csv_to_array("{}/DVP_{}_{}perigee{}apogee_lighting.csv".format(stk_data_path, study_name, i, max_apogee + r_moon), start)
        # Determine active and blackout times
        sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
        target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        # Extract blackout statistics
        max_blackouts.append(np.max(target_blackout[2]))
        mean_blackouts.append(np.mean(target_blackout[2]))

    plt.figure(1)
    plt.subplot(211)
    plt.plot([j - r_moon for j in unique_perigees], [i / 3600.0 for i in max_blackouts])
    plt.title('Max Blackout Durations')
    plt.ylabel('Duration [hrs]')
    plt.subplot(212)
    plt.plot([j - r_moon for j in unique_perigees], [i / 3600.0 for i in mean_blackouts])
    plt.title('Mean Blackout Durations')
    plt.xlabel('Perigee Altitude')
    plt.show()


main()
