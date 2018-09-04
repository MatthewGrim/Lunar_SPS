""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and blackout periods for a given configuration

"""

from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import vary_orbital_elements_incrementing_resolution
import os


def main():
    # Import data and set the start and end times of the simulation
    # Get pathway to main Lunar_SPS directory
    main_directory = os.getcwd()

    # Name of study
    study_name = 'SouthPole_IncrementedRes_ManytoOne'

    r_moon = 1737.0

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Lunar Target
    sunlight_target = parse_csv_to_array('{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name), start)
    eclipse_target = invert_events_list(sunlight_target, total_duration)

    orbit_data = vary_orbital_elements_incrementing_resolution(5000.0, 5000.0)

    # Solar Power Satellite 1
    sps_two_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_0.0meananom_access.csv'.format(stk_data_path, study_name, 750.0 + r_moon, 2800.0 + r_moon), start)
    sps_two_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_0.0meananom_lighting.csv'.format(stk_data_path, study_name, 750.0 + r_moon, 2800.0 + r_moon), start)
    sps_two_active = determine_SPS_active_time(sps_two_lighting, eclipse_target, sps_two_access)
    sps_two_blackout = determine_blackout_data(sps_two_active, eclipse_target, total_duration)

    # Solar Power Satellite 1
    sps_one_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_180.0meananom_access.csv'.format(stk_data_path, study_name, 750.0 + r_moon, 2800.0 + r_moon), start)
    sps_one_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_180.0meananom_lighting.csv'.format(stk_data_path, study_name, 750.0 + r_moon, 2800.0 + r_moon), start)
    sps_one_active = determine_SPS_active_time(sps_one_lighting, eclipse_target, sps_one_access)
    sps_one_blackout = determine_blackout_data(sps_one_active, eclipse_target, total_duration)
    sps_one_range = import_range_data_statistics('DVP_{}_{}perigee{}apogee_180.0meananom_range'.format(study_name, 750.0 + r_moon, 2800.0 + r_moon), stk_data_path)

    mean_range_one = np.sum([(i * j) / np.sum(sps_one_access[2]) for i, j in zip(sps_one_range[2], sps_one_access[2])])
    print(mean_range_one)

    combined_active = combine_events(sps_one_active, sps_two_active)
    combined_blackout = determine_blackout_data(combined_active, eclipse_target, total_duration)

    plt.figure(1)
    plt.subplot(311)
    plt.bar([i / 86400.0 for i in sps_one_active[0]], [j / 3600.0 for j in sps_one_active[2]])
    plt.title('Active Events for Each SPS, and Combined SPS')
    plt.ylabel('Hours')
    plt.subplot(312)
    plt.bar([i / 86400.0 for i in sps_two_active[0]], [j / 3600.0 for j in sps_two_active[2]])
    plt.ylabel('Hours')
    plt.subplot(313)
    plt.bar([i / 86400.0 for i in combined_active[0]], [j / 3600.0 for j in combined_active[2]])
    plt.ylabel('Hours')
    plt.xlabel('Days')

    plt.figure(2)
    plt.subplot(311)
    plt.bar([i / 86400.0 for i in sps_one_blackout[0]], [j / 3600.0 for j in sps_one_blackout[2]])
    plt.title('Blackout Events for Each SPS, and Combined SPS')
    plt.ylabel('Hours')
    plt.subplot(312)
    plt.bar([i / 86400.0 for i in sps_two_blackout[0]], [j / 3600.0 for j in sps_two_blackout[2]])
    plt.ylabel('Hours')
    plt.subplot(313)
    plt.bar([i / 86400.0 for i in combined_blackout[0]], [j / 3600.0 for j in combined_blackout[2]])
    plt.ylabel('Hours')
    plt.xlabel('Days')
    plt.show()


main()
