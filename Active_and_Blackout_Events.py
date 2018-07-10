""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and blackout periods for a given configuration

"""

from DVP_general_SPS_functions import *
import os


def main():
    # Import data and set the start and end times of the simulation
    # Get pathway to main Lunar_SPS directory
    main_directory = os.getcwd()

    # Name of study
    study_name = 'Equatorial_IncrementedRes'

    r_moon = 1737.0

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Solar Power Satellite 1
    high_sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, 1100.0 + r_moon, 2000.0 + r_moon), start)
    high_sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, 1100.0 + r_moon, 2000.0 + r_moon), start)

    # Solar Power Satellite 1
    low_sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, 1300.0 + r_moon, 2000.0 + r_moon), start)
    low_sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, 1300.0 + r_moon, 2000.0 + r_moon), start)

    # Lunar Target
    sunlight_target = parse_csv_to_array('{}/DVP_{}_Target_lighting.csv'.format(stk_data_path, study_name), start)
    eclipse_target = invert_events_list(sunlight_target, total_duration)

    low_sps_active = determine_SPS_active_time(low_sps_lighting, eclipse_target, low_sps_access)
    high_sps_active = determine_SPS_active_time(high_sps_lighting, eclipse_target, high_sps_access)

    low_blackout = determine_blackout_data(low_sps_active, eclipse_target, total_duration)
    high_blackout = determine_blackout_data(high_sps_active, eclipse_target, total_duration)

    plt.figure(1)
    plt.subplot(211)
    plt.bar([i / 86400.0 for i in low_blackout[0]], [j / 3600.0 for j in low_blackout[2]])
    plt.title('900 km Perigee Altitude')
    plt.ylabel('Blackout Duration [h]')
    plt.subplot(212)
    plt.bar([i / 86400.0 for i in high_blackout[0]], [j / 3600.0 for j in high_blackout[2]])
    plt.xlabel('Days')
    plt.title('1100 km Perigee Altitude')
    plt.ylabel('Blackout Duration [h]')

    plt.show()


main()
