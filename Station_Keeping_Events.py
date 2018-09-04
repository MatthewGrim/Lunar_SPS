""""
01/06/2018
Author: Darian van Paridon

This script evaluates the time available for performing station keeping meneouvers powered by solar radiation. Events
during which station keeping can be performed are defined as the intersection of SPS lighting events, and SPS no-access
to target events.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import *
import os


def main():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()
    r_moon = 1737.0

    # Select orbit to calculate station keeping for
    apogee_altitude = 3000.0
    perigee_altitude = 300.0

    # Get pathway to main Lunar_SPS directory
    main_directory = os.getcwd()

    # Name of study
    study_name = 'SouthPole_IncrementedRes_Inertial'

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Import access and lighting for SPS
    sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, perigee_altitude + r_moon, apogee_altitude + r_moon), start)
    sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, perigee_altitude + r_moon, apogee_altitude + r_moon), start)

    # Station keeping events
    sps_station_keeping_events = determine_battery_chargeup_events(sps_lighting, sps_access, total_duration)

    print('Total time available for station keeping maneuvers: {} hours per year'.format(round(np.sum(sps_station_keeping_events[2]) / 7200.0, 2)))

    plt.bar([i / 86400.0 for i in sps_station_keeping_events[0]], [j / 3600.0 for j in sps_station_keeping_events[2]])
    plt.xlabel('Days')
    plt.ylabel('Duration [h]')
    plt.title('Summary of Station Keeping Events')
    plt.show()


main()
