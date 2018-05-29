""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the access time for a lunar SPS in polar orbit,
targeting the south pole. Satellite flies at 90 degree inclination, with 90
degree argument of perigee

"""

from DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def main():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    sps_access_raw = 'SPS(3k)_Access_89.45-222.69.csv'
    sps_lighting_raw = 'SPS(3k)_Lighting_89.45-222.69.csv'
    target_lighting_raw = 'Target_Lighting_Terrain_89.45-222.69.csv'
    sps_access = parse_csv_to_array(sps_access_raw, start)
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    print('Maximum blackout duration without SPS: {} hrs'.format(round(max(target_eclipse[2]) / 3600.0, 2)))
    print('Percentage of time spent in blackout w/o SPS: {}%'.format(round(100.0 * np.sum(target_eclipse[2]) / (2.0 * 3600.0 * 24.0 * 365.0), 2)))

    sps_active_time = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)

    determine_blackout_data(sps_active_time, target_eclipse, total_duration)

    plt.bar([i / 86400.0 for i in target_lighting[0]], [i / 3600.0 for i in target_lighting[2]], label='Sunlit')
    plt.bar([i / 86400.0 for i in target_eclipse[0]], [i / 3600.0 for i in target_eclipse[2]], label='Eclipse')
    plt.legend()
    plt.show()


main()
