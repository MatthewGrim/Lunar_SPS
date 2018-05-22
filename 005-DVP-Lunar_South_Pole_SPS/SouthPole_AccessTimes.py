""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the access time for a lunar SPS in polar orbit,
targeting the south pole. Satellite flies at 90 degree inclination, with 90
degree argument of perigee

"""

from DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def scan_apogee_altitude(start, total_duration, scan_array):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various arguments of perigee
    total_active_time = np.zeros(len(scan_array))
    total_blackout_time = np.zeros(len(scan_array))
    max_blackout_time = np.zeros(len(scan_array))
    raw_target_lighting = "Target_Lighting.csv"
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    eclipse_times = invert_events_list(target_lighting, total_duration)
    for i in range(len(scan_array)):
        raw_sunlight = "SPS(300,{}k)_Lighting.csv".format(scan_array[i])
        raw_access = "SPS(300,{}k)_Access.csv".format(scan_array[i])

        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)

        sps_access_sunlit_times = get_event_overlaps(sunlight_times, access_times)
        active_times = get_event_overlaps(sps_access_sunlit_times, eclipse_times)

        print('\n')
        print('Altitude {}000 km'.format(scan_array[i]))
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_times, access_times)

        max_blackout_time[i] = determine_blackout_data(active_times, eclipse_times, total_duration)

    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(1)
    plt.subplot(211)
    plt.plot(scan_array * 1000.0, total_active_time / 3600.0)
    plt.xlabel('Altitude [km]')
    plt.ylabel('Total Active Time [hrs]')
    plt.title('Circular Orbits')
    plt.subplot(212)
    plt.plot(scan_array * 1000.0, max_blackout_time / 3600.0)
    plt.xlabel('Altitude [km]')
    plt.ylabel('Maximum Blackout Duration [hrs]')
    plt.show()

    return scan_array, total_active_time, total_blackout_time


def main():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    altitudes = np.asarray([2, 3, 4, 5, 6, 8, 10])
    apogees = np.array([5, 6, 7, 8, 10, 15, 20])
    scan_apogee_altitude(start, total_duration, apogees)
    # sps_access_raw = 'SPS(3k)_Access_89.45-222.69.csv'
    # sps_lighting_raw = 'SPS(3k)_Lighting_89.45-222.69.csv'
    # target_lighting_raw = 'Target_Lighting_89.45-222.69.csv'
    # sps_access = parse_csv_to_array(sps_access_raw, start)
    # sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    # target_lighting = parse_csv_to_array(target_lighting_raw, start)
    # target_eclipse = invert_events_list(target_lighting, total_duration)

    # print("\n")
    # print("For a 3000 km circular orbit:")
    # sps_active_circular_orbit = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
    # print('Maximum eclipse duration without SPS is {} hours'.format(round(max(target_eclipse[2]) / 3600.0, 2)))
    # sps_sunny_access_times = get_event_overlaps(sps_lighting, sps_access)
    # active_times_circle_orbit = get_event_overlaps(sps_sunny_access_times, target_eclipse)
    # determine_blackout_data(active_times_circle_orbit, target_eclipse, total_duration)


main()
