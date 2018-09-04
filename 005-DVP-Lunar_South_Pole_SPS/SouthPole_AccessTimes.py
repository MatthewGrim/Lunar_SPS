""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the access time for a lunar SPS in polar orbit,
targeting the south pole. Satellite flies at 90 degree inclination, with 90
degree argument of perigee

"""

from Lunar_SPS.DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def scan_apogee(start, total_duration, scan_array):
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
        print('Maximum black-out duration w/o SPS: {} hrs'.format(round(max(eclipse_times[2]) / 3600.0, 2)))
        print('Total time spent in black-out w/o SPS: {}% per year'.format(round(100.0 * np.sum(eclipse_times[2]) / (2.0 * 24.0 * 365.0 * 3600.0), 2)))
        activity = determine_SPS_active_time(sunlight_times, eclipse_times, access_times)
        total_active_time[i] = np.sum(activity[2])
        blacktivity = determine_blackout_data(active_times, eclipse_times, total_duration)
        max_blackout_time[i] = max(blacktivity[2])
    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(1)
    plt.subplot(211)
    plt.plot(scan_array * 1000.0, total_active_time / 3600.0)
    plt.ylabel('Total Active Time [hrs]')
    plt.title('Fixed Perigee 300 km')
    plt.subplot(212)
    plt.plot(scan_array * 1000.0, max_blackout_time / 3600.0)
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Maximum Blackout Duration [hrs]')
    plt.show()

    return scan_array, total_active_time, total_blackout_time


def scan_altitude(start, total_duration, scan_array):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various arguments of perigee
    max_active_time = np.zeros(len(scan_array))
    total_active_time = np.zeros(len(scan_array))
    max_blackout_time = np.zeros(len(scan_array))
    tot_blackout_time = np.zeros(len(scan_array))
    tot_original_blackout = np.zeros(len(scan_array))
    max_original_blackout = np.zeros(len(scan_array))

    raw_target_lighting = "Target_Lighting.csv"
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    eclipse_times = invert_events_list(target_lighting, total_duration)
    for i in range(len(scan_array)):
        raw_sunlight = "SPS({}k)_Lighting.csv".format(scan_array[i])
        raw_access = "SPS({}k)_Access.csv".format(scan_array[i])

        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)

        sps_access_sunlit_times = get_event_overlaps(sunlight_times, access_times)
        active_times = get_event_overlaps(sps_access_sunlit_times, eclipse_times)

        print('\n')
        print('Altitude {} km'.format(scan_array[i]*1000))
        print('Maximum black-out duration w/o SPS: {} hrs'.format(round(max(eclipse_times[2]) / 3600.0, 2)))
        print('Total time spent in black-out w/o SPS: {}% per year'.format(round(100.0 * np.sum(eclipse_times[2]) / (2.0 * 24.0 * 365.0 * 3600.0), 2)))
        activity = determine_SPS_active_time(sunlight_times, eclipse_times, access_times)
        max_active_time[i] = max(activity[2])
        total_active_time[i] = np.sum(activity[2])
        blacktivity = determine_blackout_data(active_times, eclipse_times, total_duration)
        max_blackout_time[i] = max(blacktivity[2])
        tot_blackout_time[i] = np.sum(blacktivity[2])
        tot_original_blackout[i] = np.sum(eclipse_times[2])
        max_original_blackout[i] = max(eclipse_times[2])

    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(1)
    plt.subplot(211)
    plt.plot(scan_array * 1000.0, 100.0 * total_active_time / (2.0 * 365.0 * 24.0 * 3600.0))
    plt.ylabel('Total Active Time [%/yr]')
    plt.title('Circular Orbits')
    plt.subplot(212)
    plt.plot(scan_array * 1000.0, 100.0 * tot_blackout_time / (2.0 * 365.0 * 24.0 * 3600.0), label='With SPS')
    plt.plot(scan_array * 1000.0, 100.0 * tot_original_blackout / (2.0 * 365.0 * 24.0 * 3600.0), label='No SPS', linestyle=':')
    plt.legend()
    plt.ylim([25, 55])
    plt.xlabel('Altitude [km]')
    plt.ylabel('Total Blackout Duration [%/yr]')
    plt.show()

    plt.figure(2)
    plt.subplot(211)
    plt.plot(scan_array * 1000.0, max_active_time / 3600.0)
    plt.ylabel('Max Active Duration [hrs]')
    plt.title('Circular Orbits')
    plt.subplot(212)
    plt.plot(scan_array * 1000.0, max_blackout_time / 3600.0)
    plt.legend()
    plt.xlabel('Altitude [km]')
    plt.ylabel('Max Blackout Duration [hrs]')
    plt.show()

    return scan_array, total_active_time, tot_blackout_time


def main():
    start = convert_string_to_datetime(['2018', '05', '18', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '18', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    altitudes = np.asarray([1, 1.5, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.5, 5, 6, 7, 8, 9, 10])
    # scan_altitude(start, total_duration, altitudes)

    apogees = np.array([5, 6, 7, 8, 10, 15, 20])
    scan_apogee(start, total_duration, apogees)


main()
