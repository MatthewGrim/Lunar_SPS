""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the access time for a lunar SPS in polar orbit,
targeting the south pole. Satellite flies at 90 degree inclination, with 90
degree argument of perigee

"""

from DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def scan_apogee_altitude(start, total_duration):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various arguments of perigee
    apogees = np.asarray([5, 6, 7, 8, 9, 10, 15, 20])
    total_active_time = np.zeros(len(apogees))
    total_blackout_time = np.zeros(len(apogees))
    raw_eclipse = "Target_Eclipse.csv"
    for i in range(len(apogees)):
        raw_sunlight = "SPS(300,{}k)_Lighting.csv".format(apogees[i])
        raw_access = "SPS(300,{}k)_Access.csv".format(apogees[i])

        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)
        eclipse_times = parse_csv_to_array(raw_eclipse, start)

        sps_access_sunlit_times = get_event_overlaps(sunlight_times, access_times)
        active_times = get_event_overlaps(sps_access_sunlit_times, eclipse_times)

        # q = [Q / 86400.0 for Q in active_times[0]]
        # w = [W / 3600.0 for W in active_times[2]]
        # plt.subplot(411)
        # plt.bar(eclipse_times[0] / 86400.0, eclipse_times[2] / 3600.0)
        # plt.subplot(412)
        # plt.bar(sunlight_times[0] / 86400.0, sunlight_times[2] / 3600.0)
        # plt.subplot(413)
        # plt.bar(access_times[0] / 86400.0, access_times[2] / 3600.0)
        # plt.subplot(414)
        # plt.bar(q, w)
        # plt.show()

        print('\n')
        print('Apogee Altitude {}000 km'.format(apogees[i]))
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_times, access_times)

        blackout_data = determine_blackout_data(active_times, eclipse_times, total_duration)
        total_blackout_time[i] = sum(blackout_data)

    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(1)
    plt.subplot(211)
    plt.plot(apogees * 1000.0, total_active_time / 3600.0)
    plt.xlabel('Apogee Altitude (Above Target) [km]')
    plt.ylabel('Total Access Time of SPS [hrs]')
    plt.subplot(212)
    plt.plot(apogees * 1000.0, total_blackout_time / 3600.0)
    plt.xlabel('Apogee Altitude (Above Target) [km]')
    plt.ylabel('Total Blackout Time at Target [hrs]')
    plt.show()

    return apogees, total_active_time, total_blackout_time


def main():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    sps_access_raw = 'SPS(3k)_Access.csv'
    sps_lighting_raw = 'SPS(3k)_Lighting.csv'
    target_eclipse_raw = 'Target_Eclipse.csv'
    sps_access = parse_csv_to_array(sps_access_raw, start)
    check_event_order_consistency(sps_access)
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    target_eclipse = parse_csv_to_array(target_eclipse_raw, start)

    scan_apogee_altitude(start, total_duration)

    print("\n")
    print("For a 3000 km circular orbit:")
    sps_active_circular_orbit = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
    print('Maximum eclipse duration without SPS is {} hours'.format(round(max(target_eclipse[2]) / 3600.0, 2)))
    sps_sunny_access_times = get_event_overlaps(sps_lighting, sps_access)
    active_times_cirorbit = get_event_overlaps(sps_sunny_access_times, target_eclipse)
    determine_blackout_data(active_times_cirorbit, target_eclipse, total_duration)


main()
