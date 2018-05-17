""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the access time for a lunar SPS in polar orbit,
targeting the south pole. Satellite flies at 90 degree inclination, with 90
degree argument of perigee

"""

from DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def get_times_without_coverage(access_times, duration):
    length = len(access_times[0])
    no_coverage = np.zeros((3, length))
    # Check to if first access period begins at start of simulation
    # If not, no_coverage begins at start of simulation
    if access_times[0][0] != 0:
        no_coverage[0][0] = 0
        no_coverage[1][0] = access_times[0][0]
        no_coverage[2][0] = no_coverage[1][0] - no_coverage[0][0]
        for j in range(1, length):
            no_coverage[1][j] = access_times[0][j]
            no_coverage[0][j] = access_times[1][j-1]
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]
        j = length - 1
        if access_times[1][j] == duration:
            no_coverage[0][j] = 'nan'
            no_coverage[1][j] = 'nan'
            no_coverage[2][j] = 'nan'
        elif access_times[1][j] != duration:
            no_coverage[0][j] = access_times[1][j]
            no_coverage[1][j] = duration
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]

    # If access period begins at start of simulation
    elif access_times[0][0] == 0:
        for j in range(length - 2):
            no_coverage[0][j] = access_times[1][j]
            no_coverage[1][j] = access_times[0][j+1]
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]
        j = length- 1
        if access_times[1][j] == duration:
            no_coverage[0][j] = 'nan'
            no_coverage[1][j] = 'nan'
            no_coverage[2][j] = 'nan'
        elif access_times[1][j] != duration:
            no_coverage[0][j] = access_times[1][j]
            no_coverage[1][j] = duration
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]

    return no_coverage


def determine_blackout_data(active_times, eclipse_target,  duration):

    # This function replicates the results of Figure 2 in Brandhorst's paper
    # regarding black-out periods for the lunar target (ie. times when there is
    # no SPS coverage, and the target is in eclipse.

    # Calculates the times when the SPS cannot access the target
    SPS_no_coverage = get_times_without_coverage(active_times, duration)
    check_sum(active_times, SPS_no_coverage, duration)
    # Calculates overlap of no coverage with target eclipses
    dark_periods = get_event_overlaps(SPS_no_coverage, eclipse_target)
    # Calculate number of eclipses exceeding 84 hours, and the maximum eclipse duration
    temp_array = np.array(dark_periods[2])
    long_eclipse_flag = (temp_array / 3600.0) > 5.0
    num_long_eclipse = np.sum(long_eclipse_flag)
    max_eclipse_duration = round(max([i/3600.0 for i in dark_periods[2]]), 2)
    print('Maximum black-out duration with SPS: {} hrs'.format(max_eclipse_duration))
    print("Number of blackouts exceeding 5 hrs: {}".format(num_long_eclipse))
    print("\n")
    print(temp_array)
    return sum(temp_array)

def scan_apogee_altitude(start, total_duration, eclipse_target):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various arguments of perigee
    apogees = np.asarray([5, 6, 7, 8, 9, 10, 15, 20])
    total_active_time = np.zeros(len(apogees))
    total_blackout_time = np.zeros(len(apogees))
    for i in range(len(apogees)):
        raw_sunlight = "SPS(300,{}k)_Lighting.csv".format(apogees[i])
        raw_eclipse = "Target_Eclipse.csv".format(apogees[i])
        raw_access = "SPS(300,{}k)_Access.csv".format(apogees[i])

        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)
        eclipse_times = parse_csv_to_array(raw_eclipse, start)

        sps_access_sunlit_times = get_event_overlaps(sunlight_times, access_times)
        active_times = get_event_overlaps(sps_access_sunlit_times, eclipse_target)

        print('\n')
        print('Apogee Altitude {}000 km'.format(apogees[i]))
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_target, access_times)
        print('Maximum eclipse duration without SPS is {} hours'.format(round(max(eclipse_target[2]) / 3600.0, 2)))
        total_blackout_time[i] = determine_blackout_data(active_times, eclipse_times, total_duration)

    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(1)
    plt.subplot(211)
    plt.plot(apogees, total_active_time / 3600.0)
    plt.xlabel('Apogee Altitude (Above Target) [km]')
    plt.ylabel('Total Access Time of SPS [hrs]')
    plt.subplot(212)
    plt.plot(apogees, total_blackout_time / 3600.0)
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
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    target_eclipse = parse_csv_to_array(target_eclipse_raw, start)

    scan_apogee_altitude(start, total_duration, target_eclipse)

    print("\n")
    print("For a 3000 km circular orbit:")
    sps_active_circular_orbit = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
    print('Maximum eclipse duration without SPS is {} hours'.format(round(max(target_eclipse[2]) / 3600.0, 2)))
    determine_blackout_data(sps_access, target_eclipse, total_duration)



main()
