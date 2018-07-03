""""
Author: Darian van Paridon
Date: 04/05/2018

This file processes csv data reports from STK which contain the access and lighting times for a SPS. The purpose is to
replicate Figure 3 from Brandhorst's paper " A solar electric propulsion mission for lunar power beaming " which shows
the variation in total active time as a function of relevant argument of perigee between the two satellites in the
proposed SPS constellation.

"""

from general_functions import *
from DVP_general_SPS_functions import *


def scan_perigee_angles(start, total_duration, eclipse_target):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various arguments of perigee
    relative_perigees = np.asarray([15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345])
    total_active_time = np.zeros(len(relative_perigees))
    for i in range(len(relative_perigees)):
        raw_sunlight = "SPS2-Lighting({}).csv".format(relative_perigees[i])
        raw_eclipse = "SPS2-Eclipse({}).csv".format(relative_perigees[i])
        raw_access = "SPS2-Access({}).csv".format(relative_perigees[i])

        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)
        eclipse_times = parse_csv_to_array(raw_eclipse, start)
        check_sum(sunlight_times, eclipse_times, total_duration)

        print('\n')
        print('Relative perigee {}'.format(relative_perigees[i]))
        sps_active = determine_SPS_active_time(sunlight_times, eclipse_target, access_times)
        total_active_time[i] = np.sum(sps_active[2])

    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(3)
    plt.plot(relative_perigees, total_active_time / 3600.0)
    plt.xlabel('Relative Argument of Perigee [deg]')
    plt.ylabel('Total Access Time of SPS2 [hrs]')
    plt.title("Brandhorst Figure 3 Remake")
    plt.show()

    return relative_perigees, total_active_time


def remake_brandhorst_fig3():
    # Import data and set the start and end times of the simulation
    raw_sunlight_SPS2 = "SPS2-Lighting(135).csv"
    raw_eclipse_SPS2 = "SPS2-Eclipse(135).csv"
    raw_access_SPS2 = "SPS2-Access(135).csv"
    raw_sunlight_SPS1 = "SPS1-Lighting(0)-Edited.csv"
    raw_eclipse_SPS1 = "SPS1-Eclipse(0)-Edited.csv"
    raw_access_SPS1 = "SPS1-Access(0).csv"
    raw_sunlight_target = 'Target1-Lighting-Edited.csv'
    raw_eclipse_target = 'Target1-Eclipse-Edited.csv'
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Solar Power Satellite 1
    LOS_access1 = parse_csv_to_array(raw_access_SPS1, start)
    sunlight_SPS1 = parse_csv_to_array(raw_sunlight_SPS1, start)
    eclipse_SPS1 = parse_csv_to_array(raw_eclipse_SPS1, start)

    # Solar Power Satellite 2
    LOS_access2 = parse_csv_to_array(raw_access_SPS2, start)
    sunlight_SPS2 = parse_csv_to_array(raw_sunlight_SPS2, start)
    eclipse_SPS2 = parse_csv_to_array(raw_eclipse_SPS2, start)

    # Lunar Target
    sunlight_target = parse_csv_to_array(raw_sunlight_target, start)
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)


    # Calculates the total active time for SPS, based on target access
    # and eclipses, as well as satellite illumination times
    print("\n")
    print("ACCESS AVAILABILITY for SPS1")
    determine_SPS_active_time(sunlight_SPS1, eclipse_target, LOS_access1)
    determine_blackout_data(LOS_access1, eclipse_target, total_duration)

    print("\n")
    print("ACCESS AVAILABILITY for SPS2")
    determine_SPS_active_time(sunlight_SPS2, eclipse_target, LOS_access2)
    determine_blackout_data(LOS_access2, eclipse_target, total_duration)

    # Calculates the total access time for various relative arguments of perigee
    perigees, active_times = scan_perigee_angles(start, total_duration, eclipse_target)

    return perigees, active_times


remake_brandhorst_fig3()
