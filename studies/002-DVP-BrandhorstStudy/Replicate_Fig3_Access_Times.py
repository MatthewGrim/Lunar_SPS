""""
Author: Darian van Paridon
Date: 04/05/2018

This file processes csv data reports from STK which contain the access and lighting times for a SPS. The purpose is to
replicate Figure 3 from Brandhorst's paper " A solar electric propulsion mission for lunar power beaming " which shows
the variation in total active time as a function of relevant argument of perigee between the two satellites in the
proposed SPS constellation.

"""

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *
from Lunar_SPS.pysrc.post_process_functions.general_functions import *


def scan_perigee_angles(start, total_duration, target_eclipse):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various arguments of perigee
    relative_perigees = np.asarray([15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345])
    total_active_time = np.zeros(len(relative_perigees))
    for i in range(len(relative_perigees)):
        raw_sunlight = "SPS2-Lighting({}).csv".format(relative_perigees[i])
        raw_eclipse = "SPS2-Eclipse({}).csv".format(relative_perigees[i])
        raw_access = "SPS2-Access({}).csv".format(relative_perigees[i])

        sps_access = parse_csv_to_array(raw_access, start)
        sps_lighting = parse_csv_to_array(raw_sunlight, start)
        sps_eclipse = parse_csv_to_array(raw_eclipse, start)
        check_sum(sps_lighting, sps_eclipse, total_duration)

        print('\n')
        print('Relative perigee {}'.format(relative_perigees[i]))
        sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
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
    raw_sps2_lighting = "SPS2-Lighting(135).csv"
    raw_access_SPS2 = "SPS2-Access(135).csv"
    raw_sps1_lighting = "SPS1-Lighting(0)-Edited.csv"
    raw_access_SPS1 = "SPS1-Access(0).csv"
    raw_target_eclipse = 'Target1-eclipse.csv'

    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Solar Power Satellite 1
    sps1_access = parse_csv_to_array(raw_access_SPS1, start)
    sps1_lighting = parse_csv_to_array(raw_sps1_lighting, start)

    # Solar Power Satellite 2
    sps2_access = parse_csv_to_array(raw_access_SPS2, start)
    sps2_lighting = parse_csv_to_array(raw_sps2_lighting, start)

    # Lunar Target
    target_eclipse = parse_csv_to_array(raw_target_eclipse, start)

    # Calculates the total active time for SPS, based on target access and eclipses, as well as satellite
    # illumination times
    print("\n")
    print("ACCESS AVAILABILITY for SPS1")
    determine_SPS_active_time(sps1_lighting, target_eclipse, sps1_access)
    determine_blackout_data(sps1_access, target_eclipse, total_duration)

    print("\n")
    print("ACCESS AVAILABILITY for SPS2")
    determine_SPS_active_time(sps2_lighting, target_eclipse, sps2_access)
    determine_blackout_data(sps2_access, target_eclipse, total_duration)

    # Calculates the total access time for various relative arguments of perigee
    perigees, active_times = scan_perigee_angles(start, total_duration, target_eclipse)

    return perigees, active_times


remake_brandhorst_fig3()
