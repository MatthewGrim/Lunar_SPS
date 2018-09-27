"""
Author: Rohan Ramasamy
Date: 04/03/2018

This is a data parsing script for reeding in a csv with access times from the lunar sps study. This particular analysis
compares against the results in figure 2 of:

"A solar electric propulsion mission for lunar power beaming" - H. W. Brandhorst

"""

import numpy as np
import matplotlib.pyplot as plt
import datetime

from Lunar_SPS.general_functions_v2 import *


def compute_access_to_target(access_fname, target_illumination_fname, satellite_illumination_fname, sim_start_time):
    """
    Get the access times between a target and satellite subject to a set of constraints, concerning illumination.
    Access points when the target can be illuminated are excluded. Access points when the satellite is not
    illuminated can also be excluded.

    :param access_fname: file containing access times between satellite and target
    :param target_illumination_fname: file containing illumination times of target
    :param satellite_illumination_fname: file containing illumination times of satellite
    :param sim_start_time: start time of simulation
    :return:
    """
    access_times = get_access_times(access_fname)
    check_event_order_consistency(access_times)

    satellite_illumination_times = get_illumination_event_times_from_file(satellite_illumination_fname)
    check_event_order_consistency(satellite_illumination_times)

    target_illumination_times = get_illumination_event_times_from_file(target_illumination_fname, invert=False)
    check_event_order_consistency(target_illumination_times)

    target_eclipse_times = get_illumination_event_times_from_file("Target1-eclipse.csv", invert=False)
    # target_eclipse_times = get_illumination_event_times_from_file(target_illumination_fname, invert=True)
    check_event_order_consistency(target_eclipse_times)

    # # Filter out times when satellite is not illuminated
    #     # filter_satellite_illumination = kwargs.pop("filter_satellite_illumination", False)
    #     # if filter_satellite_illumination:
    #     #     access_times = get_overlap_between_events(access_times, satellite_illumination_times)
    #     # check_event_order_consistency(access_times)
    #     #
    #     # # Get times when target power is available
    #     # target_illumination_times = combine_events(access_times, target_illumination_times)
    #     # check_event_order_consistency(target_illumination_times)
    #     #
    #     # # Filter out times when power is not needed by target
    #     # filter_target_eclipses = kwargs.pop("filter_target_illumination", False)
    #     # if filter_target_eclipses:
    #     #     access_times = get_overlap_between_events(access_times, target_eclipse_times)
    #     # check_event_order_consistency(access_times)

    sps_sunlit_and_access = get_overlap_between_events(access_times, satellite_illumination_times)
    active_times = get_overlap_between_events(sps_sunlit_and_access, target_eclipse_times)

    access_durations = list()
    for i, active_time in enumerate(active_times["Start"]):
        access_duration = (active_times["End"][i] - active_times["Start"][i]).total_seconds()
        access_durations.append(access_duration)
    # last_eclipse_start = sim_start_time
    # eclipse_durations = list()
    # for j, target_power_start in enumerate(target_illumination_times["Start"]):
    #     target_power_end = target_illumination_times["End"][j]
    #     eclipse_time = (target_power_start - last_eclipse_start).total_seconds()
    #     eclipse_durations.append(eclipse_time)
    #     last_eclipse_start = target_power_end

    # return np.sum(np.asarray(access_durations)), np.sum(np.asarray(eclipse_durations))
    print(access_durations)
    return np.sum(np.asarray(access_durations)), 1.0


def perigee_angle_scan():
    """
    Function to replicate the results from Branhorst paper
    """

    degrees = np.asarray([15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345])
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    total_access_times = np.zeros(degrees.shape)
    total_eclipse_times = np.zeros(degrees.shape)

    for i, degree in enumerate(degrees):
        access_fname = "SPS2-Access({}).csv".format(degree)
        illumination_fname = "Target1-Lighting-Edited.csv"
        sat_illumination_fname = "SPS2-Lighting({}).csv".format(degree)

        total_access, total_eclipse = compute_access_to_target(access_fname, illumination_fname, sat_illumination_fname, start)

        print(total_access)
        total_access_times[i] = total_access / 3600
        total_eclipse_times[i] = total_eclipse / 3600

    plt.figure()
    plt.plot(degrees, total_access_times)
    plt.scatter(degrees, total_access_times)
    plt.xlabel('Relative Argument of Perigee')
    plt.ylabel('Total Access Time [hrs]')
    plt.show()


if __name__ == '__main__':
    perigee_angle_scan()

