"""
Author: Rohan Ramasamy
Date: 13/03/2018

This particular analysis compares against the results in figure 2 of:

"A solar electric propulsion mission for lunar power beaming" - H. W. Brandhorst
"""

import matplotlib.pyplot as plt
import numpy as np
from Lunar_SPS.pysrc.post_process_functions.general_functions import *


def compute_access_to_target(access_fname, target_illumination_fname, satellite_illumination_fname, sim_start_time,
                             **kwargs):
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
    target_eclipse_times = get_illumination_event_times_from_file(target_illumination_fname, invert=True)
    check_event_order_consistency(target_eclipse_times)

    # Filter out times when satellite is not illuminated
    filter_satellite_illumination = kwargs.pop("filter_satellite_illumination", False)
    if filter_satellite_illumination:
        access_times = get_overlap_between_events(access_times, satellite_illumination_times)
    check_event_order_consistency(access_times)

    # Get times when target power is available
    target_illumination_times = combine_events(access_times, target_illumination_times)
    check_event_order_consistency(target_illumination_times)

    # Filter out times when power is not needed by target
    filter_target_eclipses = kwargs.pop("filter_target_illumination", False)
    if filter_target_eclipses:
        access_times = get_overlap_between_events(access_times, target_eclipse_times)
    check_event_order_consistency(access_times)

    access_durations = list()
    for i, access_start in enumerate(access_times["Start"]):
        access_end = access_times["End"][i]
        access_duration = (access_end - access_start).total_seconds()
        access_durations.append(access_duration)
    last_eclipse_start = sim_start_time
    eclipse_durations = list()
    for j, target_power_start in enumerate(target_illumination_times["Start"]):
        target_power_end = target_illumination_times["End"][j]
        eclipse_time = (target_power_start - last_eclipse_start).total_seconds()
        eclipse_durations.append(eclipse_time)
        last_eclipse_start = target_power_end

    process_results(access_durations, eclipse_durations)


def process_results(access_durations, eclipse_durations):
    """
    Process the results of analysis from computing the access and eclipse times from a simulation

    access_durations: list of time of each access period
    eclipse_durations: list of time of each eclipse period
    :return:
    """
    # Counts to make sure count the number of long eclipses - this is to compare against the statement in
    # the Brandhorst paper
    long_eclipse_time = 84.0 * 3600.0
    long_eclipse_count = 0
    for eclipse in eclipse_durations:
        long_eclipse_count = long_eclipse_count + 1 if eclipse > long_eclipse_time else long_eclipse_count

    seconds_to_hours = 1 / 3600.0

    # Convert results to arrays
    access_durations = np.asarray(access_durations) * seconds_to_hours
    eclipse_durations = np.asarray(eclipse_durations) * seconds_to_hours

    # Get salient statistics
    total = np.sum(access_durations)
    average = np.average(access_durations)
    standard_deviation = np.std(access_durations)
    maximum = np.max(access_durations)
    minimum = np.min(access_durations)

    print("Number of Eclipses longer than 84hrs: {}\n".format(long_eclipse_count))

    print("Total Access Time: {} hours".format(total))
    print("Average Access Time: {} hours".format(average))
    print("Standard deviation Access Time: {} hours".format(standard_deviation))
    print("Maximum Access Time: {} hours".format(maximum))
    print("Minimum Access Time: {} hours\n".format(minimum))

    total = np.sum(eclipse_durations)
    average = np.average(eclipse_durations)
    standard_deviation = np.std(eclipse_durations)
    maximum = np.max(eclipse_durations)
    minimum = np.min(eclipse_durations)

    print("Total Eclipse Time: {} hours".format(total))
    print("Average Eclipse Time: {} hours".format(average))
    print("Standard deviation Eclipse Time: {} hours".format(standard_deviation))
    print("Maximum Eclipse Time: {} hours".format(maximum))
    print("Minimum Eclipse Time: {} hours\n".format(minimum))

    # Plot histograms
    fig, ax = plt.subplots(2)

    ax[0].hist(access_durations, bins=200)
    ax[0].set_xlabel("Frequency")
    ax[0].set_ylabel("Time (Hours)")
    ax[0].set_title("Access Time Histogram")

    ax[1].hist(eclipse_durations, bins=200)
    ax[1].set_xlabel("Frequency")
    ax[1].set_ylabel("Time (Hours)")
    ax[1].set_title("Eclipse Time Histogram")

    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots(2)

    ax[0].plot(access_durations)
    ax[0].set_xlabel("Event number")
    ax[0].set_ylabel("Time (Hours)")
    ax[0].set_title("Access Times")

    ax[1].plot(eclipse_durations)
    ax[1].set_xlabel("Event number")
    ax[1].set_ylabel("Time (Hours)")
    ax[1].set_title("Eclipse Times")

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    access_fname = "Brandhorst_Sim_1_Low_Res.csv"
    target_illumination_fname = "Target IlluminationTimes.csv"
    satellite_illumination_fname = "Satellite IlluminationTimes.csv"
    sim_start_time = convert_string_to_datetime(['2008', '07', '01', '11', '0', '0.0'])

    compute_access_to_target(access_fname, target_illumination_fname, satellite_illumination_fname, sim_start_time,
                             filter_target_illumination=True, filter_satellite_illumination=True)