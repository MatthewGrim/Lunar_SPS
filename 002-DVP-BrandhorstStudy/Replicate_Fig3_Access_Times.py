""""
Author: Darian van Paridon
Date: 04/05/2018

This file processes csv output data from STK for evaluating the access times of a solar power satellite
to a lunar base at 45N degree latitude, for comparison to the Brandhorts paper.
"""

import numpy as np
import matplotlib.pyplot as plt
from general_functions import *


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def parse_csv_to_array(file_name, sim_start):
    """[summary]

    Arguments:
        file_name {string} -- name of the .csv file to be parsed

    This function takes in the .csv file from STK, and converts it into a 3 X N array, where N is the number of events
    in the file. The three data types are beginning of event, end of event (both measured in seconds from the
    beginning of the simulation), and duration of event. The .csv files are edited in Excel to remove an additional
    information which cause errors when running them through this script.
    """
    print("Starting to parse script...")

    file = open(file_name, "r")
    size = file_len(file_name)
    duration_array = np.zeros(size)
    start_time_sec_from_simstart = np.zeros(size)
    end_time_sec_from_simstart = np.zeros(size)
    parsed_data = np.array((3, size))
    print("Number of events = {}".format(size))

    # For loop parses data into three categories: start, end and duration of event
    # This function can be used for target to satellite access times, or target/satellite illumination events
    # .csv file must be three columns, with start time, end time, duration in that order
    # The function outputs the start and end times of the event in seconds since the beginning of the simulation
    for i, line in enumerate(file):
        if i == 0:
            continue

        # Split line into table components, assuming comma delimited
        components = line.split(",")

        # Work out sunlit times and set new time t
        start = components[0].split(":")
        start_time = convert_string_to_datetime(start)
        start_time_sec_from_simstart[i-1] = (start_time - sim_start).total_seconds()

        end = components[1].split(":")
        end_time = convert_string_to_datetime(end)
        end_time_sec_from_simstart[i-1] = (end_time - sim_start).total_seconds()

        # Save power link duration
        access_durations = float(components[2].split("\n")[0])
        duration_array[i-1] = access_durations

        parsed_data = [start_time_sec_from_simstart, end_time_sec_from_simstart, duration_array]

    return parsed_data


def get_event_overlaps(access_times, event_times):

    # This function finds overlap between events and access times between the SPS and the target.
    # "Events" should be times when the SPS would be active.
    # Therefore take satellite illumination and target eclipses.
    # Output is an array containing start, end and duration of overlap events.

    event_end_flag = np.zeros(len(event_times[0]))
    event_start_flag = np.zeros(len(event_times[0]))
    overlap_events = np.zeros((3, len(event_times[0]) + 1))

    # Iterate through all of the access period events
    for i in range(len(access_times[0])):
        # Flag all events which overlap with access periods
        for j in range(len(event_times[0])):
            event_start_flag[j] = access_times[0][i] <= event_times[0][j] <= access_times[1][i]
            event_end_flag[j] = access_times[0][i] <= event_times[1][j] <= access_times[1][i]
            # Find start, end and duration of overlap event
            # Case where entire event occurs during access period
            if event_start_flag[j] == 1 and event_end_flag[j] == 1:
                overlap_events[0][j] = event_times[0][j]
                overlap_events[1][j] = event_times[1][j]
                overlap_events[2][j] = event_times[2][j]
            # Case where event begins during access, but access period ends first
            elif event_start_flag[j] == 1 and event_end_flag[j] == 0:
                overlap_events[0][j] = event_times[0][j]
                overlap_events[1][j] = access_times[1][i]
                overlap_events[2][j] = access_times[1][i] - event_times[0][j]
            # Case where access period begins during event, but event ends before access period
            elif event_start_flag[j] == 0 and event_end_flag[j] == 1:
                overlap_events[0][j] = access_times[0][i]
                overlap_events[1][j] = event_times[1][j]
                overlap_events[2][j] = event_times[1][j] - access_times[0][i]

    # plt.subplot(311)
    # plt.bar(access_times[0] / 86400.0, access_times[2] / 3600.0)
    # plt.title('Overlapping Events Which Define SPS Active Time')
    # plt.ylabel('Access Times')
    # plt.subplot(312)
    # plt.bar(event_times[0] / 86400.0, event_times[2] / 3600.0)
    # plt.ylabel('Event Times, Condition For Active SPS)')
    # plt.subplot(313)
    # plt.bar(overlap_events[0] / 86400.0, overlap_events[2] / 3600.0)
    # plt.ylabel('Overlap of Events and Access')
    # plt.show()

    return overlap_events


def determine_SPS_active_time(sunlight_SPS, eclipse_target, access_times):

    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.

    total_availability = np.sum(access_times[2])

    print("Total possible access time based on SPS orbit and target location: {} hrs".format(total_availability / 3600.0))

    SPS_sunlit_during_access = get_event_overlaps(access_times, sunlight_SPS)

    target_eclipse_SPS_sunlit_during_access = get_event_overlaps(SPS_sunlit_during_access, eclipse_target)

    total_SPS_time = np.sum(target_eclipse_SPS_sunlit_during_access[2])

    print("Excluding times when SPS is eclipsed and/or when target is sunlit, total access time for SPS: {} hrs".format(total_SPS_time / 3600.0))


def find_long_eclipse(eclipse_SPS, eclipse_target):

    double_eclipse = get_event_overlaps(eclipse_target, eclipse_SPS)

    long_eclipse_flag = (double_eclipse[2] / 3600.0) > 84.0

    num_long_eclipse = sum(long_eclipse_flag)

    max_eclipse_duration = max(double_eclipse[2]) / 3600.0

    print('Maximum dark time: {} hrs'.format(max_eclipse_duration))

    print("Number of eclipses exceeding 84 hrs: {}".format(num_long_eclipse))

def main():
    # Import data and set the start and end times of the simulation
    raw_sunlight_SPS2 = "SPS2-Lighting-Times-Edited.csv"
    raw_eclipse_SPS2 = "SPS2-Eclipse-Times-Edited.csv"
    raw_access_SPS2 = "SPS2-Access-135-Edited.csv"
    raw_sunlight_SPS1 = "SPS1-Lighting-Times-Edited.csv"
    raw_eclipse_SPS1 = "SPS1-Lighting-Times-Edited.csv"
    raw_access_SPS1 = "SPS1-Access-0-Edited.csv"
    raw_sunlight_target = 'Target1-Lighting-Times-Edited.csv'
    raw_eclipse_target = 'Target1-Eclipse-Times-Edited.csv'
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    print("Simulation duration, in days:")
    print(total_duration / 86400.0)

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS
    sunlight_SPS2 = parse_csv_to_array(raw_sunlight_SPS2, start)
    eclipse_SPS2 = parse_csv_to_array(raw_eclipse_SPS2, start)

    check_sum = (np.sum(sunlight_SPS2[2]) + np.sum(eclipse_SPS2[2])) / (365.0*24.0*3600.0)
    if check_sum != 2:
        print('Missing time between sunlight and eclipse')

    LOS_access2 = parse_csv_to_array(raw_access_SPS2, start)
    sunlight_SPS1 = parse_csv_to_array(raw_sunlight_SPS1, start)
    eclipse_SPS1 = parse_csv_to_array(raw_eclipse_SPS1, start)
    LOS_access1 = parse_csv_to_array(raw_access_SPS1, start)
    sunlight_target = parse_csv_to_array(raw_sunlight_target, start)
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)

    # Plot which shows all access period, SPS eclipses, and target illumination periods.
    plt.figure(1)
    plt.subplot(311)
    plt.bar(LOS_access2[0] / 86400.0, LOS_access2[2] / 3600.0)
    plt.ylabel('Line-Of-Sight Access Duration [h]')
    plt.title('Illumination and Access Times for SPS and Lunar Target at 45N')
    plt.subplot(312)
    plt.bar(eclipse_SPS2[0] / 86400.0, eclipse_SPS2[2] / 3600.0)
    plt.ylabel('SPS2 Eclipse Event Duration [h]')
    plt.subplot(313)
    plt.bar(sunlight_target[0] / 86400.0, sunlight_target[2] / 3600.0)
    plt.ylabel('Target Sunlit Event Duration [h]')
    plt.xlabel('Time Since Start of Simulation (July 1, 2008) [days]')
    plt.show()

    print("\n")
    print("Access availability for SPS1")
    determine_SPS_active_time(sunlight_SPS1, eclipse_target, LOS_access1)
    find_long_eclipse(eclipse_SPS1, eclipse_target)

    print("\n")
    print("Access availability for SPS2")
    determine_SPS_active_time(sunlight_SPS2, eclipse_target, LOS_access2)
    find_long_eclipse(eclipse_SPS2, eclipse_target)


main()
