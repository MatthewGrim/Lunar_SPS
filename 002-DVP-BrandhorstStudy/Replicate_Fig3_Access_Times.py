""""
Author: Darian van Paridon
Date: 04/05/2018

This file processes csv output data from STK for evaluating the access times of a solar power satellite
to a lunar base at 45N degree latitude, for comparison to the Brandhorst paper.
"""

import numpy as np
import matplotlib.pyplot as plt
from general_functions import *


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def check_sum(lighting, eclipse, duration):
    # A quick check to ensure that every time stamp from the STK simulation is accounted for
    # and belongs to an event. Threshold for "no error" is set to 30 seconds.
    # Used to check for errors when inverting an event set (ie. find eclipses from sunlit events)
    check_sum = (np.sum(lighting[2]) + np.sum(eclipse[2]))
    if abs(check_sum - duration) > 30:
        print('Inverted and original event times do not add to total duration!')
        print('{} seconds of total time unaccounted for'.format((duration - check_sum)))


def parse_csv_to_array(file_name, sim_start):
    # This function takes in the .csv file from STK, and converts it into a 3 X N array, where N is the number of events
    # in the file. The three data types are beginning of event, end of event (both measured in seconds from the
    # beginning of the simulation), and duration of event. The .csv files are edited in Excel to remove an additional
    # information which cause errors when running them through this script.
    # print("Starting to parse script...")
    file = open(file_name, "r")
    size = file_len(file_name)
    duration_array = np.zeros(size)
    start_time_sec_from_simstart = np.zeros(size)
    end_time_sec_from_simstart = np.zeros(size)
    parsed_data = np.array((3, size))
    # print("Number of events = {}".format(size))

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
    # "Conditional events" should be times when the SPS would be active, if it can access the target
    # Therefore take satellite illumination and target eclipses as event_times.
    # Output is an array containing start, end and duration of overlap events.
    overlap_start = list()
    overlap_end = list()
    overlap_duration = list()

    access_size = len(access_times[0])
    events_size = len(event_times[0])

    # Index for access periods
    i = 0
    # Test all of the access periods for overlap with conditional events
    while i < len(access_times[0]):
        # Get start and end of access periods
        access_start = access_times[0][i]
        access_end = access_times[1][i]

        # Check each conditional event for overlap
        for j in range(events_size):
            event_start = event_times[0][j]
            event_end = event_times[1][j]

            # Determine overlapping events
            # Partial overlap, event period triggers overlap, access period ends it
            if access_start <= event_start and access_end <= event_end and event_start <= access_end:
                overlap_start.append(event_start)
                overlap_end.append(access_end)
                overlap_duration.append(access_end - event_start)

            # Partial overlap, access period triggers overlap, event period ends it
            elif event_start <= access_start and event_end <= access_end and access_start <= event_end:
                overlap_start.append(access_start)
                overlap_end.append(event_end)
                overlap_duration.append(event_end - access_start)

            # Full overlap, access period occurs within event period
            elif event_start <= access_start and access_end <= event_end:
                overlap_start.append(access_start)
                overlap_end.append(access_end)
                overlap_duration.append(access_end - access_start)

            # Full overlap, event period occurs within access period
            elif access_start <= event_start and event_end <= access_end:
                overlap_start.append(event_start)
                overlap_end.append(event_end)
                overlap_duration.append(event_end - event_start)

            # No overlap, access period ends before event
            elif access_end < event_start:
                pass

            # No overlap, access period begins after event
            elif event_end < access_start:
                pass

            else:
                raise RuntimeError("Should not be possible to get here!")
        i += 1

    overlap_events = (overlap_start, overlap_end, overlap_duration)
    return overlap_events


def determine_SPS_active_time(sunlight_SPS, eclipse_target, access_times):
    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.
    total_availability = np.sum(access_times[2])
    print("Total possible access time: {} hrs".format(round(total_availability / 3600.0, 2)))

    target_eclipse_during_access = get_event_overlaps(access_times, eclipse_target)

    target_eclipse_SPS_sunlit_during_access = get_event_overlaps(target_eclipse_during_access, sunlight_SPS)
    total_SPS_time = np.sum(target_eclipse_SPS_sunlit_during_access[2])
    print("Filtering out eclipses, total active time: {} hrs".format(round(total_SPS_time / 3600.0, 2)))
    return total_SPS_time


def get_times_without_coverage(access_times):
    length = len(access_times[0]) - 1
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
    # If access period begins at start of simulation
    elif access_times[0][0] == 0:
        for j in range(length - 2):
            no_coverage[0][j] = access_times[1][j]
            no_coverage[1][j] = access_times[0][j+1]
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]

    return no_coverage


def determine_blackout_data(access_times, eclipse_target,  duration):

    # This function replicates the results of Figure 2 in Brandhorst's paper
    # regarding black-out periods for the lunar target (ie. times when there is
    # no SPS coverage, and the target is in eclipse.

    # Calculates the times when the SPS cannot access the target
    SPS_no_coverage = get_times_without_coverage(access_times)
    check_sum(access_times, SPS_no_coverage, duration)
    # Calculates overlap of no coverage with target eclipses
    dark_periods = get_event_overlaps(SPS_no_coverage, eclipse_target)
    # Calculate number of eclipses exceeding 84 hours, and the maximum eclipse duration
    temp_array = np.array(dark_periods[2])
    long_eclipse_flag = (temp_array / 3600.0) > 84.0
    num_long_eclipse = np.sum(long_eclipse_flag)
    max_eclipse_duration = round(max([i/3600.0 for i in dark_periods[2]]), 2)
    print('Maximum black-out duration: {} hrs'.format(max_eclipse_duration))
    print("Number of blackouts exceeding 84 hrs: {}".format(num_long_eclipse))
    print("\n")

    # Plot which shows start time and duration of periods when target gets no power
    plt.figure(1)
    plt.subplot(211)
    plt.bar(SPS_no_coverage[0] / 84600.0, SPS_no_coverage[2] / 3600.0)
    plt.ylabel('Duration of No-SPS-Coverage Period')
    plt.title('Brandhorst Figure 2 Remake')
    plt.subplot(212)
    plt.bar([i / 86400.0 for i in dark_periods[0]], [i / 3600.0 for i in dark_periods[2]])
    plt.xlabel('Time Since Start of Simulation (July 1, 2008) [days]')
    plt.ylabel('Duration of Black-Out Period')
    plt.show()


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
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_target, access_times)

    # Plot which shows the variation in total access time
    # as a function of relative argument of perigee between SPS2 and SPS1
    plt.figure(3)
    plt.plot(relative_perigees, total_active_time / 3600.0)
    plt.xlabel('Relative Argument of Perigee [deg]')
    plt.ylabel('Total Access Time of SPS2 [hrs]')
    plt.title("Brandhorst Figure 3 Remake")
    plt.show()


def main():
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
    print("\n")
    print("Simulation duration in days:")
    print(total_duration / 86400.0)

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS
    # Solar Power Satellite 1
    LOS_access1 = parse_csv_to_array(raw_access_SPS1, start)
    sunlight_SPS1 = parse_csv_to_array(raw_sunlight_SPS1, start)
    eclipse_SPS1 = parse_csv_to_array(raw_eclipse_SPS1, start)
    check_sum(sunlight_SPS1, eclipse_SPS1, total_duration)
    # Solar Power Satellite 2
    LOS_access2 = parse_csv_to_array(raw_access_SPS2, start)
    sunlight_SPS2 = parse_csv_to_array(raw_sunlight_SPS2, start)
    eclipse_SPS2 = parse_csv_to_array(raw_eclipse_SPS2, start)
    check_sum(sunlight_SPS2, eclipse_SPS2, total_duration)
    # Lunar Target
    sunlight_target = parse_csv_to_array(raw_sunlight_target, start)
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)
    check_sum(sunlight_target, eclipse_target, total_duration)

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

    # Plot which shows all access period, SPS eclipses, and target illumination periods.
    # plt.figure(2)
    # plt.subplot(311)
    # plt.bar(LOS_access2[0] / 86400.0, LOS_access2[2] / 3600.0)
    # plt.ylabel('Line-Of-Sight Access Duration [h]')
    # plt.title('Illumination and Access Times for SPS and Lunar Target at 45N')
    # plt.subplot(312)
    # plt.bar(eclipse_SPS2[0] / 86400.0, eclipse_SPS2[2] / 3600.0)
    # plt.ylabel('SPS2 Eclipse Event Duration [h]')
    # plt.subplot(313)
    # plt.bar(sunlight_target[0] / 86400.0, sunlight_target[2] / 3600.0)
    # plt.ylabel('Target Sunlit Event Duration [h]')
    # plt.xlabel('Time Since Start of Simulation (July 1, 2008) [days]')
    # plt.show()

    # Calculates the total access time for various relative arguments of perigee
    scan_perigee_angles(start, total_duration, eclipse_target)


main()
