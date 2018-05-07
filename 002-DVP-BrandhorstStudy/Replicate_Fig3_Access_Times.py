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
    print("Number of events = ", size)

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


def get_eclipse_data(parsed_data):

    # This function takes in illumination event data (start, end, duration) and determines the same information
    # for eclipse events, assuming that penumbra eclipse is still eclipse, and that the simulation ends
    # with the satellite sunlit (not in eclipse)
    size = len(parsed_data[0])

    eclipse_duration = np.zeros(size)
    eclipse_start = np.zeros(size)
    eclipse_end = np.zeros(size)

    for i in range(0, size-1):
        eclipse_duration[i] = parsed_data[0][i+1] - parsed_data[1][i]
        eclipse_start[i] = parsed_data[1][i]
        eclipse_end[i] = parsed_data[0][i+1]

    eclipse_duration[size-2] = parsed_data[0][size-2] - parsed_data[1][size-3]
    eclipse_start[size-2] = parsed_data[1][size-3]
    eclipse_end[size-2] = parsed_data[0][size-2]

    eclipse_data = [eclipse_start, eclipse_end, eclipse_duration]

    # return an 3 x N array for N eclipse events
    return eclipse_data


def determine_SPS_active_time(sunlight_SPS, sunlight_target, LOS_access):

    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.

    eclipse_SPS = get_eclipse_data(sunlight_SPS)

    SPS_eclipse_end_flag = np.zeros(len(eclipse_SPS[0]))
    SPS_eclipse_start_flag = np.zeros(len(eclipse_SPS[0]))
    eclipse_time_during_access = np.zeros(len(eclipse_SPS[0]))

    target_sunlight_end_flag = np.zeros(len(sunlight_target[0]))
    target_sunlight_start_flag = np.zeros(len(sunlight_target[0]))
    target_illumination_during_access = np.zeros(len(sunlight_target[0]))

    # Check all access periods to see if the SPS would be required
    for i in range(0, len(LOS_access[0])):

        # Check to see if the target is illuminated by the sun
        for j in range(0, len(sunlight_target[0])):
            target_sunlight_start_flag[j] = LOS_access[0][i] <= sunlight_target[0][j] <= LOS_access[1][i]
            target_sunlight_end_flag[j] = LOS_access[0][i] <= sunlight_target[1][j] <= LOS_access[1][i]
            # Check to see when the target is illuminated by the sun during access periods
            # Case where entire sunlit period occurs during access period
            if target_sunlight_start_flag[j] == 1 and target_sunlight_end_flag[j] == 1:
                target_illumination_during_access[j] = sunlight_target[2][j]
            # Case where eclipse begins during access, but access period ends before eclipse
            elif target_sunlight_start_flag[j] == 1 and target_sunlight_end_flag[j] == 0:
                target_illumination_during_access[j] = LOS_access[1][i] - sunlight_target[0][j]
            # Case where access period begins during eclipse, and eclipse ends before access period
            elif target_sunlight_start_flag[j] == 0 and target_sunlight_end_flag[j] == 1:
                target_illumination_during_access[j] = sunlight_target[1][j] - LOS_access[0][i]

        total_target_illumination_during_access = np.sum(target_illumination_during_access) / 3600.0

        # Check to see if the SPS would be in eclipse
        for j in range(0, len(eclipse_SPS[0])):
            SPS_eclipse_start_flag[j] = LOS_access[0][i] <= eclipse_SPS[0][j] <= LOS_access[1][i]
            SPS_eclipse_end_flag[j] = LOS_access[0][i] <= eclipse_SPS[1][j] <= LOS_access[1][i]
            # Check to see when satellite is in eclipse during access periods
            # Case where entire eclipse occurs during access period
            if SPS_eclipse_start_flag[j] == 1 and SPS_eclipse_end_flag[j] == 1:
                eclipse_time_during_access[j] = eclipse_SPS[2][j]
            # Case where eclipse begins during access, but access period ends before eclipse
            elif SPS_eclipse_start_flag[j] == 1 and SPS_eclipse_end_flag[j] == 0:
                eclipse_time_during_access[j] = LOS_access[1][i] - eclipse_SPS[0][j]
            # Case where access period begins during eclipse, and eclipse ends before access period
            elif SPS_eclipse_start_flag[j] == 0 and SPS_eclipse_end_flag[j] == 1:
                eclipse_time_during_access[j] = eclipse_SPS[1][j] - LOS_access[0][i]

        total_SPS_eclipse_during_access = np.sum(eclipse_time_during_access) / 3600.0

    print("Total time SPS is eclipsed during access period [hrs]: ")
    print(total_SPS_eclipse_during_access)
    print("Total time target is illuminated during access periods [hrs]: ")
    print(total_target_illumination_during_access)

    total_SPS_time = (sum(LOS_access[2]) / 3600.0) - total_target_illumination_during_access - total_SPS_eclipse_during_access

    print("Total active time for SPS [hrs]: ")
    print(total_SPS_time)

def main():
    # Import data and set the start and end times of the simulation
    sunlight_data_SPS = "SPS2-Lighting-Times-Edited.csv"
    access_data_SPS = "SPS2-Access-135-Edited.csv"
    sunlight_data_target = 'Target1-Lighting-Times-Edited.csv'
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    print("Simulation duration, in days:")
    print(total_duration / 86400)

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS
    sunlight_SPS = parse_csv_to_array(sunlight_data_SPS, start)
    LOS_access = parse_csv_to_array(access_data_SPS, start)
    sunlight_target = parse_csv_to_array(sunlight_data_target, start)

    eclipse_SPS = get_eclipse_data(sunlight_SPS)

    plt.figure(1)
    plt.subplot(311)
    plt.bar(LOS_access[0] / 86400, LOS_access[2] / 3600)
    plt.ylabel('Line-Of-Sight Access Duration [h]')
    plt.title('Illumination and Access Times for SPS and Lunar Target at 45N')
    plt.subplot(312)
    plt.bar(eclipse_SPS[0] / 86400, sunlight_SPS[2] / 3600)
    plt.ylabel('SPS Eclipse Event Duration [h]')
    plt.subplot(313)
    plt.bar(sunlight_target[0] / 86400, sunlight_target[2] / 3600)
    plt.ylabel('Target Sunlit Event Duration [h]')
    plt.xlabel('Time Since Start of Simulation (July 1, 2008) [days]')
    plt.show()

    total_availability = np.sum(LOS_access[2])

    print("Total possible access time [hrs]: ")
    print(total_availability / 3600.0)

    determine_SPS_active_time(sunlight_SPS, sunlight_target, LOS_access)


main()
