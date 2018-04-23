"""
Author: Rohan Ramasamy
Date: 04/03/2018

This is a data parsing script for reeding in a csv with access times from the lunar sps study.
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime


def get_datetime(time_components):
    year = int(time_components[0])
    month = int(time_components[1])
    day = int(time_components[2])
    hour = int(time_components[3])
    minute = int(time_components[4])
    seconds = int(float(time_components[5]))

    t = datetime.datetime(year, month, day, hour, minute, seconds)

    return t


def parse_csv_to_array(file_name, sim_start):
    """[summary]
    
    Arguments:
        file_name {string} -- name of the .csv file to be parsed
    """
    print("Starting to parse script...")

    eclipse_start_time = get_datetime(sim_start)
    access_durations = list()
    eclipse_durations = list()
    file = open(file_name, "r")
    for i, line in enumerate(file):
        if i == 0:
            continue

        # Split line into table components, assuming comma delimited
        components = line.split(",")

        # For the last line, break it table elements do not exist
        if len(components) < 4:
            break

        # Work out eclipse time and set new time t
        start = components[1].split(":")
        start_time = get_datetime(start)
        eclipse_time = (start_time - eclipse_start_time).total_seconds()
        end = components[2].split(":")
        eclipse_start_time = get_datetime(end)

        # Save power link duration and eclipse durations
        eclipse_durations.append(eclipse_time)
        access_durations.append(float(components[3].split("\n")[0]))

    # Convert results to arrays
    seconds_to_hours = 1 / 3600.0
    access_durations = np.asarray(access_durations) * seconds_to_hours
    eclipse_durations = np.asarray(eclipse_durations) * seconds_to_hours

    # Get salient statistics
    average = np.average(access_durations)
    standard_deviation = np.std(access_durations)
    maximum = np.max(access_durations)
    minimum = np.min(access_durations)

    print("Average Access Time: {} hours".format(average))
    print("Standard deviation Access Time: {} hours".format(standard_deviation))
    print("Maximum Access Time: {} hours".format(maximum))
    print("Minimum Access Time: {} hours\n".format(minimum))

    average = np.average(eclipse_durations)
    standard_deviation = np.std(eclipse_durations)
    maximum = np.max(eclipse_durations)
    minimum = np.min(eclipse_durations)

    print("Average Eclipse Time: {} hours".format(average))
    print("Standard deviation Eclipse Time: {} hours".format(standard_deviation))
    print("Maximum Eclipse Time: {} hours".format(maximum))
    print("Minimum Eclipse Time: {} hours\n".format(minimum))

    # Plot histograms
    fig, ax = plt.subplots(2)
    ax[0].hist(access_durations, bins=200)
    ax[1].hist(eclipse_durations, bins=200)
    plt.show()

    fig, ax = plt.subplots(2)
    ax[0].plot(access_durations)
    ax[1].plot(eclipse_durations)
    plt.show()


if __name__ == '__main__':
    # file_name = "Target-Target1-Sensor-Sensor1-To-Satellite-Satellite1 Access.csv"
    file_name = "high_res.csv"
    start = ['2017', '03', '01', '11', '0', '0.0']
    parse_csv_to_array(file_name, start)

