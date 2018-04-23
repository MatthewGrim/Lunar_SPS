"""
Author: Rohan Ramasamy
Date: 04/03/2018

This is a data parsing script for reeding in a csv with access times from the lunar sps study.
"""

import numpy as np
import matplotlib.pyplot as plt


def parse_csv_to_array(file_name):
    """[summary]
    
    Arguments:
        file_name {string} -- name of the .csv file to be parsed
    """
    print("Starting to parse script...")

    durations = list()
    file = open(file_name, "r")
    for i, line in enumerate(file):
        if i == 0:
            continue
        components = line.split(",")
        durations.append(float(components[3].split("\n")[0]))
    
    durations = np.asarray(durations)
    print(durations.shape)

    # Get salient statistics    
    seconds_to_hours = 1 / 3600.0
    average = np.average(durations) * seconds_to_hours
    standard_deviation = np.std(durations) * seconds_to_hours
    maximum = np.max(durations) * seconds_to_hours
    minimum = np.min(durations) * seconds_to_hours

    print("Average Access Time: {} hours".format(average))
    print("Average Access Time: {} hours".format(standard_deviation))
    print("Maximum Access Time: {} hours".format(maximum))
    print("Minimum Access Time: {} hours".format(minimum))

    plt.figure()
    plt.hist(durations, bins=500)
    plt.show()

if __name__ == '__main__':
    parse_csv_to_array("Target-North-Sensor-Sensor1-To-Satellite-Satellite1 Access.csv")

