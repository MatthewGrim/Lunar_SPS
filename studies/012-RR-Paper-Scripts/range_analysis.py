"""
Author: Rohan Ramasamy
Date: 14/09/2018

This script is used to consider the range data for the optimised orbit in more detail
"""

import os
import numpy as np
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import convert_string_to_datetime, file_len


def import_range_data(file_name, sim_start=convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])):
    processed_file_name = "{}_processed".format(file_name)
    file_name = "{}.csv".format(file_name)
    if not os.path.exists(processed_file_name):
        the_file = open(file_name, "r")
        size = file_len(file_name)
        range_data = np.zeros((4, size - 2))

        # For loop parses data into three categories: start, end and duration of event
        # This function can be used for target to satellite access times, or target/satellite illumination events
        # .csv file must be three columns, with start time, end time, duration in that order
        # The function outputs the start and end times of the event in seconds since the beginning of the simulation
        for i, line in enumerate(the_file):
            if i == 0:
                continue
            # Split line into table components, assuming comma delimited
            components = line.split(",")
            # Break once the final line is reached (which is a blank 'return')
            if components[0] == '\n':
                break

            # Work out sunlit times and set new time t
            start = components[0].split(":")
            start_time = convert_string_to_datetime(start)
            range_data[0, i - 1] = (start_time - sim_start).total_seconds()
            range_data[1, i - 1] = components[1]
            range_data[2, i - 1] = components[2]
            range_data[3, i - 1] = components[3]
        range_data[0, :] /= 3600 * 24
        np.savetxt(processed_file_name, range_data)
    else:
        range_data = np.loadtxt(processed_file_name)

    return range_data


def power_delivered(range_data, laser_power, laser_radius, laser_wavelength, receiver_area):
    beam_radius = laser_radius * np.sqrt(1.0 + (range_data[1, :] * 1e3 * laser_wavelength / (np.pi * laser_radius ** 2)) ** 2)
    beam_area = np.pi * beam_radius ** 2
    receiver_power = laser_power * receiver_area / beam_area

    return receiver_power


def process_data(rover):
    if rover == "sorato":
        data = import_range_data("2500_2500_range_data")
        power = power_delivered(data, 29.77e3, 0.317, 1070e-9, 0.175)
    elif rover == "AMALIA":
        data = import_range_data("1700_1700_range_data")
        power = power_delivered(data, 19.43e3, 0.2392, 850e-9, 0.366)
    else:
        raise ValueError("Rover name does not exist!")

    receiver_efficiency = 0.5
    power = receiver_efficiency * power[:]
    plt.figure()
    plt.plot(data[0, :], power[:])
    plt.xlim((0, 1))
    plt.show()


if __name__ == '__main__':
    rover_name = "AMALIA"
    process_data(rover_name)
