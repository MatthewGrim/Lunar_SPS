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
        laser_power = 4.29e3
        laser_radius = 1.0562
        wavelength = 1070e-9
        receiver_area = 0.07863935625457205
        data = import_range_data(os.path.join("2300", "2300_2300_range_data"))
    elif rover == "AMALIA_1700":
        laser_power = 19.43e3
        laser_radius = 0.2392
        wavelength = 850e-9
        receiver_area = 0.366
        data = import_range_data(os.path.join("1700", "1700_1700_range_data"))
    elif rover == "AMALIA_1300":
        laser_power = 2.98e3
        laser_radius = 0.8844
        wavelength = 1070e-9
        receiver_area = 0.366
        data = import_range_data(os.path.join("1300", "1300_1300_range_data"))
    elif rover == "AMALIA_1300_low_eff":
        laser_power = 7.79e3 * (1 - np.exp(-2))
        laser_radius = 0.2985
        wavelength = 1070e-9
        receiver_area = 1.83
        data = import_range_data(os.path.join("1300", "1300_1300_range_data"))
    elif rover == "AMALIA_800":
        laser_power = 9.61e3 * (1 - 1 / np.exp(2))
        laser_radius = 0.2888
        wavelength = 1070e-9
        receiver_area = 0.366
        data = import_range_data(os.path.join("800", "800_800_range_data"))
    else:
        raise ValueError("Rover name does not exist!")

    power = power_delivered(data, laser_power, laser_radius, wavelength, receiver_area)

    power_density = power / receiver_area
    receiver_efficiency = 0.5
    power = receiver_efficiency * power[:]
    fig, ax = plt.subplots(2, figsize=(12, 7), sharex=True)
    ax[0].plot(data[0, :], power_density[:])
    ax[0].set_xlim([data[0, 0], 1])
    ax[0].set_ylabel("Surface flux [$Wm^{-2}$]")
    ax[1].plot(data[0, :], power[:])
    ax[1].set_ylabel("Received power [$W$]")
    ax[1].set_xlim([data[0, 0], 1])
    fig.suptitle("Power received plots for {}".format(rover))
    plt.savefig("{}_power_delivered".format(rover))
    plt.show()


if __name__ == '__main__':
    rover_name = "AMALIA_1300_low_eff"
    process_data(rover_name)
