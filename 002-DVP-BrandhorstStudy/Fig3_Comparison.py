""""
11/05/2018
Author: Darian van Paridon

This script compares the results of the the replication of the Brandhorst figure 3 for SPS access time as a function
of argument of perigee to the results published in the paper
"""

from Replicate_Fig3_Access_Times import *


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx


def parse_fig3_data(file_name):

    the_file = open(file_name, "r")
    size = file_len(file_name)
    perigees = np.zeros(size)
    total_access_time = np.zeros(size)
    parsed_data = np.array((2, size))

    # For loop parses data into three categories: start, end and duration of event
    # This function can be used for target to satellite access times, or target/satellite illumination events
    # .csv file must be three columns, with start time, end time, duration in that order
    # The function outputs the start and end times of the event in seconds since the beginning of the simulation
    for i, line in enumerate(the_file):
        # Split line into table components, assuming comma delimited
        components = line.split(",")

        # Copy data into array
        perigees[i] = components[0]
        total_access_time[i] = components[1]

    parsed_data = [perigees, total_access_time]

    return parsed_data


def main():
    # Process the STK data, import Brandhorst data, initialize array for comparison of data
    perigees, active_times = remake_brandhorst_fig3()
    brandhorst_data = parse_fig3_data("Brandhorst_Fig3_Data.csv")

    # i = 0
    # for i in range(len(brandhorst_data[0])):
    #     if brandhorst_data[0][i] < 0:
    #         brandhorst_data[0][i] += 360
    #     i += 1

    brandhorst_data[0] += abs(min(brandhorst_data[0]))

    nearest_brandhorst_perigee = np.zeros(len(perigees))
    brandhorst_active_time = np.zeros(len(perigees))
    delta_active_time = np.zeros(len(perigees))
    percent_difference = np.zeros(len(perigees))
    active_times_flip = active_times[::-1]

    for i in range(len(perigees)):
        nearest_brandhorst_perigee[i], match_index = find_nearest(brandhorst_data[0], perigees[i])
        brandhorst_active_time[i] = brandhorst_data[1][match_index]
        delta_active_time[i] = (active_times_flip[i] / 3600.0) - brandhorst_active_time[i]
        percent_difference[i] = (delta_active_time[i] / brandhorst_active_time[i]) * 100.0

    print("\n")
    print("Average error in the remake data: {} hours ".format(round(np.mean(delta_active_time), 2)))
    print("Average percentage error in the remake data: {} % ".format(round(np.mean(percent_difference), 2)))

    plt.figure(1)
    plt.subplot(211)
    plt.plot(nearest_brandhorst_perigee, brandhorst_active_time, label='Actual')
    plt.plot(perigees, active_times_flip / 3600.0, label='Recreation')
    plt.xlabel('Angle of Perigee [deg]')
    plt.ylabel('Total Access Time [hrs]')
    plt.title('Brandhorst Figure 3 Recreation')
    plt.legend()
    plt.subplot(212)
    plt.plot(nearest_brandhorst_perigee, percent_difference, label='Recreation - Actual')
    plt.xlabel('Angle of Perigee [deg]')
    plt.ylabel('Error in Remake Data [%]')
    plt.title('Difference')
    plt.legend()
    plt.show()


main()
