""""
Date: 14/05/2018
Author: Darian van Paridon

This file is for verifying the results of Brandhorst (shown in Figures 8 and 9 in his lunar SPS paper ___).

"""

from DVP_general_SPS_functions import *
import matplotlib.pylab as plt


def import_range_data(file_name, sim_start):
    the_file = open(file_name, "r")
    size = file_len(file_name)
    sps_range = np.zeros(size - 1)
    start_time_sec_from_simstart = np.zeros(size - 1)
    parsed_data = np.array((3, size))

    # For loop parses data into three categories: start, end and duration of event
    # This function can be used for target to satellite access times, or target/satellite illumination events
    # .csv file must be three columns, with start time, end time, duration in that order
    # The function outputs the start and end times of the event in seconds since the beginning of the simulation
    for i, line in enumerate(the_file):
        if i == 0:
            continue
        # Split line into table components, assuming comma delimited
        components = line.split(",")

        # Work out sunlit times and set new time t
        start = components[0].split(":")
        start_time = convert_string_to_datetime(start)
        start_time_sec_from_simstart[i - 1] = (start_time - sim_start).total_seconds()

        sps_range[i - 1] = components[1]

        parsed_data = [start_time_sec_from_simstart, sps_range]

    return parsed_data


start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
total_duration = (end - start).total_seconds()
raw_range_data = 'SPS-to-Target-Range.csv'
range_data = import_range_data(raw_range_data, start)

plt.figure(1)
plt.plot(range_data[0] / 86400.0, range_data[1])
plt.show()
