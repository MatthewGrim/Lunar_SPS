"""
Author: Rohan Ramasamy
Date: 04/03/2018

This is a data parsing script for reeding in a csv with access times from the lunar sps study. This particular analysis
compares against the results in figure 2 of:

"A solar electric propulsion mission for lunar power beaming" - H. W. Brandhorst

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


def get_illumination_times(solar_illumination_fname, is_eclipse=False):
    """
    Generates two lists containing the start and end illumination points respectively
    """
    event_start = list()
    event_end = list()

    num_lines = sum(1 for line in open(solar_illumination_fname, "r"))
    file = open(solar_illumination_fname, "r")
    for i, line in enumerate(file):
        if i == 0:
            continue

        # Split line into table components, assuming comma delimited
        components = line.split(",")

        # For the last line, break it table elements do not exist
        if len(components) < 4:
            break

        # Get start and end time of illumination period
        start_time = get_datetime(components[0].split(":"))
        end_time = get_datetime(components[1].split(":"))

        if is_eclipse:
            if i == 1:
                event_start.append(end_time)
            elif i == num_lines - 1:
                event_end.append(start_time)
            else:
                event_start.append(end_time)
                event_end.append(start_time)
        else:
            event_start.append(start_time)
            event_end.append(end_time)

    return event_start, event_end


def update_access_and_eclipse_times(illumination_index, illumination_start, illumination_end,
                                    access_start_time, access_end_time, last_eclipse_start_time,
                                    access_durations, access_duration, eclipse_time):
    """
    Function to recursively update the access and eclipse times with the solar illumination periods

    :param illumination_index: current index of list of solar illumination periods
    :param illumination_start: list of start points of illumination periods
    :param illumination_end: list of end points of illuminations periods
    :param access_start_time: start point of current access period
    :param access_end_time: end point of current access period
    :param last_eclipse_start_time: start of last eclipse, according to only LSPS access times
    :param access_durations: list of access durations
    :param access_duration: current estimate of access duration
    :param eclipse_time: current estimate of eclipse time
    :return:
    """
    # Check illumination index is within array
    if illumination_index >= len(illumination_start):
        return access_duration, eclipse_time, illumination_index
    else:
        start_illumination = illumination_start[illumination_index]
        end_illumination = illumination_end[illumination_index]

    # Eliminate solar illuminated eclipses
    if end_illumination < access_start_time:
        # Next illumination period ends before the next eclipse
        illumination_index += 1
        return update_access_and_eclipse_times(illumination_index, illumination_start, illumination_end,
                                               access_start_time, access_end_time, last_eclipse_start_time,
                                               access_durations, access_duration, eclipse_time)
    elif start_illumination > access_end_time:
        # Next illumination period begins after this access period
        pass
    elif start_illumination < access_start_time and end_illumination > access_end_time:
        # Target is illuminated by the sun through whole access period
        eclipse_time = (start_illumination - last_eclipse_start_time).total_seconds() \
            if last_eclipse_start_time < start_illumination else None

        return None, eclipse_time, illumination_index
    elif start_illumination < access_start_time and end_illumination <= access_end_time:
        # Partial overlap from start point
        eclipse_time = (start_illumination - last_eclipse_start_time).total_seconds()
        eclipse_time = 0.0 if eclipse_time < 0.0 else eclipse_time
        access_duration = (access_end_time - end_illumination).total_seconds()

        illumination_index += 1
    elif start_illumination > access_start_time and end_illumination > access_end_time:
        access_duration = (start_illumination - access_start_time).total_seconds()
    elif start_illumination > access_start_time and end_illumination <= access_end_time:
        access_duration = (start_illumination - access_start_time).total_seconds()
        access_durations.append(access_duration)
        access_duration = (access_end_time - end_illumination).total_seconds()

        illumination_index += 1
    else:
        raise NotImplementedError()

    return access_duration, eclipse_time, illumination_index


def parse_csv_to_array(access_file_name, solar_eclipse_fname, sim_start):
    """[summary]
    
    Arguments:
        file_name {string} -- name of the .csv file to be parsed
    """
    print("Starting to parse script...")

    # Counts to make sure count the number of long eclipses - this is to compare against the statement in
    # the Brandhorst paper
    long_eclipse_time = 84.0 * 3600.0
    long_eclipse_count = 0

    illumination_start, illumination_end = get_illumination_times(solar_eclipse_fname)
    illumination_index = 0

    # Lists to store eclipse and access time data
    last_eclipse_start_time = get_datetime(sim_start)
    access_end_time = get_datetime(sim_start)
    access_durations = list()
    eclipse_durations = list()
    file = open(access_file_name, "r")
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
        access_start_time = get_datetime(start)
        eclipse_time = (access_start_time - access_end_time).total_seconds()
        end = components[2].split(":")
        access_end_time = get_datetime(end)

        # Get access time
        access_duration = float(components[3].split("\n")[0])

        # Eliminate solar illuminated eclipses
        access_duration, eclipse_time, illumination_index = update_access_and_eclipse_times(illumination_index,
                                                                                            illumination_start,
                                                                                            illumination_end,
                                                                                            access_start_time,
                                                                                            access_end_time,
                                                                                            last_eclipse_start_time,
                                                                                            access_durations,
                                                                                            access_duration,
                                                                                            eclipse_time)

        # Count long eclipses
        if eclipse_time is not None and eclipse_time > long_eclipse_time:
            long_eclipse_count += 1

        # Save power link duration and eclipse durations
        if eclipse_time is not None:
            eclipse_durations.append(eclipse_time)
        if access_duration is not None:
            access_durations.append(access_duration)

        # Update last eclipse time for next iteration
        last_eclipse_start_time = access_end_time

    # Convert results to arrays
    seconds_to_hours = 1 / 3600.0
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
    ax[1].hist(eclipse_durations, bins=200)
    plt.show()

    fig, ax = plt.subplots(2)
    ax[0].plot(access_durations)
    ax[1].plot(eclipse_durations)
    plt.show()


if __name__ == '__main__':
    use_eclipse_data = False

    # access_fname = "Brandhorst_Sim_1_Low_Res.csv"
    access_fname = "Brandhorst_Sim_1_High_Res.csv"
    illumination_fname = "Target Umbra Times.csv" if use_eclipse_data else "Target IlluminationTimes.csv"
    start = ['2008', '07', '01', '11', '0', '0.0']
    parse_csv_to_array(access_fname, illumination_fname, start)

