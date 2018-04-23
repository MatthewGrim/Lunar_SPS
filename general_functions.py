"""
Author: Rohan Ramasamy
Date: 14/03/2018

This file contains general functions for the analysis carried out in this project
"""

import datetime


def convert_string_to_datetime(time_components):
    year = int(time_components[0])
    month = int(time_components[1])
    day = int(time_components[2])
    hour = int(time_components[3])
    minute = int(time_components[4])
    seconds = int(float(time_components[5]))

    t = datetime.datetime(year, month, day, hour, minute, seconds)

    return t


def get_illumination_event_times_from_file(solar_illumination_fname, invert=False):
    """
    Generates two lists containing the start and end illumination points respectively
    of a target in simulations

    solar_illumination_fname: file name to of csv to be read with start and end illumination times
                            the first and second elements of the file must be the start and end times
    invert: boolean to change the behaviour of code if the umbra times are passed instead of
                illumination times. This should change results as illumination times include
                periods in penumbra. This can also be used to get the eclipse times from an illumination
                file
    """
    event_times = dict()
    event_times["Start"] = list()
    event_times["End"] = list()

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
        start_time = convert_string_to_datetime(components[0].split(":"))
        end_time = convert_string_to_datetime(components[1].split(":"))

        if invert:
            if i == 1:
                event_times["Start"].append(end_time)
            elif i == num_lines - 2:
                event_times["End"].append(start_time)
            else:
                event_times["Start"].append(end_time)
                event_times["End"].append(start_time)
        else:
            event_times["Start"].append(start_time)
            event_times["End"].append(end_time)

    return event_times


def get_access_times(access_fname):
    """
    Generate two lists of the start and end points of access events between two objects

    access_fname: file name of a csv containing computed access times from STK. The second, and
                  third entries of the csv should be the start and end times respectively
    """
    access_times = dict()
    access_times["Start"] = list()
    access_times["End"] = list()

    file = open(access_fname, "r")
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
        access_start_time = convert_string_to_datetime(start)
        end = components[2].split(":")
        access_end_time = convert_string_to_datetime(end)

        access_times["Start"].append(access_start_time)
        access_times["End"].append(access_end_time)

    return access_times


def combine_events(events_1, events_2):
    """
    Takes two event dicts containing start and end events and compines them into a single list of events.

    :param events_1: dict of first events
    :param events_2: dict of second events
    :return:
    """
    assert len(events_1["Start"]) == len(events_1["End"])
    assert len(events_2["Start"]) == len(events_2["End"])
    new_event_list = dict()
    new_event_list["Start"] = list()
    new_event_list["End"] = list()

    idx_1 = 0
    idx_2 = 0
    while idx_1 < len(events_1["Start"]) and idx_2 < len(events_2["Start"]):
        start_1 = events_1["Start"]
        end_1 = events_1["End"]
        start_2 = events_2["Start"]
        end_2 = events_2["End"]

        if idx_1 >= len(events_1["Start"]) and idx_2 < len(events_2["Start"]):
            new_event_list["Start"] = start_2
            new_event_list["End"] = end_2

            idx_2 += 1
            continue
        if idx_2 >= len(events_2["Start"]):
            new_event_list["Start"] = start_1
            new_event_list["End"] = end_1

            idx_1 += 1
            continue

        if start_1 < end_1 < start_2 < end_2:
            # First event happens before second event with no overlap
            new_event_list["Start"] = start_1
            new_event_list["End"] = end_1

            idx_1 += 1
        elif start_2 < end_2 < start_1 < end_1:
            # Second event happens before first event with no overlap
            new_event_list["Start"] = start_2
            new_event_list["End"] = end_2

            idx_2 += 1
        elif start_1 < start_2 < end_2 < end_1:
            # First event contains second event
            new_event_list["Start"] = start_1
            new_event_list["End"] = end_1

            idx_1 += 1
            idx_2 += 1
        elif start_2 < start_1 < end_1 < end_2:
            # Second event contains first event
            new_event_list["Start"] = start_2
            new_event_list["End"] = end_2

            idx_1 += 1
            idx_2 += 1
        elif start_1 < start_2 < end_1 < end_2:
            # Events overlap starting at first, and ending in second
            new_event_list["Start"] = start_1
            new_event_list["End"] = end_2

            idx_1 += 1
            idx_2 += 1
        elif start_2 < start_1 < end_2 < end_1:
            # Events overlap starting at second, and ending in second
            new_event_list["Start"] = start_2
            new_event_list["End"] = end_1

            idx_1 += 1
            idx_2 += 1
        else:
            raise RuntimeError("Shouldn't be possible to get here!")

    return new_event_list


def get_overlap_between_events(first_event, second_event):
    """
    Takes two event dicts, containing start and end times of events, and computes the overlap between
    them

    :param first_event: dict of first events
    :param second_event: dict of second events
    :return:
    """
    assert len(first_event["Start"]) == len(first_event["End"])
    assert len(second_event["Start"]) == len(second_event["End"])
    new_event_times = dict()
    new_event_times["Start"] = list()
    new_event_times["End"] = list()

    first_idx = 0
    second_idx = 0
    while first_idx < len(first_event["Start"]):
        # Get start and end times
        first_start = first_event["Start"][first_idx]
        first_end = first_event["End"][first_idx]

        # If there are no more events in second list
        if second_idx >= len(second_event["Start"]):
            # Go straight to next event in first list
            new_event_times["Start"].append(first_start)
            new_event_times["End"].append(first_end)
            first_idx += 1

            continue
        else:
            second_start = second_event["Start"][second_idx]
            second_end = second_event["End"][second_idx]

        # Find Overlap times
        if second_end < first_start:
            # No overlap - move to next event pair with the chance of overlap
            second_idx += 1
        elif second_start > first_end:
            # No overlap guaranteed - move to next event pair with the chance of overlap
            first_idx += 1
        elif second_start <= first_start and second_end >= first_end:
            # Complete overlap - move to next event without adding to new events list
            new_event_times["Start"].append(first_start)
            new_event_times["End"].append(first_end)

            first_idx += 1
        elif second_start <= first_start and second_end < first_end:
            # Partial overlap from start
            new_event_times["Start"].append(first_start)
            new_event_times["End"].append(second_end)

            second_idx += 1
        elif second_start > first_start and second_end >= first_end:
            # Partial overlap from end
            new_event_times["Start"].append(second_start)
            new_event_times["End"].append(first_end)

            first_idx += 1
        elif second_start > first_start and second_end < first_end:
            # Second event is within the first event
            new_event_times["Start"].append(second_start)
            new_event_times["End"].append(second_end)

            second_idx += 1
        else:
            raise NotImplementedError("Shouldn't be possible to get here!")

    return new_event_times


def check_event_order_consistency(events):
    """
    Checks that the event order is equential. This is necessary for results to be considered valid.
    """
    assert len(events["Start"]) == len(events["End"])

    idx = 0
    last_end = None
    while idx < len(events["Start"]):
        if idx == 0:
            last_end = events["End"][idx]
            idx += 1
            continue

        start = events["Start"][idx]
        end = events["End"][idx]

        if start <= last_end:
            raise RuntimeError("Inconsistent event times!")
        if end <= start:
            raise RuntimeError("Inconsistent event times!")

        last_end = end

        idx += 1
