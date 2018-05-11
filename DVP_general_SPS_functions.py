import numpy as np
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


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def parse_csv_to_array(file_name, sim_start):
    # This function takes in the .csv file from STK, and converts it into a 3 X N array, where N is the number of events
    # in the file. The three data types are beginning of event, end of event (both measured in seconds from the
    # beginning of the simulation), and duration of event. The .csv files are edited in Excel to remove an additional
    # information which cause errors when running them through this script.
    the_file = open(file_name, "r")
    size = file_len(file_name)
    duration_array = np.zeros(size)
    start_time_sec_from_simstart = np.zeros(size)
    end_time_sec_from_simstart = np.zeros(size)
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