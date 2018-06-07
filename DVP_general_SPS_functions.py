import numpy as np
import datetime
import matplotlib.pyplot as plt


def convert_string_to_datetime(time_components):
    year = int(time_components[0])
    month = int(time_components[1])
    day = int(time_components[2])
    hour = int(time_components[3])
    minute = int(time_components[4])
    seconds = int(float(time_components[5]))

    if seconds == 60:
        seconds = 59

    t = datetime.datetime(year, month, day, hour, minute, seconds)

    return t


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def check_sum(lighting, eclipse, duration):
    # A quick check to ensure that every time stamp from the STK simulation is accounted for
    # and belongs to an event. Threshold for "no error" is set to 30 seconds.
    # Used to check for errors when inverting an event set (ie. find eclipses from sunlit events)
    check_sum = (np.sum(lighting[2]) + np.sum(eclipse[2]))
    if abs(check_sum - duration) > 30:
        print('Inverted and original event times do not add to total duration! '
              '{} hours of total time unaccounted for.'.format(round((duration - check_sum) / 3600.0, 2)))


def check_event_order_consistency(events):
    """
    Checks that the event order is sequential. This is necessary for results to be considered valid.
    """
    assert len(events[0]) == len(events[1])

    idx = 1
    while idx < len(events[0]):
        last_end = events[1][idx - 1]
        start = events[0][idx]
        end = events[1][idx]

        if start < last_end:
            print("Inconsistent event times!")
        if end < start:
            print("Inconsistent event times!")

        idx += 1


def parse_csv_to_array(file_name, sim_start):
    # This function takes in the .csv file from STK, and converts it into a 3 X N array, where N is the number of events
    # in the file. The three data types are beginning of event, end of event (both measured in seconds from the
    # beginning of the simulation), and duration of event. The .csv files are edited in Excel to remove an additional
    # information which cause errors when running them through this script.
    the_file = open(file_name, "r")
    duration_array = list()
    start_time_sec_from_simstart = list()
    end_time_sec_from_simstart = list()

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
        start_time_sec_from_simstart.append((start_time - sim_start).total_seconds())

        end = components[1].split(":")
        end_time = convert_string_to_datetime(end)
        end_time_sec_from_simstart.append((end_time - sim_start).total_seconds())

        # Save power link duration
        access_durations = float(components[2].split("\n")[0])
        duration_array.append(access_durations)

    parsed_data = [start_time_sec_from_simstart, end_time_sec_from_simstart, duration_array]

    return parsed_data


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
        # Break once the final line is reached (which is a blank 'return')
        if components[0] == '\n':
            break

        # Work out sunlit times and set new time t
        start = components[0].split(":")
        start_time = convert_string_to_datetime(start)
        start_time_sec_from_simstart[i - 1] = (start_time - sim_start).total_seconds()

        sps_range[i - 1] = components[1]

        parsed_data = [start_time_sec_from_simstart, sps_range]

    return parsed_data


def import_range_data_statistics(file_name, stk_data_path):

    import re

    the_file = open("{}/{}.txt".format(stk_data_path, file_name), "r")
    new_file = open("{}/{}_Semi_Parsed.txt".format(stk_data_path, file_name), "w+")
    mean_range = []
    max_range = []
    min_range = []

    # Remove lines which do not contain data, and sort columns
    # by removing white space
    for i, line in enumerate(the_file):
        if line[0] == "M":
            removed_white_space = re.split('\s{2,}', line.strip())
            new_file.write("{}\n".format(removed_white_space))
    the_file.close()
    new_file.close()

    last_file = open("{}/{}_Semi_Parsed.txt".format(stk_data_path, file_name), "r")
    # Sort data into minimum maximum and mean arrays
    for j, new_line in enumerate(last_file):
        new_line_split = new_line.split(",")
        if new_line_split[0][2:-1] == "Min Range":
            min_range.append(float(new_line_split[2][2:-3]))
        elif new_line_split[0][2:-1] == "Max Range":
            max_range.append(float(new_line_split[2][2:-3]))
        elif new_line_split[0][2:-1] == "Mean Range":
            mean_range.append(float(new_line_split[1][2:-3]))
    new_file.close()

    range_stats = [min_range, max_range, mean_range]

    return range_stats


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

    # Test all of the access periods for overlap with conditional events
    for i in range(access_size):
        # Get start and end of access periods
        access_start = access_times[0][i]
        access_end = access_times[1][i]

        # Check each conditional event for overlap
        for j in range(events_size):
            event_start = event_times[0][j]
            event_end = event_times[1][j]

            # Determine overlapping events
            # Partial overlap, event period triggers overlap, access period ends it
            if access_start <= event_start and access_end <= event_end and access_end >= event_start:
                overlap_start.append(event_start)
                overlap_end.append(access_end)
                overlap_duration.append(access_end - event_start)
                # i += 1
            # Partial overlap, access period triggers overlap, event period ends it
            elif event_start <= access_start and event_end <= access_end and event_end >= access_start:
                overlap_start.append(access_start)
                overlap_end.append(event_end)
                overlap_duration.append(event_end - access_start)
                # j += 1
            # Full overlap, access period occurs within event period
            elif event_start <= access_start and access_end <= event_end:
                overlap_start.append(access_start)
                overlap_end.append(access_end)
                overlap_duration.append(access_end - access_start)
                # i += 1
            # Full overlap, event period occurs within access period
            elif access_start <= event_start and event_end <= access_end:
                overlap_start.append(event_start)
                overlap_end.append(event_end)
                overlap_duration.append(event_end - event_start)
                # j += 1
            # No overlap, access period ends before event
            elif access_end < event_start:
                # i += 1
                pass
            # No overlap, access period begins after event
            elif event_end < access_start:
                # j += 1
                pass
            else:
                raise RuntimeError("Should not be possible to get here!")

    overlap_events = [overlap_start, overlap_end, overlap_duration]

    return overlap_events


def invert_events_list(active_times, duration):

    assert len(active_times[0]) == len(active_times[1])
    final_idx = len(active_times[1]) - 1
    nc_start = list()
    nc_end = list()
    nc_duration = list()

    # For the case where no_coverage coincides with neither beginning nor end of simulation
    # Iterate through access periods, and calculate leading no_coverage period
    if active_times[0][0] != 0 and active_times[1][final_idx] != duration:
        nc_start.append(0)
        nc_end.append(active_times[0][0])
        nc_duration.append(active_times[0][0])
        for j in range(0, final_idx):
            nc_start.append(active_times[1][j])
            nc_end.append(active_times[0][j+1])
            nc_duration.append(active_times[0][j+1] - active_times[1][j])
        nc_start.append(active_times[1][final_idx])
        nc_end.append(duration)
        nc_duration.append(duration - active_times[1][final_idx])

    # For the case where no_coverage coincides with beginning but not end of simulation
    # Iterate through access periods and calculate leading no_coverage period
    elif active_times[0][0] == 0 and active_times[1][final_idx] != duration:
        for j in range(0, final_idx):
            nc_start.append(active_times[1][j])
            nc_end.append(active_times[0][j+1])
            nc_duration.append(active_times[0][j+1] - active_times[1][j])
        nc_start.append(active_times[1][final_idx])
        nc_end.append(duration)
        nc_duration.append(duration - active_times[1][final_idx])

    # For the case where no_coverage coincides with end but not beginning of simulation
    # Iterate through access periods and calculate leading no_coverage period
    elif active_times[0][0] != 0 and active_times[1][final_idx] == duration:
        nc_start.append(0)
        nc_end.append(active_times[0][0])
        nc_duration.append(active_times[0][0])
        for j in range(0, final_idx):
            nc_start.append(active_times[1][j])
            nc_end.append(active_times[0][j + 1])
            nc_duration.append(active_times[0][j + 1] - active_times[1][j])

    # For the case where no_coverage coincides with beginning and end of simulation
    # Iterate through access periods and calculate leading no_coverage period
    elif active_times[0][0] == 0 and active_times[1][final_idx] == duration:
        for j in range(0, final_idx):
            nc_start.append(active_times[1][j])
            nc_end.append(active_times[0][j + 1])
            nc_duration.append(active_times[0][j + 1] - active_times[1][j])

    no_coverage = [nc_start, nc_end, nc_duration]

    check_sum(active_times, no_coverage, duration)

    dummy_check = get_event_overlaps(active_times, no_coverage)
    if np.sum(dummy_check[2]) > 60:
        print('Original and inverted events overlapping!')

    return no_coverage


def combine_events(events_1, events_2):
    """
    Takes two event dicts containing start and end events and compines them into a single list of events.

    :param events_1: dict of first events
    :param events_2: dict of second events
    :return:
    """
    assert len(events_1[0]) == len(events_1[1])
    assert len(events_2[0]) == len(events_2[1])
    new_start = []
    new_end = []
    new_dur = []

    idx_1 = 0
    idx_2 = 0
    while idx_1 < len(events_1[0]) and idx_2 < len(events_2[0]):
        start_1 = events_1[0][idx_1]
        end_1 = events_1[1][idx_1]
        start_2 = events_2[0][idx_2]
        end_2 = events_2[1][idx_2]

        if idx_1 >= len(events_1[0]) and idx_2 < len(events_2[0]):
            new_start.append(start_2)
            new_end.append(end_2)
            new_dur.append(end_2 - start_2)

            idx_2 += 1
            continue
        if idx_2 >= len(events_2[0]):
            new_start.append(start_1)
            new_end.append(end_1)
            new_dur.append(end_1 - start_1)

            idx_1 += 1
            continue
        if start_1 < end_1 < start_2 < end_2:
            # First event happens before second event with no overlap
            new_start.append(start_1)
            new_end.append(end_1)
            new_dur.append(end_1 - start_1)

            idx_1 += 1
        elif start_2 < end_2 < start_1 < end_1:
            # Second event happens before first event with no overlap
            new_start.append(start_2)
            new_end.append(end_2)
            new_dur.append(end_2 - start_2)

            idx_2 += 1
        elif start_1 <= start_2 < end_2 <= end_1:
            # First event contains second event
            new_start.append(start_1)
            new_end.append(end_1)
            new_dur.append(end_1 - start_1)

            idx_1 += 1
            idx_2 += 1
        elif start_2 <= start_1 < end_1 <= end_2:
            # Second event contains first event
            new_start.append(start_2)
            new_end.append(end_2)
            new_dur.append(end_2 - start_2)

            idx_1 += 1
            idx_2 += 1
        elif start_1 <= start_2 < end_1 <= end_2:
            # Events overlap starting at first, and ending in second
            new_start.append(start_1)
            new_end.append(end_2)
            new_dur.append(end_2 - start_1)

            idx_1 += 1
            idx_2 += 1
        elif start_2 <= start_1 < end_2 <= end_1:
            # Events overlap starting at second, and ending in second
            new_start.append(start_2)
            new_end.append(end_1)
            new_dur.append(end_1 - start_2)

            idx_1 += 1
            idx_2 += 1
        else:
            raise RuntimeError("Shouldn't be possible to get here!")

    new_event_list = [new_start, new_end, new_dur]

    return new_event_list


def determine_SPS_active_time(sunlight_sps, eclipse_target, access_times):
    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.
    print("\n")
    # Get events which are intersection of SPS access periods, and SPS sunlit periods
    sps_available = get_event_overlaps(access_times, sunlight_sps)
    # Get events which are intersection of previously mentioned events with SPS sunlit periods
    sps_active = get_event_overlaps(sps_available, eclipse_target)
    # Calculate total active time for SPS (sum of durations)
    total_sps_time = np.sum(sps_active[2])
    print("Total time which SPS is sunlit and beaming (active): {} hrs".format(round(total_sps_time / 3600.0, 2)))
    print("Maximum active event duration: {} hrs".format(round(max(sps_active[2]) / 3600.0, 2)))
    return sps_active


def determine_SPS_storedpower_time(sps_eclipse, eclipse_target, sps_access):
    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.
    print("\n")
    # Get events which are intersection of SPS access periods, and SPS eclipse periods
    sps_available = get_event_overlaps(sps_access, eclipse_target)
    sps_uses_stored_power = get_event_overlaps(sps_available, sps_eclipse)
    total_availability = np.sum(sps_uses_stored_power[2])
    print("Total time which SPS could use stored power: {} hrs".format(round(total_availability / 3600.0, 2)))
    sps_stored_power = get_event_overlaps(sps_available, sps_eclipse)
    print("Maximum duration for which stored power would be required: {} hrs".format(round(max(sps_stored_power[2]) / 3600.0, 2)))
    return sps_stored_power


def determine_battery_chargeup_events(sps_sunlit, sps_access, duration):

    # Get times when SPS cannot access target - thus it will not be beaming power
    sps_no_access = invert_events_list(sps_access, duration)
    # Get times during this when the SPS is sunlit - thus it can charge a battery
    charge_events = get_event_overlaps(sps_sunlit, sps_no_access)
    return charge_events


def determine_blackout_data(active_times, eclipse_target,  duration):

    # This function determines black-out periods for the lunar target (ie. times when there is
    # no active SPS coverage, and the target is in eclipse).

    # Calculates the times when the SPS cannot beam to the target
    sps_no_coverage = invert_events_list(active_times, duration)

    # Get events where no coverage overlaps with target eclipses
    dark_events = get_event_overlaps(sps_no_coverage, eclipse_target)

    # Calculate number of eclipses exceeding 84 hours, and the maximum eclipse duration
    dark_durations = np.array(dark_events[2])
    long_eclipse_flag = (dark_durations / 3600.0) > 6.0
    num_long_eclipse = np.sum(long_eclipse_flag)
    max_eclipse_duration = round(max([i/3600.0 for i in dark_events[2]]), 2)
    print('Maximum black-out duration with SPS: {} hrs'.format(max_eclipse_duration))
    print('Percent per year spent in black-out with SPS: {}% '.format(round(100.0 * np.sum(dark_durations) / duration, 2)))
    print("Number of times black-out duration exceeds six hours: {}".format(num_long_eclipse))

    return dark_events


def determine_field_of_view(altitude):

    # This function calculates the field of view of an SPS for a given altitude. FOV is calculated
    # as maximum change in latitude/longitude from directly below the SPS which is in view.

    R_moon = 1737.0

    theta_max = np.cos(R_moon / (R_moon + altitude))

    return theta_max
