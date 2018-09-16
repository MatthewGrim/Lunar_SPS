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
              '{} minutes of total time unaccounted for.'.format(round((duration - check_sum) / 60.0, 2)))


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
            raise RuntimeError("Inconsistent event times!")
        if end < start:
            raise RuntimeError("Inconsistent event times!")

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
    import os
    # Check if any access periods occur
    if os.stat(file_name).st_size == 1:
        # If not, write nan
        duration_array = np.nan
        start_time_sec_from_simstart = np.nan
        end_time_sec_from_simstart = np.nan
    else:
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

            if start == end:
                raise RuntimeError("No event should have zero time")
            if end_time_sec_from_simstart[-1] == start_time_sec_from_simstart[-1]:
                # Remove events that are shorter than a second long
                end_time_sec_from_simstart.pop()
                start_time_sec_from_simstart.pop()
                continue

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
    parsed_data = np.array((2, size))

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

    # This function imports range data statistics reports from STK and
    # returns the data in a tuple.

    # First entry to is minimum range for each access event
    # Second entry is maximum range for each access event
    # Third entry is mean range for each access event

    import re
    import os

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

    # Check if any access events exist
    if os.path.getsize("{}/{}_Semi_Parsed.txt".format(stk_data_path, file_name)) == 0:
        # if not, pass
        pass
    # If there are access events, parse out mean, min and max range for each
    else:
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
        last_file.close()

    # The last statistic is a section statistic. This is the average value over the whole simulation
    min_range.pop()
    max_range.pop()
    mean_range.pop()

    # The number of each range stat should be the same
    assert len(min_range) == len(max_range)
    assert len(max_range) == len(mean_range)
    range_stats = [min_range, max_range, mean_range]

    os.remove("{}/{}_Semi_Parsed.txt".format(stk_data_path, file_name))

    return range_stats


def get_event_overlaps(access_times, event_times):
    # This function finds overlap between events and access times between the SPS and the target.
    # "Conditional events" should be times when the SPS would be active, if it can access the target
    # Therefore take satellite illumination and target eclipses as event_times.
    # Output is an array containing start, end and duration of overlap events.
    overlap_start = []
    overlap_end = []
    overlap_duration = []

    # Check to see if a list of event is empty
    if not hasattr(access_times[0], "__len__") or not hasattr(event_times[0], "__len__"):
        pass
    # Else find overlap of events
    else:
        access_size = len(access_times[0])
        events_size = len(event_times[0])

        # Test all of the access periods for overlap with conditional events
        for i in range(access_size):
            # Get start and end of access periods
            access_start = access_times[0][i]
            access_end = access_times[1][i]
            assert access_start < access_end

            # Check each conditional event for overlap
            for j in range(events_size):
                event_start = event_times[0][j]
                event_end = event_times[1][j]
                assert event_start <= event_end

                # Determine overlapping events
                # Partial overlap, event period triggers overlap, access period ends it
                if access_start <= event_start and access_end <= event_end and access_end > event_start:
                    overlap_start.append(event_start)
                    overlap_end.append(access_end)
                    overlap_duration.append(access_end - event_start)
                # Partial overlap, access period triggers overlap, event period ends it
                elif event_start <= access_start and event_end <= access_end and event_end > access_start:
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
                elif access_end <= event_start:
                    pass
                # No overlap, access period begins after event
                elif event_end <= access_start:
                    pass
                else:
                    raise RuntimeError("Should not be possible to get here!")

                if len(overlap_duration) > 0 and overlap_duration[-1] == 0.0:
                    print("WARNING = This should NEVER happen!")
                    overlap_start.pop()
                    overlap_end.pop()
                    overlap_duration.pop()

    overlap_events = [overlap_start, overlap_end, overlap_duration]

    return overlap_events


def get_event_overlaps_fast(access_times, event_times):
    """
    This function finds overlap between events and access times between the SPS and the target.
    "Conditional events" should be times when the SPS would be active, if it can access the target
    Therefore take satellite illumination and target eclipses as event_times.
    Output is an array containing start, end and duration of overlap events.

    WARNING: This function only works for data that is ordered in time.

    :param access_times:
    :param event_times:
    :return:
    """
    assert len(access_times[0]) == len(access_times[1])
    assert len(event_times[0]) == len(event_times[1])
    # Check events are ordered in time
    check_event_order_consistency(access_times)
    check_event_order_consistency(event_times)

    overlap_start = []
    overlap_end = []
    overlap_duration = []
    overlap_original_idx = []

    original_exists = len(access_times) == 4
    if not original_exists:
        assert len(access_times) == 3
    access_idx = 0
    event_idx = 0
    while access_idx < len(access_times[0]):
        access_start = access_times[0][access_idx]
        access_end = access_times[1][access_idx]
        if access_start >= access_end:
            raise RuntimeError("It should not be possible to get here!")

        if event_idx >= len(event_times[0]):
            # There are no more events
            access_idx += 1

            continue
        else:
            event_start = event_times[0][event_idx]
            event_end = event_times[1][event_idx]
            if event_start >= event_end:
                raise RuntimeError("It should not be possible to get here: {} {}".format(event_start, event_end))

        # Look for no overlap cases - move to next event pair with the chance of overlap
        if event_end < access_start:
            event_idx += 1
            continue
        elif event_start > access_end:
            access_idx += 1
            continue

        # Find event overlap
        start = None
        end = None
        if event_start <= access_start and event_end >= access_end:
            # Complete overlap
            start = access_start
            end = access_end

            if original_exists:
                overlap_original_idx.append(access_times[3][access_idx])
            else:
                overlap_original_idx.append(access_idx)
            access_idx += 1
        elif event_start <= access_start and event_end < access_end:
            # Partial overlap from start
            start = access_start
            end = event_end

            if original_exists:
                overlap_original_idx.append(access_times[3][access_idx])
            else:
                overlap_original_idx.append(access_idx)
            event_idx += 1
        elif event_start > access_start and event_end >= access_end:
            # Partial overlap from end
            start = event_start
            end = access_end

            if original_exists:
                overlap_original_idx.append(access_times[3][access_idx])
            else:
                overlap_original_idx.append(access_idx)
            access_idx += 1
        elif event_start > access_start and event_end < access_end:
            # Second event is within the first event
            start = event_start
            end = event_end

            if original_exists:
                overlap_original_idx.append(access_times[3][access_idx])
            else:
                overlap_original_idx.append(access_idx)
            event_idx += 1
        else:
            raise NotImplementedError("Should not be possible to get here!")

        assert start is not None and end is not None
        # Don't add events with zero overall time
        if end == start:
            continue
        # Throw error if there's an inconsistent event
        if end < start:
            raise RuntimeError("Should not get here!")

        overlap_start.append(start)
        overlap_end.append(end)
        overlap_duration.append(end - start)

    overlap_events = [overlap_start, overlap_end, overlap_duration, overlap_original_idx]

    return overlap_events


def invert_events_list(active_times, duration):

    nc_start = []
    nc_end = []
    nc_duration = []

    if len(active_times[0]) == 0 and np.sum(active_times) == 0.0:
        no_coverage = [nc_start, nc_end, nc_duration]
        pass
    else:
        final_idx = len(active_times[1]) - 1
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

        dummy_check = get_event_overlaps_fast(active_times, no_coverage)
        if np.sum(dummy_check[2]) > 60:
            print('Original and inverted events overlapping!')

    return no_coverage


def combine_events(events_1, events_2):
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

        # If one event list is finished, add element from other list
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

        # If there is no overlap, add the first event
        if start_1 < end_1 < start_2 < end_2:
            # First event happens before second event with no overlap
            new_start.append(start_1)
            new_end.append(end_1)
            new_dur.append(end_1 - start_1)

            idx_1 += 1
            continue
        if start_2 < end_2 < start_1 < end_1:
            # Second event happens before first event with no overlap
            new_start.append(start_2)
            new_end.append(end_2)
            new_dur.append(end_2 - start_2)

            idx_2 += 1
            continue
        else:
            # Combine overlapped events recursively
            if start_1 <= start_2 and end_2 <= end_1:
                # First event contains second event
                new_start.append(start_1)
                new_end.append(end_1)
                new_dur.append(end_1 - start_1)

                idx_1 += 1
                idx_2 += 1
            elif start_2 <= start_1 and end_1 <= end_2:
                # Second event contains first event
                new_start.append(start_2)
                new_end.append(end_2)
                new_dur.append(end_2 - start_2)

                idx_1 += 1
                idx_2 += 1
            elif start_1 <= start_2 and end_1 <= end_2:
                # Events overlap starting at first, and ending in second
                new_start.append(start_1)
                new_end.append(end_2)
                new_dur.append(end_2 - start_1)

                idx_1 += 1
                idx_2 += 1
            elif start_2 <= start_1 and end_2 <= end_1:
                # Events overlap starting at second, and ending in second
                new_start.append(start_2)
                new_end.append(end_1)
                new_dur.append(end_1 - start_2)

                idx_1 += 1
                idx_2 += 1
            else:
                raise RuntimeError("Shouldn't be possible to get here!")

            # Recursively loop through lists to get new events
            while idx_1 < len(events_1[0]) and idx_2 < len(events_2[0]):
                if events_1[0][idx_1] < events_2[0][idx_2]:
                    first_idx = idx_1
                    first_list = events_1
                else:
                    first_idx = idx_2
                    first_list = events_2

                first_start_time = first_list[0][first_idx]

                if first_start_time < new_end[-1]:
                    new_end[-1] = first_list[1][first_idx]
                    new_dur[-1] = new_end[-1] - new_start[-1]

                    if events_1[0][idx_1] < events_2[0][idx_2]:
                        idx_1 += 1
                    else:
                        idx_2 += 1
                else:
                    break

    new_event_list = [new_start, new_end, new_dur]

    if np.sum(new_event_list[2]) > np.sum(events_1[2]) + np.sum(events_2[2]):
        print('More combined events than sum of individual events list!')
        raise RuntimeError

    check_event_order_consistency(new_event_list)

    return new_event_list


def determine_SPS_active_time(sunlight_sps, eclipse_target, access_times):
    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.

    # Get events which are intersection of SPS access periods, and SPS sunlit periods
    sps_available = get_event_overlaps_fast(access_times, sunlight_sps)

    # Get events which are intersection of previously mentioned events with SPS sunlit periods
    sps_active = get_event_overlaps_fast(sps_available, eclipse_target)

    print("Total SPS active time: {} hrs".format(round(np.sum(sps_active[2]) / 3600.0, 2)))

    return sps_active


def determine_SPS_storedpower_time(sps_eclipse, eclipse_target, sps_access):

    # The total active time for solar power satellite is assumed to be the total available access time, minus
    # the times when the satellite is in eclipse, and minus the times that the target is illuminated by the sun.
    sps_available = get_event_overlaps_fast(sps_access, eclipse_target)

    # Get events which are intersection of SPS access periods, and SPS eclipse periods
    sps_stored_power = get_event_overlaps_fast(sps_available, sps_eclipse)

    return sps_stored_power


def determine_battery_chargeup_events(sps_sunlit, sps_access, duration):

    # Get times when SPS cannot access target - thus it will not be beaming power
    sps_no_access = invert_events_list(sps_access, duration)

    # Get times during this when the SPS is sunlit - thus it can charge a battery
    charge_events = get_event_overlaps_fast(sps_sunlit, sps_no_access)

    return charge_events


def determine_blackout_data(active_times, eclipse_target,  duration):

    # This function determines black-out periods for the lunar target (ie. times when there is
    # no active SPS coverage, and the target is in eclipse).

    # Calculates the times when the SPS cannot beam to the target
    sps_no_coverage = invert_events_list(active_times, duration)

    # Get events where no coverage overlaps with target eclipses
    blackout_events = get_event_overlaps_fast(sps_no_coverage, eclipse_target)

    print('Total target blackout time: {} hrs'.format(round(np.sum(blackout_events[2]) / 3600.0, 2)))

    return blackout_events


def determine_field_of_view(altitude):

    # This function calculates the field of view of an SPS for a given altitude. FOV is calculated
    # as maximum change in latitude/longitude from directly below the SPS which is in view.

    R_moon = 1737.0

    theta_max = np.cos(R_moon / (R_moon + altitude))

    return theta_max


def determine_rover_battery_storage(sps_active, eclipse_times, rover_battery_capacity, rover_operation_power, rover_hibernation_power):
    """
    This function is used to carry out a power balance on the energy provided by the SPS compared to that consumed during
    the lunar night. The rover is assumed to have a full battery at the beginning of the first blackout.

    sps_active: List of SPS active events. 0: Start time, 1: End time, 2: Duration
    eclipse_times: List of target eclipses. 0: Start time, 1: End time, 2: Duration
    rover_battery_capacity: in Joules
    rover_operation_power: in Watts
    rover_hibernation_power: in Watts
    """
    assert rover_operation_power > rover_hibernation_power

    access_idx = 0
    eclipse_idx = 0
    battery_energy = []
    prev_energy_battery = rover_battery_capacity
    battery_energy.append(prev_energy_battery)
    times = [0.0]
    # Loop through all access and blackout events in temporal order
    energy_delivered = 0.0
    energy_used = 0.0
    while access_idx < len(sps_active[0]):
        while eclipse_idx < len(eclipse_times[0]):
            # Defensive checks - making sure durations are not zero and event start and finish are not equal
            if access_idx < len(sps_active[0]):
                assert not sps_active[2][access_idx] <= 0.0
                assert not sps_active[0][access_idx] == sps_active[1][access_idx], sps_active[2][access_idx]
            if eclipse_idx < len(eclipse_times[0]):
                assert not eclipse_times[2][eclipse_idx] <= 0.0
                assert not eclipse_times[0][eclipse_idx] == eclipse_times[1][eclipse_idx], eclipse_times[2][eclipse_idx]

            if access_idx < len(sps_active[0]):
                if eclipse_times[0][eclipse_idx] < sps_active[0][access_idx]:
                    # Check that current event does not run into next event
                    if not sps_active[0][access_idx] >= eclipse_times[1][eclipse_idx]:
                        string = "{}, {}, {}, {}".format(eclipse_times[0][eclipse_idx], eclipse_times[1][eclipse_idx], sps_active[0][access_idx], sps_active[1][access_idx])
                        raise RuntimeError(string)

                    energy_diff = eclipse_times[2][eclipse_idx] * rover_hibernation_power
                    prev_energy_battery -= energy_diff
                    energy_used += energy_diff
                    times.append(eclipse_times[1][eclipse_idx])
                    assert times[-1] > times[-2]
                    eclipse_idx += 1
                else:
                    # Check that current event does not run into next event
                    if not eclipse_times[0][eclipse_idx] >= sps_active[1][access_idx]:
                        string = "{}, {}, {}, {}".format(eclipse_times[0][eclipse_idx], eclipse_times[1][eclipse_idx], sps_active[0][access_idx], sps_active[1][access_idx])
                        raise RuntimeError(string)

                    energy_diff = sps_active[2][access_idx] * (rover_operation_power - rover_hibernation_power)
                    prev_energy_battery += energy_diff
                    energy_delivered += energy_diff
                    times.append(sps_active[1][access_idx])
                    assert times[-1] > times[-2]
                    access_idx += 1
            else:
                energy_diff = eclipse_times[2][eclipse_idx] * rover_hibernation_power
                prev_energy_battery -= energy_diff
                energy_used += energy_diff
                times.append(eclipse_times[1][eclipse_idx])
                assert times[-1] > times[-2]
                eclipse_idx += 1

            # If rover is going through a lunar day, reset the battery to full capacity
            if times[-1] - times[-2] > 5 * 24 * 3600.0:
                prev_energy_battery = rover_battery_capacity

            prev_energy_battery = min(prev_energy_battery, rover_battery_capacity)
            assert prev_energy_battery > 0.0
            # prev_energy_battery = max(prev_energy_battery, 0.0)
            battery_energy.append(prev_energy_battery)

        if access_idx < len(sps_active[0]):
            energy_diff = sps_active[2][access_idx] * (rover_operation_power - rover_hibernation_power)
            prev_energy_battery += energy_diff
            energy_delivered += energy_diff
            times.append(sps_active[1][access_idx])
            assert times[-1] > times[-2]
            access_idx += 1

            # If rover is going through a lunar day, reset the battery to full capacity
            if times[-1] - times[-2] > 5 * 24 * 3600.0:
                prev_energy_battery = rover_battery_capacity

            prev_energy_battery = min(prev_energy_battery, rover_battery_capacity)
            assert prev_energy_battery > 0.0
            # prev_energy_battery = max(prev_energy_battery, 0.0)
            battery_energy.append(prev_energy_battery)

    assert energy_delivered == np.sum(sps_active[2]) * (rover_operation_power - rover_hibernation_power)
    assert energy_used == np.sum(eclipse_times[2]) * rover_hibernation_power

    return times, battery_energy

