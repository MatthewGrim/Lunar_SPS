""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and blackout periods for a given configuration

"""

from DVP_general_SPS_functions import *


def main():
    # Import data and set the start and end times of the simulation
    raw_sunlight_sps = "SPS(4k)_Lighting.csv"
    raw_access_sps = "SPS(4k)_Access.csv"
    raw_sunlight_target = 'Target_Lighting.csv'
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Solar Power Satellite 1
    sps_access = parse_csv_to_array(raw_access_sps, start)
    sps_sunlit = parse_csv_to_array(raw_sunlight_sps, start)

    # Lunar Target
    sunlight_target = parse_csv_to_array(raw_sunlight_target, start)
    eclipse_target = invert_events_list(sunlight_target, total_duration)

    # Calculates the total active time for SPS, based on target access
    # and eclipses, as well as satellite illumination times
    print("\n")
    print("ACCESS AVAILABILITY for SPS1")
    sps_active = determine_SPS_active_time(sps_sunlit, eclipse_target, sps_access)
    sps_blackout = determine_blackout_data(sps_active, eclipse_target, total_duration)

    sps_days = [i / 86400.0 for i in sps_blackout[0]]
    sps_dur = [i / 3600.0 for i in sps_blackout[2]]

    plt.figure(1)
    plt.bar(sps_days, sps_dur)
    plt.title('Blackout Durations - Brandhorst Configuration')
    plt.ylabel('Duration [h]')
    plt.show()


main()