""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and blackout periods for a given configuration

"""

from DVP_general_SPS_functions import *

def main():
    # Import data and set the start and end times of the simulation
    raw_sunlight_SPS2 = "SPS2-Lighting(135).csv"
    raw_eclipse_SPS2 = "SPS2-Eclipse(135).csv"
    raw_access_SPS2 = "SPS2-Access(135).csv"
    raw_sunlight_SPS1 = "SPS1-Lighting(0)-Edited.csv"
    raw_eclipse_SPS1 = "SPS1-Eclipse(0)-Edited.csv"
    raw_access_SPS1 = "SPS1-Access(0).csv"
    raw_sunlight_target = 'Target1-Lighting-Edited.csv'
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Solar Power Satellite 1
    LOS_access1 = parse_csv_to_array(raw_access_SPS1, start)
    sunlight_SPS1 = parse_csv_to_array(raw_sunlight_SPS1, start)
    eclipse_SPS1 = parse_csv_to_array(raw_eclipse_SPS1, start)

    # Solar Power Satellite 2
    LOS_access2 = parse_csv_to_array(raw_access_SPS2, start)
    sunlight_SPS2 = parse_csv_to_array(raw_sunlight_SPS2, start)
    eclipse_SPS2 = parse_csv_to_array(raw_eclipse_SPS2, start)

    # Lunar Target
    sunlight_target = parse_csv_to_array(raw_sunlight_target, start)
    eclipse_target = invert_events_list(sunlight_target, total_duration)

    # Calculates the total active time for SPS, based on target access
    # and eclipses, as well as satellite illumination times
    print("\n")
    print("ACCESS AVAILABILITY for SPS1")
    sps1_active = determine_SPS_active_time(sunlight_SPS1, eclipse_target, LOS_access1)
    sps1_blackout = determine_blackout_data(sps1_active, eclipse_target, total_duration)

    print("\n")
    print("ACCESS AVAILABILITY for SPS2")
    sps2_active = determine_SPS_active_time(sunlight_SPS2, eclipse_target, LOS_access2)
    sps2_blackout = determine_blackout_data(sps2_active, eclipse_target, total_duration)

    sps1_days = [i / 86400.0 for i in sps1_blackout[0]]
    sps1_dur = [i / 3600.0 for i in sps1_blackout[2]]
    sps2_days = [i / 86400.0 for i in sps2_blackout[0]]
    sps2_dur = [i / 3600.0 for i in sps2_blackout[2]]

    plt.figure(1)
    plt.subplot(211)
    plt.bar(sps1_days, sps1_dur)
    plt.title('Blackout Durations - Brandhorst Configuration')
    plt.ylabel('Duration [h]')
    plt.subplot(212)
    plt.bar(sps2_days, sps2_dur)
    plt.xlabel('Days')
    plt.ylabel('Duration [h]')
    plt.show()


main()