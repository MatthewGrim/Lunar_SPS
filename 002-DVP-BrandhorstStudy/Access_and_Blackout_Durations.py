""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and blackout periods in order to replicate the results
of the Brandhorst study: " A solar electric propulsion mission for lunar power beaming "

"""

from Lunar_SPS.DVP_general_SPS_functions import *


def main():
    # Import data and set the start and end times of the simulation
    raw_sps2_lighting = "SPS2-Lighting(135).csv"
    raw_sps2_access = "SPS2-Access(135).csv"
    raw_sps1_lighting = "SPS1-Lighting(0)-Edited.csv"
    raw_sps1_access = "SPS1-Access(0).csv"
    raw_target_lighting = 'Target1-Lighting-Edited.csv'
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Parse data into start and end time of events, in seconds with respect to the start of the simulation
    # and total duration of the event. Events refer to illumination of the SPS or target by the Sun,
    # and line-of-sight access between the target and the SPS

    # Solar Power Satellite 1
    sps1_access = parse_csv_to_array(raw_sps1_access, start)
    sps1_lighting = parse_csv_to_array(raw_sps1_lighting, start)

    # Solar Power Satellite 2
    sps2_access = parse_csv_to_array(raw_sps2_access, start)
    sps2_lighting = parse_csv_to_array(raw_sps2_lighting, start)

    # Lunar Target
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Calculates the total active time for SPS, based on target access
    # and eclipses, as well as satellite illumination times
    print("\n")
    print("ACCESS AVAILABILITY for SPS1")
    # Calculate blackout periods and active periods for SPS1 solo
    sps1_active = determine_SPS_active_time(sps1_lighting, target_eclipse, sps1_access)
    sps1_blackout = determine_blackout_data(sps1_active, target_eclipse, total_duration)

    print("\n")
    print("ACCESS AVAILABILITY for SPS2")
    # Calculate blackout periods and active periods for SPS2 solo
    sps2_active = determine_SPS_active_time(sps2_lighting, target_eclipse, sps2_access)
    sps2_blackout = determine_blackout_data(sps2_active, target_eclipse, total_duration)

    print("\n")
    print('ACCESS AVAILABILITY for PAIR')
    # Calculate blackout periods when SPS are both actively servicing target
    sps1_inactive = invert_events_list(sps1_active, total_duration)
    sps2_inactive = invert_events_list(sps2_active, total_duration)
    # Times when neither satellite can service target
    sps_pair_inactive = get_event_overlaps(sps1_inactive, sps2_inactive)
    # Need to invert to call function determine_blackout_data, which takes in active times not inactive times
    sps_pair_active = invert_events_list(sps_pair_inactive, total_duration)
    sps_pair_blackout = determine_blackout_data(sps_pair_active, target_eclipse, total_duration)

    sps1_days = [i / 86400.0 for i in sps1_blackout[0]]
    sps1_dur = [i / 3600.0 for i in sps1_blackout[2]]
    sps2_days = [i / 86400.0 for i in sps2_blackout[0]]
    sps2_dur = [i / 3600.0 for i in sps2_blackout[2]]
    sps_days = [i / 86400.0 for i in sps_pair_blackout[0]]
    sps_dur = [i / 3600.0 for i in sps_pair_blackout[2]]

    plt.figure(1)
    plt.subplot(311)
    plt.bar(sps1_days, sps1_dur)
    plt.title('Blackout Durations - Brandhorst Configuration')
    plt.ylabel('SPS 1 Solo [h]')
    plt.xlim([0, 718])
    plt.subplot(312)
    plt.bar(sps2_days, sps2_dur)
    plt.ylabel('SPS 2 Solo [h]')
    plt.xlim([0, 718])
    plt.subplot(313)
    plt.bar(sps_days, sps_dur)
    plt.xlabel('Days')
    plt.ylabel('SPS Pair [h]')
    plt.xlim([0, 718])
    plt.ylim([0, 175])
    plt.show()


main()
