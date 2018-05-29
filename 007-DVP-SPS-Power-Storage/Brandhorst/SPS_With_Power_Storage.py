""""
27/05/2018
Author: Darian van Paridon

This script takes calculates the active and blackout durations for an SPS configuration, including the
events when the SPS is eclipsed during an access period. This simulates SPS uses stored power. As a result,
the script also determines the available charging events (ie. when the SPS is sunlit, but not in range of the target)
"""

from DVP_general_SPS_functions import *


def main():
    # Set simulation start and end dates
    # Data is from Brandhorst orbital study (002-DVP-SPS_Orbit_Study)
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Import SPS illumination events
    sps_lighting_raw = 'Lighting(30k).csv'
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    sps_eclipse = invert_events_list(sps_lighting, total_duration)

    # Import access data for sps
    raw_access = "Access(30k).csv"
    sps_access = parse_csv_to_array(raw_access, start)
    no_access = invert_events_list(sps_access, total_duration)

    # Import target lighting data, and invert to find eclipse events
    raw_target_lighting = "Target1-Lighting-Edited.csv"
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    print("\n")
    print('Maximum black-out duration w/o SPS: {} hrs'.format(round(max(target_eclipse[2]) / 3600.0, 2)))
    print('Total time spent in black-out w/o SPS: {}%'.format(
        round(100.0 * np.sum(target_eclipse[2]) / total_duration, 2)))

    sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
    sps_stored_power = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)
    all_access_times = combine_events(sps_active, sps_stored_power)

    print("\n")
    print("For case where SPS only active while sunlit:")
    determine_blackout_data(sps_active, target_eclipse, total_duration)
    print("\n")
    print("For case where SPS uses stored power when eclipsed:")
    determine_blackout_data(all_access_times, target_eclipse, total_duration)

    # Calculate events when SPS could charge a battery
    # Overlap of events when target not within line of sight, and when sps is sunlit
    charge_events = determine_battery_chargeup_events(sps_lighting, sps_access, total_duration)
    print("\n")
    print("Total time available for charging batteries: {} hrs".format(round(np.sum(charge_events[2]) / 3600.0, 2)))
    print("Maximum charging event duration: {} hrs".format(round(max(charge_events[2]) / 3600.0, 2)))

    plt.subplot(311)
    plt.bar([i / 86400.0 for i in sps_active[0]], [i / 3600.0 for i in sps_active[2]])
    plt.ylabel("SPS Active Event Duration [h]")
    plt.title('Brandhorst Configuration (Single Satellite)')
    plt.subplot(312)
    plt.bar([i / 86400.0 for i in sps_stored_power[0]], [i / 3600.0 for i in sps_stored_power[2]])
    plt.ylabel("Access-While-Eclipsed Duration [h]")
    plt.subplot(313)
    plt.bar([i / 86400.0 for i in charge_events[0]], [i / 3600.0 for i in charge_events[2]])
    plt.ylabel("Charging Event Duration [h]")
    plt.xlabel('Days')
    plt.show()


main()
