""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and blackout periods for a given configuration

"""

from DVP_general_SPS_functions import *


def main():
    # Import data and set the start and end times of the simulation
    raw_sunlight_sps = "SPS(4.0k)_Lighting.csv"
    raw_access_sps = "SPS(4.0k)_Access.csv"
    raw_sunlight_target = 'Target_Lighting_89.45-222.69.csv'
    raw_sunlight_target_terrain = 'Target_Lighting_Terrain_89.45-222.69.csv'
    raw_sunlight_target_terrain_elevated = 'Target_Lighting_Terrain_Elevated_89.45-222.69.csv'
    raw_sunlight_target_terrain_v2 = 'Target_Lighting_Terrain_89.45-222.69_HannaTerrain.csv'
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

    sunlight_target_terrain = parse_csv_to_array(raw_sunlight_target_terrain, start)
    eclipse_target_terrain = invert_events_list(sunlight_target_terrain, total_duration)

    sunlight_target_terrain_elevated = parse_csv_to_array(raw_sunlight_target_terrain_elevated, start)
    eclipse_target_terrain_elevated = invert_events_list(sunlight_target_terrain_elevated, total_duration)

    sunlight_target_terrain_v2 = parse_csv_to_array(raw_sunlight_target_terrain_v2, start)
    eclipse_target_terrain_v2 = invert_events_list(sunlight_target_terrain_v2, total_duration)


    # Calculates the total active time for SPS, based on target access
    # and eclipses, as well as satellite illumination times
    print("\n")
    print("ACCESS AVAILABILITY for SPS1")
    sps_active = determine_SPS_active_time(sps_sunlit, eclipse_target, sps_access)
    sps_blackout = determine_blackout_data(sps_active, eclipse_target, total_duration)

    es = [i / 86400.0 for i in eclipse_target[0]]
    ed = [i / 3600.0 for i in eclipse_target[2]]
    ets = [i / 86400.0 for i in eclipse_target_terrain[0]]
    etd = [i / 3600.0 for i in eclipse_target_terrain[2]]
    etes = [i / 86400.0 for i in eclipse_target_terrain_elevated[0]]
    eted = [i / 3600.0 for i in eclipse_target_terrain_elevated[2]]
    et2s = [i / 86400.0 for i in eclipse_target_terrain_v2[0]]
    et2d = [i / 3600.0 for i in eclipse_target_terrain_v2[2]]


    print("Max eclipse w/o terrain = {}".format(max(ed)))
    print("Max eclipse w terrain = {}".format(max(etd)))
    print("Max eclipse w terrain and elevation = {}".format(max(eted)))
    print("Max eclipse w terrain from Rohan = {}".format(max(et2d)))

    plt.figure(1)
    plt.subplot(411)
    plt.bar(es, ed)
    plt.ylabel('Duration [h]')
    plt.xlim([0, 718])
    plt.subplot(412)
    plt.bar(ets, etd)
    plt.ylabel('Duration [h]')
    plt.xlim([0, 718])
    plt.subplot(413)
    plt.bar(etes, eted)
    plt.ylabel('Duration [h]')
    plt.xlim([0, 718])
    plt.subplot(414)
    plt.bar(et2s, et2d)
    plt.ylabel('Duration [h]')
    plt.xlim([0, 718])
    plt.show()


main()
