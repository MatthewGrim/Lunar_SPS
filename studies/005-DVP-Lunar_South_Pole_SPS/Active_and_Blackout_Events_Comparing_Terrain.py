""""
21/05/2018
Author: Darian van Paridon

This script is used to plot a bar graph of the SPS active times and target blackout periods for a given SPS configuration.
It is currently set up to compare the eclipse events at a South pole target with and without a lunar terrain map being
included in the STK simulations.

"""

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def main():
    # Import data and set the start and end times of the simulation
    # Initialize simulation times
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Data for solar power satellite access to south pole
    raw_sps_lighting = "SPS(4.0k)_Lighting.csv"
    raw_access_sps = "SPS(4.0k)_Access.csv"

    # Data for simulation with terrain information
    raw_target_lighting = 'Target_Lighting_89.45-222.69.csv'
    raw_target_lighting_terrain = 'Target_Lighting_Terrain_89.45-222.69.csv'

    # Data for simulation with terrain information and target elevated
    raw_target_lighting_terrain_elevated = 'Target_Lighting_Terrain_Elevated_89.45-222.69.csv'
    raw_target_lighting_terrain_v2 = 'Target_Lighting_Terrain_89.45-222.69_HannaTerrain.csv'

    # Process SPS data
    sps_access = parse_csv_to_array(raw_access_sps, start)
    sps_sunlit = parse_csv_to_array(raw_sps_lighting, start)

    # Process target data without terrain
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Process target data with terrain
    target_lighting_terrain = parse_csv_to_array(raw_target_lighting_terrain, start)
    target_eclipse_terrain = invert_events_list(target_lighting_terrain, total_duration)

    # Process target data with terrain and elevation
    target_lighting_terrain_elevated = parse_csv_to_array(raw_target_lighting_terrain_elevated, start)
    target_eclipse_terrain_elevated = invert_events_list(target_lighting_terrain_elevated, total_duration)

    # Calculates the total active time for SPS, and blackout events at target
    sps_active = determine_SPS_active_time(sps_sunlit, target_eclipse, sps_access)
    target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)

    plt.figure(1)

    plt.subplot(311)
    plt.title('Simulation without Terrain')
    plt.bar([i / 86400.0 for i in target_eclipse[0]], [i / 86400.0 for i in target_eclipse[2]])
    plt.ylabel('Duration [days]')
    plt.xlim([0, 718])
    plt.ylim([0, 180])

    plt.subplot(312)
    plt.title('Simulation with Terrain')
    plt.bar([i / 86400.0 for i in target_eclipse_terrain[0]], [i / 86400.0 for i in target_eclipse_terrain[2]])
    plt.ylabel('Duration [days]')
    plt.xlim([0, 718])
    plt.ylim([0, 180])

    plt.subplot(313)
    plt.title('Simulation with Terrain and Elevated Target')
    plt.bar([i / 86400.0 for i in target_eclipse_terrain_elevated[0]], [i / 86400.0 for i in target_eclipse_terrain_elevated[2]])
    plt.ylabel('Duration [days]')
    plt.xlim([0, 718])
    plt.ylim([0, 180])
    plt.xlabel('Days')

    plt.show()


main()
