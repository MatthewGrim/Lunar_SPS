"""
Author: Rohan Ramasamy
Date: 29/10/2018

This script containts code to assess the variation in access of different 3 satellite constellations. This is done
in order to understand the strange behaviour in constellation results
"""

import os
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def get_access_times(study_name, data_path, target_eclipse, perigee, apogee, start, end):
    arg_perigees = [0.0, 120.0, 240.0]
    sps_active_total = None
    for j, arg in enumerate(arg_perigees):
        sim_name = 'DVP_{}_{}perigee{}apogee_{}meananom'.format(study_name, perigee, apogee, arg)

        # Import access and lighting for SPS
        sps_lighting = parse_csv_to_array('{}/{}_lighting.csv'.format(data_path, sim_name), start)
        sps_access = parse_csv_to_array('{}/{}_access.csv'.format(data_path, sim_name), start)

        # Get active time of individual satellite and combine events with total
        sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
        if j == 0:
            sps_active_total = sps_active
        else:
            sps_active_total = combine_events(sps_active_total, sps_active)

    return sps_active_total

def plot_access_times():    
    # Define start and end time
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Get pathway to SPS data directory
    study_name = "NorthPole_IncrementedRes_Generic"
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)
    stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

    # Import target illumination events
    target_lighting_raw = os.path.join(stk_data_path, 'DVP_{}_Target_Lighting.csv'.format(study_name))
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    plt.figure()
    for peri_apo_pair in [(800.0, 5000.0), (2000.0, 5000.0), (3000.0, 5000.0), (4000.0, 5000.0), (5000.0, 5000.0)]:
        print(peri_apo_pair)
        active_times = get_access_times(study_name, stk_data_path, target_eclipse, peri_apo_pair[0] + 1737, peri_apo_pair[1] + 1737, start, end)

        print("Total access time: {}".format(np.sum(active_times[2]) / 86400))
        day = [i / 86400.0 for i in active_times[0]]
        duration = [i / 86400.0 for i in active_times[2]]

        plt.plot(day, duration, label="peri-{}".format(peri_apo_pair[0]))
    plt.title("Plots of access durations for 5000km apogee polar orbits")
    plt.ylabel("Access Duration [days]")
    plt.xlabel("Simulation Day")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    plot_access_times()

