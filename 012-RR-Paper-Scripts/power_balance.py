"""
Author: Rohan Ramasamy
Date: 17/08/2018

This script is used to assess the energy balance of a solar power satellite link
"""

import os

from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def get_energy_balance():
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Name of study
    study_name = 'Equatorial_IncrementedRes'

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # File path
    stk_data_path = '{}\STK Data\{}'.format(main_directory, study_name)

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    assert os.path.isdir(stk_data_path), stk_data_path
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Rover parameters
    Whr_to_J = 3600.0
    rover_names = ["Sorato", "AMALIA"]
    rover_battery_capacity = [38.0 * Whr_to_J, 100.0 * Whr_to_J]
    rover_operation_power = [17.0, 93.0]
    rover_hibernation_power = [4.5, 7.0]
    # Import access and lighting for SPS
    perigees = [2500.0 + 1737.0, 1800.0 + 1737.0]
    apogees = [2500.0 + 1737.0, 1800.0 + 1737.0]

    fig, ax = plt.subplots(2, sharex=False)
    for i, rover_name in enumerate(rover_names):
        sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, perigees[i], apogees[i]), start)
        sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, perigees[i], apogees[i]), start)

        sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
        check_event_order_consistency(sps_active)
        target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
        check_event_order_consistency(target_blackout)

        times, battery_energy = determine_rover_battery_storage(sps_active, target_blackout, rover_battery_capacity[i], rover_operation_power[i], rover_hibernation_power[i])
        times = np.asarray(times)
        times /= (3600.0 * 24.0)
        battery_energy = np.asarray(battery_energy)
        battery_energy /= Whr_to_J

        ax[i].plot(times, battery_energy)
        ax[i].set_ylabel("{}\n Energy in battery (Whr)".format(rover_name))
        # ax[i].set_xlim((19, 35))
        ax[i].set_ylim((0.95 * np.min(battery_energy), 1.05 * np.max(battery_energy)))
    ax[1].set_xlabel("Time (days)")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    get_energy_balance()

