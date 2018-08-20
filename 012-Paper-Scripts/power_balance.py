"""
Author: Rohan Ramasamy
Date: 17/08/2018

This script is used to assess the energy balance of a solar power satellite link
"""

from DVP_general_SPS_functions import *
import os


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

    # Import access and lighting for SPS
    perigee = 2200.0 + 1737.0
    apogee = 2200.0 + 1737.0
    sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, perigee, apogee), start)
    sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, perigee, apogee), start)

    sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
    target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)

    # Rover parameters
    Whr_to_J = 3600.0
    rover_battery_capacity = 100.0 * Whr_to_J
    rover_operation_power = 100.0
    rover_hibernation_power = 7.0
    battery_energy = determine_sps_power_balance(sps_access, target_blackout, rover_battery_capacity, rover_operation_power, rover_hibernation_power)
    battery_energy = np.asarray(battery_energy)
    if np.min(battery_energy) < 0.0:
        raise RuntimeError("Battery is depleted!")

    # Plot energy balance
    plt.figure()
    plt.plot(battery_energy / Whr_to_J)
    plt.ylabel("Energy in battery (Whr)")
    plt.title("Battery storage for AMALIA rover")
    plt.show()

if __name__ == "__main__":
    get_energy_balance()

