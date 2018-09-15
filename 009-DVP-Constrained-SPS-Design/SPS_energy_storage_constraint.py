"""
Author: Rohan Ramasamy
Date: 15/09/2018

This script contains a function to test that the battery in the rover is never depleted during the lunar night.
"""

import numpy as np

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *
from SPS_Constrained_DesignFunctions import rover_metrics


def energy_storage_constraint(rover_name):
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    rover = rover_metrics(rover_name)

    # Difference in subsequent apogee/perigee radii
    r_moon = 1737.0
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Name of study
    study_name = 'Equatorial_IncrementedRes'

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Get orbital data
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # READ AND PROCESS DATA REPORTS
    ####################################################################################################################
    # Initialize lists
    rover_survives = []

    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    # Cycle through available orbit configurations and calculate active/blackout durations
    file_name = "survival_constraint_{}".format(rover_name)
    Whr_to_J = 3600.0
    if os.path.exists(file_name):
        sorted_survival_grid = np.loadtxt(file_name)
    else:
        for i in range(1, len(orbit_data)):
            print('Progress: {}%'.format(round(100.0 * (i-1) / (len(orbit_data) - 2), 2)))
            print("Perigee altitude: {} km, Apogee altitude: {} km".format(orbit_data[i][0] - r_moon, orbit_data[i][1] - r_moon))

            # Import access and lighting for SPS
            sps_lighting = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)
            sps_access = parse_csv_to_array('{}/DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, orbit_data[i][0], orbit_data[i][1]), start)

            # If no access periods exist, insert nan into each list - infeasible orbit design
            if not hasattr(sps_access[0], "__len__"):
                rover_survives.append(0)

            # Otherwise insert results into list
            else:
                # Calculate active times
                sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
                # Check if there is active times, if not insert nan
                if np.sum(sps_active) == 0.0:
                    rover_survives.append(0)
                else:
                    # Blackout times
                    target_blackout = determine_blackout_data(sps_active, target_eclipse, total_duration)
                    times, battery_energy = determine_rover_battery_storage(sps_active, target_blackout,
                                                                            rover['battery_capacity'] * Whr_to_J,
                                                                            rover['operation_pwr'] - rover['hibernation_pwr'],
                                                                            rover['hibernation_pwr'])
                    if np.any(np.asarray(battery_energy) < 0.0):
                        rover_survives.append(0)
                    else:
                        rover_survives.append(1)

        sorted_survival_grid = sort_incremented_resolution_data(orbit_data, rover_survives)

        # save to file
        np.savetxt(file_name, sorted_survival_grid)

    # Extract unique perigees and apogees tested for plotting
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    r_moon = 1737.0
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    # Reduce perigee and apogee to altitudes instead of radii
    perigee_altitudes = [i - r_moon for i in unique_perigees]
    apogee_altitudes = [i - r_moon for i in unique_apogees]

    plt.figure()
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_survival_grid, 100)
    plt.title('Rover Survival')
    plt.ylabel('Perigee Altitude [km]')
    plt.xlabel('Apogee Altitude [km]')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    for rover_name in ['sorato', 'amalia']:
        energy_storage_constraint(rover_name)
