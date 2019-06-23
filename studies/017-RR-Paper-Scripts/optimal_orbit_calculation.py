"""
Author: Rohan Ramasamy
Date: 06/11/2018

This file contains optimisation results for constellations of solar power satellites targeting a lunar lander at 
45 degrees North of the Equator
"""

import numpy as np
from matplotlib import pyplot as plt

from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignTool import generate_design_space
from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignFunctions import rover_metrics
from Lunar_SPS.pysrc.utils.physical_constants import PhysicalConstants
from Lunar_SPS.pysrc.utils.unit_conversions import UnitConversions


def main():
    # --- INITIALIZATION ---
    # Select scenario/configuration
    study_name = 'Equatorial_IncrementedRes_Generic'
    study_type = study_name.split('_')[0]

    # Select transmitter
    transmitter_selection = 'L-BEAM'

    # Select receiver
    rover_selection = 'rover_concept_five'

    # Defines the number of satellites in the constellation being studied
    num_sps = 1

    # --- DEFINE CONSTRAINTS ---
    # Initialize dictionaries
    constraints = {}
    active_constraints = {}

    # Retrieve rover metrics for constraint formulation
    rover = rover_metrics(rover_selection)

    # Minimum pointing error of SPS system in radians
    constraints['point_error'] = 8e-7
    # Minimum reduction in overall blackout time in percent
    constraints['min_active_time'] = 49.76 * rover["hibernation_pwr"] / (rover["operation_pwr"])
    # Minimum power requirement at target in Watts
    constraints['min_power'] = rover['operation_pwr']
    # Maximum time rover can survive without recharging in hours
    constraints['max_blackout'] = rover['battery_capacity'] / rover['hibernation_pwr']
    # Minimum allowable delta v margin related to station keeping in km/s
    constraints['min_delta_v_margin'] = 0.0
    constraints['min_altitude'] = 1250.0

    # Specify which constraints are active
    # 1 = active, anything else = inactive
    active_constraints['point_error'] = 1
    active_constraints['min_active_time'] = 1
    active_constraints['min_power'] = 1
    active_constraints['max_blackout'] = 1
    active_constraints['min_delta_v_margin'] = 0
    active_constraints['min_altitude'] = 0

    active_constraints['transmitter_pwr_optimization'] = 1
    ####################################################################################################################

    kwargs = {
        "resolutions": np.array((50.0, 100.0, 100.0, 250.0)),
        "thresholds": np.array((1000.0, 1500.0, 2500.0)),
        "min_perigee": 800.0,
        "max_perigee": 5000.0,
        "max_apogee": 5000.0
    }
    apogee_altitudes, perigee_altitudes, sorted_data_set, best_orbit, transmitter = generate_design_space(study_name, rover_selection, transmitter_selection, 
                                                                                             constraints, active_constraints, num_sps, 
                                                                                             use_storage=True, 
                                                                                             **kwargs)
    # Plot constrained design variables
    plt.figure(1, figsize=(12, 7))
    plt.subplot(221)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['total_active_time'], 500)
    plt.title('Total Active Time [%]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.subplot(222)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['max_blackout_duration'], 500)
    plt.title('Max Blackout Time [hrs]')
    plt.colorbar()
    plt.subplot(223)
    plt.contourf(apogee_altitudes, perigee_altitudes, transmitter['power'] * UnitConversions.W_to_kW, 500)
    plt.title('Laser Power [kW]')
    plt.ylabel('Perigee Altitude [km]')
    plt.xlabel('Apogee Altitude [km]')
    plt.colorbar()
    plt.subplot(224)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['mean_link_efficiency'] * 100.0, 500)
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.scatter(best_orbit[0], best_orbit[1], marker='x')
    suptitle = "{} {} results".format(study_type, rover_selection)
    plt.suptitle(suptitle)
    plt.savefig('{}_{}_results'.format(study_type, rover_selection))
    plt.show()

main()

