""""
25/06/2018
Author: Rohan Ramasamy

This script serves for simplified user interaction with the SPS constrained design tool.
"""

import numpy as np

from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignTool import generate_design_space
from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignFunctions import rover_metrics


def main():
    # --- INITIALIZATION ---
    # Select scenario/configuration
    # study_name = 'NorthPole_IncrementedRes_Generic'
    study_name = 'Equatorial_IncrementedRes_Generic'

    # Select transmitter
    transmitter_selection = '100kW'

    # Select receiver
    # Fleet size must be integer, rover separation distance must be float (with decimal)
    # Fleet size comes first, then separation
    rover_selection = 'sorato'

    # Defines the number of satellites in the constellation being studied
    num_sps = 1

    # --- DEFINE CONSTRAINTS ---
    # Initialize dictionaries
    constraints = {}
    active_constraints = {}

    # Retrieve rover metrics for constraint formulation
    rover = rover_metrics(rover_selection)

    # Minimum pointing error of SPS system in radians
    constraints['point_error'] = 1e-6
    # Minimum reduction in overall blackout time in percent
    constraints['min_active_time'] = 10.0
    # Minimum power requirement at target in Watts
    constraints['min_power'] = rover['operation_pwr']
    # Maximum time rover can survive without recharging in hours
    constraints['max_blackout'] = rover['battery_capacity'] / rover['hibernation_pwr']
    # Minimum allowable delta v margin related to station keeping in km/s
    constraints['min_delta_v_margin'] = 0.0

    # Specify which constraints are active
    # 1 = active, anything else = inactive
    active_constraints['point_error'] = 1
    active_constraints['min_active_time'] = 1
    active_constraints['min_power'] = 1
    active_constraints['max_blackout'] = 1
    active_constraints['min_delta_v_margin'] = 0

    active_constraints['transmitter_pwr_optimization'] = 1
    ####################################################################################################################

    kwargs = {
        "resolutions": np.array((50.0, 100.0, 100.0, 250.0)),
        "thresholds": np.array((1000.0, 1500.0, 2500.0)),
        "min_perigee": 800.0,
        "max_perigee": 5000.0,
        "max_apogee": 5000.0
    }
    generate_design_space(study_name, rover_selection, transmitter_selection, constraints, active_constraints, num_sps, **kwargs)

main()
