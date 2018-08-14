""""
25/06/2018
Author: Darian van Paridon

This script serves for simplified user interaction with the SPS constrained design tool.

"""

from SPS_Constellation_DesignTool import generate_design_space, default_design_plots
from SPS_Constellation_DesignFunctions import rover_metrics


def main():

    # INITIALIZATION
    ####################################################################################################################
    # Select scenario/configuration
    study_name = 'SouthPole_IncrementedRes_ManytoOne'
    # study_name = 'Equatorial_IncrementedRes_ManytoOne'

    # Select transmitter
    transmitter_selection = '100kW'

    # Select receiver
    # Fleet size must be integer, rover separation distance must be float (with decimal)
    # Fleet size comes first, then separation
    rover_selection = 'sorato'
    ####################################################################################################################

    # DEFINE CONSTRAINTS
    ####################################################################################################################
    # Initialize dictionaries
    constraints = {}
    active_constraints = {}

    # Retrieve rover metrics for constraint formulation
    rover = rover_metrics(rover_selection)

    # Minimum pointing error of SPS system in radians
    constraints['point_error'] = 1e-6
    # Minimum reduction in overall blackout time in percent
    constraints['min_active_time'] = 10.0
    # Minimum allowable single active event duration in hours
    constraints['min_active_duration'] = rover['battery_capacity'] / rover['operation_pwr']
    # Minimum power requirement at target in Watts
    constraints['min_power'] = rover['operation_pwr']
    # Maximum time rover can survive without recharging in hours
    constraints['max_blackout'] = rover['battery_capacity'] / rover['hibernation_pwr']
    # Minimum allowable delta v margin related to station keeping in km/s
    constraints['min_delta_v_margin'] = 0.0

    # Specify which constraints are active
    # 1 = active, anything else = inactive
    active_constraints['point_error'] = 0
    active_constraints['min_active_time'] = 0
    active_constraints['min_active_duration'] = 0
    active_constraints['min_power'] = 0
    active_constraints['max_blackout'] = 0
    active_constraints['min_delta_v_margin'] = 0

    active_constraints['transmitter_pwr_optimization'] = 0
    ####################################################################################################################

    apogees, perigees, sorted_data, best_orbit = generate_design_space(study_name, rover_selection, transmitter_selection, constraints, active_constraints)
    default_design_plots(apogees, perigees, sorted_data, best_orbit)

main()
