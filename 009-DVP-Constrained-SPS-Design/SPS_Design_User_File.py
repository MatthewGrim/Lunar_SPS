""""
25/06/2018
Author: Darian van Paridon

This script serves for simple user interaction with the SPS constrained design tool.

"""

from SPS_Constrained_DesignTool import generate_design_space
from SPS_Constrained_DesignFunctions import rover_metrics


def main():

    # SELECT TRANSMITTER
    ####################################################################################################################
    transmitter_selection = '4kW'
    ####################################################################################################################

    # SELECT RECEIVER
    ####################################################################################################################
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
    constraints['min_active_time'] = 0.0
    # Minimum power requirement at target in Watts
    constraints['min_power'] = rover['operation_pwr']
    # Maximum time rover can survive without recharging in hours
    constraints['max_blackout'] = rover['battery_capacity'] / rover['hibernation_pwr']
    # Maximum allowable skew in argument of perigee, in degrees per year
    constraints['max_arg_perigee_skew'] = 113.0

    # Specify which constraints are active
    # 1 = active, anything else = inactive
    active_constraints['point_error'] = 1
    active_constraints['min_active_time'] = 0
    active_constraints['min_power'] = 1
    active_constraints['max_blackout'] = 1
    active_constraints['max_arg_perigee_skew'] = 0
    ####################################################################################################################

    generate_design_space(rover_selection, transmitter_selection, constraints, active_constraints)


main()
