""""
25/06/2018
Author: Darian van Paridon

This script serves for simplified user interaction with the SPS constrained design tool.

"""

from SPS_Constrained_DesignTool import generate_design_space
from SPS_Constrained_DesignFunctions import rover_metrics
import numpy as np

def main():

    # INITIALIZATION
    ####################################################################################################################
    # Select scenario/configuration
    study_name = 'SouthPole_IncrementedRes_Inertial'
    # study_name = 'Brandhorst_1000.0kmRes'

    # Select transmitter
    transmitter_selection = '100kW'

    # Select receiver
    rover_selection = 'amalia'
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
    constraints['min_active_time'] = 34.0
    # Minimum power requirement at target in Watts
    constraints['min_power'] = rover['operation_pwr']
    # Maximum time rover can survive without recharging in hours
    constraints['max_blackout'] = rover['battery_capacity'] / rover['hibernation_pwr']
    # Maximum allowable skew in argument of perigee, in degrees per year
    constraints['max_arg_perigee_skew'] = 150.0

    # Specify which constraints are active
    # 1 = active, anything else = inactive
    active_constraints['point_error'] = 1
    active_constraints['min_active_time'] = 1
    active_constraints['min_power'] = 1
    active_constraints['max_blackout'] = 1
    active_constraints['max_arg_perigee_skew'] = 1
    ####################################################################################################################

    generate_design_space(study_name, rover_selection, transmitter_selection, constraints, active_constraints)


main()
