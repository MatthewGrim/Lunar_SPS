""""
06/08/2018
Author: Rohan Ramasamy

This script uses the one satellite design tool to generate results. To use this script, copy it into:
009-DVP-Constrained-SPS-Design/SPS_Design_User_file.py
"""

from SPS_Constrained_DesignTool import generate_design_space, default_design_plots
from SPS_Constrained_DesignFunctions import rover_metrics
from matplotlib import pyplot as plt


def plot_results(apogee_altitudes, perigee_altitudes, sorted_data_set, best_orbit, rover_metric):
    # Plot constrained design variables
    num_contours = 100
    fig, ax = plt.subplots(2, 4, figsize=(15, 8))

    rovers = ["Sorato", "AMALIA"]
    for i, rover in enumerate(rovers):
        im = ax[i, 0].contourf(apogee_altitudes[rover], perigee_altitudes[rover], sorted_data_set[rover]['total_active_time'], num_contours)
        ax[i, 0].set_title('Total Active Time [%]')
        ax[i, 0].set_ylabel('{}\n Perigee Altitude [km]'.format(rover))
        fig.colorbar(im, ax=ax[i, 0])

        im = ax[i, 1].contourf(apogee_altitudes[rover], perigee_altitudes[rover], sorted_data_set[rover]['max_blackout_time'], num_contours)
        ax[i, 1].set_title('Max Blackout Time [hrs]')
        fig.colorbar(im, ax=ax[i, 1])

        im = ax[i, 2].contourf(apogee_altitudes[rover], perigee_altitudes[rover], sorted_data_set[rover]['mean_power_received'], num_contours)
        ax[i, 2].set_title('Mean Power [W]')
        fig.colorbar(im, ax=ax[i, 2])

        im = ax[i, 3].contourf(apogee_altitudes[rover], perigee_altitudes[rover], sorted_data_set[rover]['mean_link_efficiency'] * 100.0, num_contours)
        ax[i, 3].set_title('Mean Link Efficiency [%]')
        ax[i, 3].scatter(best_orbit[rover][0], best_orbit[rover][1], marker='x')
        fig.colorbar(im, ax=ax[i, 3])
    for i in range(4):
        ax[1, i].set_xlabel("Apogee Altitude [km]")
    fig.tight_layout()
    plt.show()

    fig, ax = plt.subplots(2, sharex=True)
    for i, rover in enumerate(rovers):
        power_delivered_to_consumed = (rover_metric[rover]["operation_pwr"] - rover_metric[rover]["hibernation_pwr"]) / rover_metric[rover]["hibernation_pwr"]
        active_to_eclipsed_ratio = sorted_data_set[rover]["mean_active_time"] / sorted_data_set[rover]["mean_blackout_time"]
        power_balance = active_to_eclipsed_ratio * power_delivered_to_consumed
        im = ax[i].contourf(apogee_altitudes[rover], perigee_altitudes[rover], power_balance, num_contours)
        fig.colorbar(im, ax=ax[i])

        ax[i].set_ylabel("{}\n Perigee Altitude [km]".format(rover))
    ax[1].set_xlabel("Apogee Altitude [km]")
    plt.show()


def main():

    # INITIALIZATION
    ####################################################################################################################
    # Select scenario/configuration
    study_name = 'Equatorial_IncrementedRes'
    # study_name = 'Equatorial_IncrementedRes'

    # Select transmitter
    transmitter_selection = '100kW'

    # Select receiver
    # Fleet size must be integer, rover separation distance must be float (with decimal)
    # Fleet size comes first, then separation
    rover_apogees = dict()
    rover_perigees = dict()
    rover_data = dict()
    rover_orbit = dict()
    rover_metric = dict()
    rovers = ['Sorato', 'AMALIA']
    for rover_selection in rovers:
        ####################################################################################################################

        # DEFINE CONSTRAINTS
        ####################################################################################################################
        # Initialize dictionaries
        constraints = {}
        active_constraints = {}

        # Retrieve rover metrics for constraint formulation
        rover_name = 'sorato' if rover_selection is 'Sorato' else 'amalia'
        rover = rover_metrics(rover_name)
        rover_metric[rover_selection] = rover

        # Minimum pointing error of SPS system in radians
        constraints['point_error'] = 1e-6
        # Minimum reduction in overall blackout time in percent
        constraints['min_active_time'] = 49.76 * rover["hibernation_pwr"] / (rover["operation_pwr"])
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
        active_constraints['point_error'] = 1
        active_constraints['min_active_time'] = 1
        active_constraints['min_active_duration'] = 0
        active_constraints['min_power'] = 1
        active_constraints['max_blackout'] = 1
        active_constraints['min_delta_v_margin'] = 0

        active_constraints['transmitter_pwr_optimization'] = 1
        ####################################################################################################################

        apogees, perigees, sorted_data, best_orbit = generate_design_space(study_name, rover_name, transmitter_selection, constraints, active_constraints, include_tracking=False)
        rover_apogees[rover_selection] = apogees
        rover_perigees[rover_selection] = perigees
        rover_data[rover_selection] = sorted_data
        rover_orbit[rover_selection] = best_orbit

    plot_results(rover_apogees, rover_perigees, rover_data, rover_orbit, rover_metric)

main()

