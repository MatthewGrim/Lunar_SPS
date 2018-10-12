""""
01/06/2018
Author: Darian van Paridon

This script takes the processed data for the Lunar SPS for south pole target orbit parametric scan, and applies
constraints relevant to the SPS and target systems to remove data points which do no the meet the requirements. Then
a weighted objective function is applied to determine the highest performing orbit within the feasible design space.
The resulting system performance metrics for the highest performing design are then printed.

The user can select from a small variety of parametrized transmitters and receivers to analyze the design space for
different combinations.

The user can also activate or deactivate constraints, and set the thresholds for the constraints to explore the design
space. Constraints and their activity are defined in SPS_Constrained_DesignFunctions.py

"""

from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignFunctions import *
from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignOptimizer import optimize_link_efficiency
from numpy import unravel_index
from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *


def generate_design_space(study_name, rover_selection, transmitter_selection, constraints, active_constraints, num_sps,
                          **kwargs):
    # --- INITIALIZATION ---
    study = study_initialization(study_name, **kwargs)
    transmitter = trans_metrics(transmitter_selection)
    rover = rover_metrics(rover_selection)

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    parent_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(parent_folder)
    stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

    # --- OPTIMIZE TRANSMITTER APERTURE SIZE ---
    # Based on current constraints, determine transmitter aperture size which provides highest possible link efficiency
    # within constrained design space
    optimum = optimize_link_efficiency(num_sps, transmitter_selection, rover_selection, constraints,
                                       active_constraints, study_name, study['duration'])
    transmitter['radius'] = optimum.x

    # --- READ IN DATA FILES ---
    data_set = read_in_processed_data_reports(stk_data_path, study_name, num_sps)

    # --- ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS and STATION KEEPING ---
    # Orbital perturbations on argument of perigee [0], eccentricity [1], and inclination [2]
    perturbations = calculate_orbital_perturbations(study['semi-maj-axis'], study['eccentricity'], study['inclination'], study['arg_perigee'])
    # Calculate skew in argument of perigee in degrees per year
    data_set['arg_perigee_drift'] = [abs(i * (365.0 * 24.0 * 3600.0) * 180.0 / np.pi) for i in perturbations[0]]
    data_set['eccentricity_dt'] = perturbations[1]
    data_set['inclination_dt'] = perturbations[2]

    # --- ESTIMATE DELTA V REQUIRED TO MAINTAIN ORBIT ---
    mu_moon = 6.674e-11 * 7.347673e22
    # delta v required to adjust argument of perigee to negate approximate drift rate
    data_set['delta_v_to_maintain'] = [2 * i * np.sqrt(mu_moon / (j * (1 - i ** 2))) * np.sin((k / 2) * (np.pi / 180.0)) for i, j, k in zip(study['eccentricity'], study['semi-maj-axis'], data_set['arg_perigee_drift'])]

    # --- ENFORCE CONSTRAINTS ---
    # Remove infeasible designs which do not have any active events
    for i in data_set['total_active_time']:
        if math.isnan(i):
            idx = data_set['total_active_time'].index(i)
            for j in data_set:
                data_set[j][idx] = np.nan

    # Convert data to appropriate units for applying constraints
    data_set['max_blackout_duration'] = [i / 3600.0 for i in data_set['max_blackout_duration']]
    data_set['min_active_duration'] = [i / 3600.0 for i in data_set['min_active_duration']]
    data_set['total_active_time'] = [100.0 * i / study['duration'] for i in data_set['total_active_time']]

    # Remove data points for which blackout durations exceed the limit
    if active_constraints['max_blackout'] == 1:
        data_set = enforce_constraints(data_set, 'max_blackout_duration', constraints, 'max_blackout', 'max')
    else:
        pass
    # Remove data points for which the overall blackout time is not sufficiently reduced
    if active_constraints['min_active_time'] == 1:
        data_set = enforce_constraints(data_set, 'total_active_time', constraints, 'min_active_time', 'min')
    else:
        pass

    # --- CALCULATE LINK EFFICIENCY AND POWER/ENERGY DELIVERED, APPLY POINTING CONSTRAINT ---
    data_set = calculate_link_efficiency_and_power_delivered_for_single_rover(rover, data_set, transmitter, constraints, active_constraints, study["duration"])

    # Remove data points for which not enough power is delivered on average
    if active_constraints['min_power'] == 1:
        data_set = enforce_constraints(data_set, 'min_power_received', constraints, 'min_power', 'min')
    else:
        pass

    # --- PARSING/REORGANIZING DATA ---
    # Reorganize the data lists into 2D arrays
    sorted_data_set, perigee_altitudes, unique_perigees, apogee_altitudes, unique_apogees = sort_data_lists(data_set, study['orbits'], study_name, **kwargs)

    # --- SELECT SOLUTION WITH HIGHEST LINK EFFICIENCY ---
    # Find best orbit according to weighted objective function
    best_orbit_idx = unravel_index(np.nanargmax(sorted_data_set['mean_link_efficiency']), sorted_data_set['mean_link_efficiency'].shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    eccentricity = (best_apogee - best_perigee)/(best_apogee + best_perigee)
    semi_maj_axis = best_perigee / (1 - eccentricity)
    orbit_period = 2 * np.pi * np.sqrt((semi_maj_axis * 1000.0) ** 3 / mu_moon)

    # --- OPTIMIZE TRANSMITTER POWER ---
    sorted_data_set, surf_flux, transmitter = optimize_transmitter_power(transmitter, rover, sorted_data_set, best_orbit_idx, constraints, active_constraints)
    transmitter['radius'] = optimum.x

    # # --- OPTIMIZE SOLAR ARRAY (GENERATOR) SIZE ---
    # # Power generator parameters - Stretched lens array SquareRigger platform
    solar_array_eff = 0.3
    # solar_array_spec_pwr = 300.0
    # # Metric which measures force applied by electric propulsion per unit power
    # electric_propulsion_specific_thrust = 1.0 / 20e3
    # # Add desired delta-v margin if constraint is active
    # if active_constraints['min_delta_v_margin'] == 1:
    #     delta_v_requirement = sorted_data_set['delta_v_to_maintain'][best_orbit_idx] + (constraints['min_delta_v_margin'] * 1e3)
    # else:
    #     delta_v_requirement = sorted_data_set['delta_v_to_maintain'][best_orbit_idx]
    # # Determine minimum allowable transmitter power based on ability to maintain orbit with electric propulsion
    # min_allowable_generator_pwr_delta_v = (transmitter['mass'] * delta_v_requirement) / (electric_propulsion_specific_thrust * sorted_data_set['total_station_keeping'][best_orbit_idx] - (delta_v_requirement / solar_array_spec_pwr))
    # # Minimum allowable power required from generator to operate transmitter
    # min_allowable_generator_pwr_transmitter = transmitter['power'] / transmitter['efficiency']
    # # Select largest of two minimum generator sizes
    # generator_pwr = np.max([min_allowable_generator_pwr_delta_v, min_allowable_generator_pwr_transmitter])
    # # Mass of generator
    # generator_mass = generator_pwr / solar_array_spec_pwr
    #
    # # --- DETERMINE DELTA V MARGIN FOR STATION KEEPING ---
    # # Available delta v integrated across all events, assuming 1N thrust for 20kW power
    # # Estimate delta v available based on electric propulsion system, integrated across all station keeping events
    # data_set['delta_v_available'] = [(electric_propulsion_specific_thrust * generator_pwr * i) / (generator_mass + transmitter['mass']) for i in data_set['total_station_keeping']]
    # sorted_data_set['delta_v_available'] = sort_incremented_resolution_data(study['orbits'], data_set['delta_v_available'])
    #
    # # Determine delta v margin between available and required delta v for orbit maintenance
    # data_set['delta_v_margin'] = [(i - j) / 1000.0 for i, j in zip(data_set['delta_v_available'], data_set['delta_v_to_maintain'])]
    # sorted_data_set['delta_v_margin'] = sort_incremented_resolution_data(study['orbits'], data_set['delta_v_margin'])
    #
    # if active_constraints['min_delta_v_margin'] == 1:
    #     # Remove design points for which the margin is negative
    #     for i in range(len(sorted_data_set['delta_v_margin'][0])):
    #         for j in range(len(sorted_data_set['delta_v_margin'][1])):
    #             if sorted_data_set['delta_v_margin'][i][j] < 0.0:
    #                 for k in sorted_data_set:
    #                     sorted_data_set[k][i][j] = np.nan
    # else:
    #     pass

    # --- EVALUATE STEADY STATE TEMPERATURE ---
    # Steady state temperature
    solar_irradiance = 1367.0
    emissivity = 0.8
    steady_state_temp = (solar_irradiance * (1 - solar_array_eff * transmitter['efficiency']) / (2 * emissivity * 5.67e-8)) ** 0.25

    # --- EVALUATE FLUX AND HEAT LOAD AT RECEIVER ---
    surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * (sorted_data_set['mean_range'][best_orbit_idx] * 1000.0) / (np.pi * transmitter['radius'] ** 2)) ** 2)
    target_flux = transmitter['power'] / (np.pi * surf_beam_radius ** 2)
    target_heat_load = target_flux * (1 - rover['rec_efficiency']) * np.pi * rover['rec_radius'] ** 2

    # --- DISPLAY RESULTS ---
    # Print out performance results for optimal orbit
    r_moon = 1737.0
    print('-----------------------------------------------------------------------------------------------------------')
    print('Number of SPS: {}'.format(num_sps))
    print('Optimal orbit altitudes --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Orbital period --> {} minutes'.format(round(orbit_period / 60.0, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Minimum allowable transmitter power --> {} kW'.format(round(transmitter['power'] / 1000.0, 2)))
    print('Mean flux at receiver --> {} AM0'.format(round(surf_flux['mean'] / solar_irradiance, 2)))
    print('Maximum flux at receiver --> {} AM0'.format(round(surf_flux['max'] / solar_irradiance, 2)))
    print('Minimum flux at receiver --> {} AM0'.format(round(surf_flux['min'] / solar_irradiance, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Mean link efficiency --> {} %'.format(round(sorted_data_set['mean_link_efficiency'][best_orbit_idx] * 100.0, 5)))
    print('Mean power delivered --> {} W'.format(round(sorted_data_set['mean_power_received'][best_orbit_idx], 2)))
    print('Mean heat load on receiver --> {} W'.format(round(target_heat_load, 2)))
    print('Total energy transferred --> {} MJ per year'.format(round(sorted_data_set['total_energy'][best_orbit_idx] / 2e6, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Transmitter aperture radius: {} cm'.format(round(transmitter['radius'] * 100.0, 2)))
    print('Steady state temperature --> {} Celsius'.format(round(steady_state_temp - 273.15, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Estimated argument of perigee drift rate --> {} deg/yr'.format(round(sorted_data_set['arg_perigee_drift'][best_orbit_idx], 2)))
    if best_apogee - best_perigee == 0.0:
        print('NOTE: Circular orbit, correction of argument of perigee drift practically not necessary')
    print('-----------------------------------------------------------------------------------------------------------')
    print('Total active time (blackout reduction) --> {} %'.format(round(sorted_data_set['total_active_time'][best_orbit_idx], 2)))
    print('Total blackout time --> {} %'.format(round(100.0 * sorted_data_set['total_blackout_time'][best_orbit_idx] / study['duration'], 2)))
    print('Max active period duration --> {} hours'.format(round(sorted_data_set['max_active_time'][best_orbit_idx] / 3600.0, 2)))
    print('Max blackout period duration --> {} hours'.format(round(sorted_data_set['max_blackout_duration'][best_orbit_idx], 2)))
    print('-----------------------------------------------------------------------------------------------------------')

    # Plot constrained design variables
    plt.figure(1, figsize=(15, 8))
    plt.subplot(221)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['total_active_time'], 500)
    plt.title('Total Active Time [%]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.subplot(222)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['max_blackout_duration'], 500)
    plt.title('Max Blackout Time [hrs]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.subplot(223)
    plt.contourf(apogee_altitudes, perigee_altitudes, 1e-3 * rover['operation_pwr'] / (rover['rec_efficiency'] * sorted_data_set['min_link_efficiency']), 500)
    plt.title('Laser Power [kW]')
    plt.xlabel('Apogee Altitude [km]')
    plt.colorbar()
    plt.subplot(224)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['delta_v_to_maintain'], 500)
    plt.title('$\Delta V$ Maintenence Requirement [$kms^{-1}$]')
    plt.xlabel('Apogee Altitude [km]')
    plt.colorbar()

    # Plot mean link efficiency
    r_moon = 1737.0
    plt.figure(2, figsize=(15, 8))
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['mean_link_efficiency'] * 100.0, 500)
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.scatter(best_apogee - r_moon, best_perigee - r_moon, marker='x')
    plt.show()

