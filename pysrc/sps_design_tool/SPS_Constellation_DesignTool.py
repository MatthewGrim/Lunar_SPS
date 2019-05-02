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
from Lunar_SPS.pysrc.utils.physical_constants import PhysicalConstants


def generate_design_space(study_name, rover_selection, transmitter_selection, constraints, active_constraints, num_sps,
                          use_storage=False, include_misalignment=False,
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
    # optimum = optimize_link_efficiency(num_sps, transmitter_selection, rover_selection, constraints,
    #                                    active_constraints, study_name, study['duration'], use_storage)
    # transmitter['radius'] = optimum.x

    # --- READ IN DATA FILES ---
    data_set = read_in_processed_data_reports(stk_data_path, study_name, num_sps, use_storage)

    # --- OPTIMIZE TRANSMITTER APERTURE SIZE ---
    # Get optimum transmitter size for orbit
    transmitter = optimize_transmitter_radius(data_set, constraints, active_constraints, transmitter, rover)

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
    data_set = calculate_link_efficiency_and_power_delivered_for_single_rover(rover, data_set, transmitter, constraints, 
                                                                              active_constraints, study["duration"], 
                                                                              include_misalignment)

    # Remove data points for which not enough power is delivered on average
    if active_constraints['min_power'] == 1:
        data_set = enforce_constraints(data_set, 'min_power_received', constraints, 'min_power', 'min')
    else:
        pass

    # Remove data points below a certain altitude if specified
    r_moon = PhysicalConstants.r_moon
    if active_constraints['min_altitude'] == 1:
        for idx, orbit in enumerate(study['orbits']):
            altitude = orbit - r_moon
            if altitude[0] == 1300.0 and altitude[1] == 1300.0:
                print(altitude[0] < constraints['min_altitude'])
            if altitude[0] < constraints['min_altitude']:
                for j in data_set:
                    data_set[j][idx] = np.nan

    # --- PARSING/REORGANIZING DATA ---
    # Reorganize the data lists into 2D arrays
    sorted_data_set, transmitter, perigee_altitudes, unique_perigees, apogee_altitudes, unique_apogees = sort_data_lists(data_set, transmitter, study['orbits'], study_name, **kwargs)

    # --- SELECT SOLUTION WITH HIGHEST LINK EFFICIENCY ---
    # Find best orbit according to weighted objective function
    best_orbit_idx = unravel_index(np.nanargmax(sorted_data_set['mean_link_efficiency']), sorted_data_set['mean_link_efficiency'].shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    eccentricity = (best_apogee - best_perigee)/(best_apogee + best_perigee)
    semi_maj_axis = best_perigee / (1 - eccentricity)
    
    mu_moon = PhysicalConstants.mu_moon
    orbit_period = 2 * np.pi * np.sqrt((semi_maj_axis * 1000.0) ** 3 / mu_moon)

    # --- OPTIMIZE TRANSMITTER POWER ---
    # sorted_data_set, surf_flux, transmitter = optimize_transmitter_power(transmitter, rover, sorted_data_set, best_orbit_idx, constraints, active_constraints)
    # transmitter['radius'] = optimum.x
    sorted_data_set, surf_flux, transmitter = get_transmitter_power_and_surf_flux(transmitter, rover, sorted_data_set, constraints, active_constraints)

    # --- EVALUATE STEADY STATE TEMPERATURE ---
    # Steady state temperature
    solar_irradiance = PhysicalConstants.solar_irradiance
    emissivity = 0.8
    solar_array_eff = 0.3
    steady_state_temp = (solar_irradiance * (1 - solar_array_eff * transmitter['efficiency']) / (2 * emissivity * 5.67e-8)) ** 0.25

    # --- EVALUATE FLUX AND HEAT LOAD AT RECEIVER ---
    surf_beam_radius = transmitter['radius'][best_orbit_idx] * np.sqrt(1 + (transmitter['wavelength'] * (sorted_data_set['mean_range'][best_orbit_idx] * 1000.0) / (np.pi * transmitter['radius'][best_orbit_idx] ** 2)) ** 2)
    target_flux = transmitter['power'][best_orbit_idx] / (np.pi * surf_beam_radius ** 2)
    target_heat_load = target_flux * (1 - rover['rec_efficiency']) * np.pi * rover['rec_radius'] ** 2

    # ESTIMATE SPS BATTERY MASS
    # Get full transmitter power - including power beyond 1 / e2 drop off
    transmitter['power'] /= 1 - np.exp(-2)
    sps_battery_capacity = transmitter['power'][best_orbit_idx] * sorted_data_set['max_stored_power_time'][best_orbit_idx] / (3600.0 * transmitter['efficiency'])
    lipo_specific_power = 140.0
    sps_battery_mass = sps_battery_capacity / lipo_specific_power
    number_of_cycles = 24 / orbit_period * 3600.0 * 365 * 10
    charge_time = orbit_period / 2

    # --- DISPLAY RESULTS ---
    # Print out performance results for optimal orbit
    num_lines_in_header = 80

    header = ' SATELLITE ORBIT '
    num_lines = (num_lines_in_header - len(header)) // 2
    print('-' * num_lines + header + '-' * num_lines)
    print('Number of SPS: {}'.format(num_sps))
    print('Optimal orbit altitudes --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Orbital period --> {} minutes'.format(round(orbit_period / 60.0, 2)))
    print('Total active time (blackout reduction) --> {} %'.format(round(sorted_data_set['total_active_time'][best_orbit_idx], 2)))
    print('Total blackout time --> {} %'.format(round(100.0 * sorted_data_set['total_blackout_time'][best_orbit_idx] / study['duration'], 2)))
    print('Max active period duration --> {} hours'.format(round(sorted_data_set['max_active_time'][best_orbit_idx] / 3600.0, 2)))
    print('Max blackout period duration --> {} hours'.format(round(sorted_data_set['max_blackout_duration'][best_orbit_idx], 2)))
    print('Min range to target --> {} km'.format(round(sorted_data_set['min_range'][best_orbit_idx], 2)))
    print('Max range to target --> {} km'.format(round(sorted_data_set['max_range'][best_orbit_idx], 2)))

    header = ' LASER TRANSMITTER '
    num_lines = (num_lines_in_header - len(header)) // 2
    print('-' * num_lines + header + '-' * num_lines)
    print('Minimum allowable transmitter power --> {} kW'.format(round(transmitter['power'][best_orbit_idx] / 1000.0, 2)))
    print('Transmitter aperture radius: {} cm'.format(round(transmitter['radius'][best_orbit_idx] * 100.0, 2)))
    print('Min link efficiency --> {} %'.format(round(sorted_data_set['min_link_efficiency'][best_orbit_idx] * 100.0, 5)))
    print('Min power delivered --> {} W'.format(round(sorted_data_set['min_power_received'][best_orbit_idx], 2)))
    print('Mean link efficiency --> {} %'.format(round(sorted_data_set['mean_link_efficiency'][best_orbit_idx] * 100.0, 5)))
    print('Mean power delivered --> {} W'.format(round(sorted_data_set['mean_power_received'][best_orbit_idx], 2)))
    print('Steady state temperature --> {} Celsius'.format(round(steady_state_temp - 273.15, 2)))
    
    header = ' BATTERY CHARACTERISTICS '
    num_lines = (num_lines_in_header - len(header)) // 2
    print('-' * num_lines + header + '-' * num_lines)
    print('Battery capacity --> {} Whr'.format(round(sps_battery_capacity, 2)))
    print('Battery mass --> {} kg'.format(round(sps_battery_mass, 2)))
    print('Battery cycles --> {}'.format(round(number_of_cycles, 2)))
    print('Battery charge time --> {} hr'.format(round(charge_time / 3600, 2)))
    

    header = ' RECEIVER CHARACTERISTICS '
    num_lines = (num_lines_in_header - len(header)) // 2
    print('-' * num_lines + header + '-' * num_lines)
    print('Receiver Area --> {} $m^2$'.format(np.pi * rover['rec_radius'] ** 2))
    print('Mean flux at receiver --> {} AM0'.format(round(surf_flux['mean'][best_orbit_idx] / solar_irradiance, 2)))
    print('Maximum flux at receiver --> {} AM0'.format(round(surf_flux['max'][best_orbit_idx] / solar_irradiance, 2)))
    print('Minimum flux at receiver --> {} AM0'.format(round(surf_flux['min'][best_orbit_idx] / solar_irradiance, 2)))
    print('Mean heat load on receiver --> {} W'.format(round(target_heat_load, 2)))
    print('Total energy transferred --> {} MJ per year'.format(round(sorted_data_set['total_energy'][best_orbit_idx] / 2e6, 2)))
    
    best_orbit = [best_apogee - r_moon, best_perigee - r_moon]
    return apogee_altitudes, perigee_altitudes, sorted_data_set, best_orbit, transmitter

