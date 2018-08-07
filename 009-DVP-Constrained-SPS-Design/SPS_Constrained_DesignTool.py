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

from DVP_Programmatic_Functions import *
from SPS_Constrained_DesignOptimizer import optimize_link_efficiency
from SPS_Constrained_DesignFunctions import *
from numpy import unravel_index


def generate_design_space(study_name, rover_selection, transmitter_selection, constraints, active_constraints):

    # INITIALIZATION
    ####################################################################################################################
    # Specific info for selected study
    study = study_initialization(study_name)

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    # Retrieve chosen transmitter design metrics
    transmitter = trans_metrics(transmitter_selection)

    # Retrieve chosen rover design metrics
    rover = rover_metrics(rover_selection)
    ####################################################################################################################

    # OPTIMIZE TRANSMITTER APERTURE SIZE
    ####################################################################################################################
    # Based on current constraints, determine transmitter aperture size which provides highest possible link efficiency
    # within constrained design space
    optimum = optimize_link_efficiency(transmitter_selection, rover_selection, constraints, active_constraints, study_name)
    transmitter['radius'] = optimum.x
    ####################################################################################################################

    # READ IN DATA FILES
    ####################################################################################################################
    data_set = {}
    data_set['total_active_time'] = read_data_from_file(stk_data_path, study_name, "TotalActive")
    data_set['min_active_time'] = read_data_from_file(stk_data_path, study_name, "MinActive")
    data_set['max_active_time'] = read_data_from_file(stk_data_path, study_name, "MaxActive")
    data_set['mean_active_time'] = read_data_from_file(stk_data_path, study_name, "MeanActive")
    data_set['std_active_time'] = read_data_from_file(stk_data_path, study_name, "StdActive")

    data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackout")
    data_set['max_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MaxBlackout")
    data_set['mean_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MeanBlackout")
    data_set['min_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MinBlackout")
    data_set['std_blackout_time'] = read_data_from_file(stk_data_path, study_name, "StdBlackout")

    data_set['mean_range'] = read_data_from_file(stk_data_path, study_name, "MeanRange")
    data_set['max_range'] = read_data_from_file(stk_data_path, study_name, "MeanMaxRange")
    data_set['min_range'] = read_data_from_file(stk_data_path, study_name, "MeanMinRange")

    data_set['max_stored_power_time'] = read_data_from_file(stk_data_path, study_name, 'MaxStoredPowerEvent')
    data_set['total_stored_power_time'] = read_data_from_file(stk_data_path, study_name, 'TotalStoredPowerEvent')
    data_set['total_station_keeping'] = read_data_from_file(stk_data_path, study_name, 'TotalStationKeeping')
    ####################################################################################################################

    # ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS and STATION KEEPING
    ####################################################################################################################
    # Orbital perturbations on argument of perigee [0], eccentricity [1], and inclination [2]
    perturbations = calculate_orbital_perturbations(study['semi-maj-axis'], study['eccentricity'], study['inclination'], study['arg_perigee'])
    # Calculate skew in argument of perigee in degrees per year
    data_set['arg_perigee_drift'] = [abs(i * study['duration']) for i in perturbations[0]]
    data_set['eccentricity_dt'] = perturbations[1]
    data_set['inclination_dt'] = perturbations[2]
    ####################################################################################################################

    # ESTIMATE DELTA V REQUIRED TO MAINTAIN ORBIT
    ####################################################################################################################
    mu_moon = 6.674e-11 * 7.347673e22
    # delta v required to adjust argument of perigee to negate approximate drift rate
    data_set['delta_v_to_maintain'] = [2 * i * np.sqrt(mu_moon / (j * (1 - i ** 2))) * np.sin((k / 2) * (np.pi / 180.0)) for i, j, k in zip(study['eccentricity'], study['semi-maj-axis'], data_set['arg_perigee_drift'])]
    ####################################################################################################################

    # ENFORCE CONSTRAINTS
    ####################################################################################################################
    # Remove infeasible designs which do not have any active events
    if 'Equatorial' in study_name:
        for i in data_set['total_active_time']:
            if math.isnan(i):
                idx = data_set['total_active_time'].index(i)
                for j in data_set:
                    data_set[j][idx] = np.nan

    # Convert data to appropriate units for applying constraints
    data_set['max_blackout_time'] = [i / 3600.0 for i in data_set['max_blackout_time']]
    data_set['mean_blackout_time'] = [i / 3600.0 for i in data_set['mean_blackout_time']]
    data_set['min_blackout_time'] = [i / 3600.0 for i in data_set['min_blackout_time']]
    data_set['std_blackout_time'] = [i / 3600.0 for i in data_set['std_blackout_time']]

    data_set['max_active_time'] = [i / 3600.0 for i in data_set['max_active_time']]
    data_set['mean_active_time'] = [i / 3600.0 for i in data_set['mean_active_time']]
    data_set['min_active_time'] = [i / 3600.0 for i in data_set['min_active_time']]
    data_set['std_active_time'] = [i / 3600.0 for i in data_set['std_active_time']]

    data_set['total_active_time'] = [100.0 * i / study['duration'] for i in data_set['total_active_time']]

    # Remove data points for which blackout durations exceed the limit
    if active_constraints['max_blackout'] == 1:
        data_set = enforce_constraints(data_set, 'max_blackout_time', constraints, 'max_blackout', 'max')
    else:
        pass
    # Remove data points for which the overall blackout time is not sufficiently reduced
    if active_constraints['min_active_time'] == 1:
        data_set = enforce_constraints(data_set, 'total_active_time', constraints, 'min_active_time', 'min')
    else:
        pass
    # Remove data points for which the the battery could not be charged sufficiently
    if active_constraints['min_active_duration'] == 1:
        data_set = enforce_constraints(data_set, 'min_active_time', constraints, 'min_active_duration', 'min')
    else:
        pass
    # Calculate link efficiency and power delivered, applying pointing constraint
    if "fleet" in rover_selection:
        data_set = calculate_link_efficiency_and_power_delivered_for_fleet(rover, data_set, transmitter, constraints, active_constraints)
    else:
        data_set = calculate_link_efficiency_and_power_delivered_for_single_rover(rover, data_set, transmitter, constraints, active_constraints)

    # Remove data points for which not enough power is delivered on average
    if active_constraints['min_power'] == 1:
        data_set = enforce_constraints(data_set, 'min_power_received', constraints, 'min_power', 'min')
    else:
        pass
    ####################################################################################################################

    # PARSING/REORGANIZING DATA
    ####################################################################################################################
    # Reorganize the data lists into 2D arrays
    sorted_data_set = sort_data_lists(data_set, study['orbits'], study_name)

    # Extract unique perigees and apogees tested for plotting
    unique_perigees = [study['orbits'][1][0]]
    unique_apogees = [study['orbits'][1][1]]
    r_moon = 1737.0
    for i in range(1, len(study['orbits'])):
        if study['orbits'][i][0] > max(unique_perigees):
            unique_perigees.append(study['orbits'][i][0])
        if study['orbits'][i][1] > max(unique_apogees):
            unique_apogees.append(study['orbits'][i][1])

    # Reduce perigee and apogee to altitudes instead of radii
    perigee_altitudes = [i - r_moon for i in unique_perigees]
    apogee_altitudes = [i - r_moon for i in unique_apogees]

    # Plot Event statistics
    get_event_statistics(apogee_altitudes, perigee_altitudes, sorted_data_set, rover, 'Event Statistics')
    ####################################################################################################################

    # SELECT SOLUTION WITH HIGHEST LINK EFFICIENCY
    ####################################################################################################################
    # Find best orbit according to weighted objective function
    best_orbit_idx = unravel_index(np.nanargmax(sorted_data_set['mean_link_efficiency']), sorted_data_set['mean_link_efficiency'].shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    eccentricity = (best_apogee - best_perigee)/(best_apogee + best_perigee)
    semi_maj_axis = best_perigee / (1 - eccentricity)
    orbit_period = 2 * np.pi * np.sqrt((semi_maj_axis * 1000.0) ** 3 / mu_moon)
    ####################################################################################################################

    # OPTIMIZE TRANSMITTER POWER
    ####################################################################################################################
    sorted_data_set, surf_flux, transmitter = optimize_transmitter_power(transmitter, rover, sorted_data_set, best_orbit_idx, constraints, active_constraints)
    transmitter['radius'] = optimum.x
    ####################################################################################################################

    # OPTIMIZE SOLAR ARRAY (GENERATOR) SIZE
    ####################################################################################################################
    # Power generator parameters - Stretched lens array SquareRigger platform
    solar_array_eff = 0.3
    solar_array_spec_pwr = 300.0
    # Metric which measures force applied by electric propulsion per unit power
    thrust_to_power_ratio = 1.0 / 20e3
    # Add desired delta-v margin if constraint is active
    if active_constraints['min_delta_v_margin'] == 1:
        delta_v_requirement = sorted_data_set['delta_v_to_maintain'][best_orbit_idx] + (constraints['min_delta_v_margin'] * 1e3)
    else:
        delta_v_requirement = sorted_data_set['delta_v_to_maintain'][best_orbit_idx]
    # Determine minimum allowable transmitter power based on ability to maintain orbit with electric propulsion
    min_allowable_generator_pwr_delta_v = (transmitter['mass'] * delta_v_requirement) / (thrust_to_power_ratio * sorted_data_set['total_station_keeping'][best_orbit_idx] - (delta_v_requirement / solar_array_spec_pwr))
    # Minimum allowable power required from generator to operate transmitter
    min_allowable_generator_pwr_transmitter = transmitter['power'] / transmitter['efficiency']
    # Select largest of two minimum generator sizes
    generator_pwr = np.max([min_allowable_generator_pwr_delta_v, min_allowable_generator_pwr_transmitter])
    # Mass of generator
    generator_mass = generator_pwr / solar_array_spec_pwr
    ####################################################################################################################

    # DETERMINE DELTA V MARGIN FOR STATION KEEPING
    ####################################################################################################################
    # Available delta v integrated across all events, assuming 1N thrust for 20kW power
    # Estimate delta v available based on electric propulsion system, integrated across all station keeping events
    number_impulses = 2 * study['duration'] / orbit_period
    duration_of_impulse = 5 * 60.0
    delta_v_available = (number_impulses * duration_of_impulse * thrust_to_power_ratio * generator_pwr) / (generator_mass + transmitter['mass'])

    # Determine delta v margin between available and required delta v for orbit maintenance
    data_set['delta_v_margin'] = [(delta_v_available - j) / 1000.0 for j in data_set['delta_v_to_maintain']]
    sorted_data_set['delta_v_margin'] = sort_incremented_resolution_data(study['orbits'], data_set['delta_v_margin'])

    if active_constraints['min_delta_v_margin'] == 1:
        # Remove design points for which the margin is negative
        for i in range(len(sorted_data_set['delta_v_margin'][0])):
            for j in range(len(sorted_data_set['delta_v_margin'][1])):
                if sorted_data_set['delta_v_margin'][i][j] < constraints['min_delta_v_margin']:
                    for k in sorted_data_set:
                        sorted_data_set[k][i][j] = np.nan
    else:
        pass
    ####################################################################################################################

    # EVALUATE HEAT GENERATED AND STEADY STATE TEMPERATURE
    ####################################################################################################################
    # Heat absorbed by generator
    generator_heat = (1 - solar_array_eff) * generator_pwr
    # Heat generated by transmitter
    transmitter_heat = (1 - transmitter['efficiency']) * transmitter['power']
    # Total heat in system
    total_heat = transmitter_heat + generator_heat
    # Steady state temperature
    solar_irradiance = 1367.0
    emissivity = 0.8
    steady_state_temp = (solar_irradiance * (1 - solar_array_eff * transmitter['efficiency']) / (2 * emissivity * 5.67e-8)) ** 0.25
    ####################################################################################################################

    # ESTIMATE SPS BATTERY MASS
    ####################################################################################################################
    sps_battery_capacity = transmitter['power'] * sorted_data_set['max_stored_power_time'][best_orbit_idx] / 3600.0 * transmitter['efficiency']
    # Lithium polymer battery
    lipo_specific_power = 950.0
    # Fuel cell generator
    fuel_cell_specific_pwr = 1500.0
    sps_battery_mass = sps_battery_capacity / lipo_specific_power
    sps_fuel_cell_mass = sps_battery_capacity / fuel_cell_specific_pwr
    ####################################################################################################################

    # EVALUATE FLUX AND HEAT LOAD AT RECEIVER
    ####################################################################################################################
    surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * (sorted_data_set['mean_range'][best_orbit_idx] * 1000.0) / (np.pi * transmitter['radius'] ** 2)) ** 2)
    target_flux = transmitter['power'] / (np.pi * surf_beam_radius ** 2)
    target_heat_load = target_flux * (1 - rover['rec_efficiency']) * np.pi * rover['rec_radius'] ** 2
    ####################################################################################################################

    # DISPLAY RESULTS
    ####################################################################################################################
    # Print out performance results for optimal orbit
    print('-----------------------------------------------------------------------------------------------------------')
    print('Optimal orbit altitudes --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Orbital period --> {} minutes'.format(round(orbit_period / 60.0, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Minimum allowable transmitter power --> {} kW'.format(round(transmitter['power'] / 1000.0, 2)))
    print('Mean flux at receiver --> {} AM0'.format(round(surf_flux['mean'] / solar_irradiance, 2)))
    print('Maximum flux at receiver --> {} AM0'.format(round(surf_flux['max'] / solar_irradiance, 2)))
    print('Minimum flux at receiver --> {} AM0'.format(round(surf_flux['min'] / solar_irradiance, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    if "fleet" in rover_selection:
        print('Mean link efficiency (per rover) --> {} %'.format(round(sorted_data_set['mean_link_efficiency'][best_orbit_idx] * 100.0, 5)))
        print('Mean power delivered (per rover) --> {} W'.format(round(sorted_data_set['mean_power_received'][best_orbit_idx], 2)))
    else:
        print('Mean link efficiency --> {} %'.format(round(sorted_data_set['mean_link_efficiency'][best_orbit_idx] * 100.0, 5)))
        print('Mean power delivered --> {} W'.format(round(sorted_data_set['mean_power_received'][best_orbit_idx], 2)))
    print('Mean heat load on receiver --> {} W'.format(round(target_heat_load, 2)))
    print('Total energy transferred --> {} MJ per year'.format(round(sorted_data_set['total_energy'][best_orbit_idx] / 2e6, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Transmitter aperture radius: {} cm'.format(round(transmitter['radius'] * 100.0, 2)))
    print('Combined mass of generator and transmitter --> {} kg'.format(round(generator_mass + transmitter['mass'], 2)))
    print('Combined heat load of generator and transmitter --> {} kW'.format(round(total_heat / 1000.0, 2)))
    print('Steady state temperature --> {} Celsius'.format(round(steady_state_temp - 273.15, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Estimated argument of perigee drift rate --> {} deg/yr'.format(round(sorted_data_set['arg_perigee_drift'][best_orbit_idx] * (180.0 / np.pi), 2)))
    print('Estimated delta v available to correct drift with EP --> {} km/s'.format(round(delta_v_available / 1000.0, 2)))
    print('Estimated delta v required to negate drift --> {} km/s'.format(round(sorted_data_set['delta_v_to_maintain'][best_orbit_idx] / 1000.0, 2)))
    if best_apogee - best_perigee == 0.0:
        print('NOTE: Circular orbit, correction of argument of perigee drift practically not necessary')
    print('-----------------------------------------------------------------------------------------------------------')
    print('Total active time (blackout reduction) --> {} %'.format(round(sorted_data_set['total_active_time'][best_orbit_idx], 2)))
    print('Total blackout time --> {} %'.format(round(100.0 * sorted_data_set['total_blackout_time'][best_orbit_idx] / study['duration'], 2)))
    print('Max active period duration --> {} hours'.format(round(sorted_data_set['max_active_time'][best_orbit_idx] / 3600.0, 2)))
    print('Max blackout period duration --> {} hours'.format(round(sorted_data_set['max_blackout_time'][best_orbit_idx], 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    if math.isnan(sorted_data_set['max_stored_power_time'][best_orbit_idx]):
        print('No events when SPS could use stored power for transmission link')
    else:
        print('Max event duration during which SPS requires stored power --> {} hours'.format(round(sorted_data_set['max_stored_power_time'][best_orbit_idx] / 3600.0, 2)))
        print('Approximate battery mass required to eliminate max duration event --> {} kg'.format(round(sps_battery_mass, 2)))
        print('Approximate fuel cell mass required to eliminate max duration event --> {} kg'.format(round(sps_fuel_cell_mass, 2)))
        print('Total time blackout time which could be eliminated with battery --> {} %'.format(round(100.0 * sorted_data_set['total_stored_power_time'][best_orbit_idx] / study['duration'], 2)))

    # Plot constrained design variables
    plt.figure(1)
    plt.subplot(221)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['total_active_time'], 500)
    plt.title('Total Active Time [%]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.subplot(222)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['max_blackout_time'], 500)
    plt.title('Max Blackout Time [hrs]')
    plt.colorbar()
    plt.subplot(223)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['delta_v_margin'], 500)
    plt.title('$\Delta v$ Margin [km/s]')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.subplot(224)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['mean_power_received'], 500)
    plt.title('Mean Power [W]')
    plt.xlabel('Apogee Altitude [km]')
    plt.colorbar()

    # Plot mean link efficiency
    plt.figure(2)
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['mean_link_efficiency'] * 100.0, 500)
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.scatter(best_apogee - r_moon, best_perigee - r_moon, marker='x')
    plt.show()

