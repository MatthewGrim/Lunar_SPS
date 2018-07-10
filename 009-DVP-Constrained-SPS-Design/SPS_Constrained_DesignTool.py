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

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *
from SPS_Constrained_DesignOptimizer import optimize_link_efficiency
from SPS_Constrained_DesignFunctions import *


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
    data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackout")
    data_set['max_active_time'] = read_data_from_file(stk_data_path, study_name, "MaxActive")
    data_set['max_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MaxBlackout")
    data_set['mean_range'] = read_data_from_file(stk_data_path, study_name, "MeanRange")
    data_set['max_range'] = read_data_from_file(stk_data_path, study_name, "MeanMaxRange")
    data_set['min_range'] = read_data_from_file(stk_data_path, study_name, "MeanMinRange")
    data_set['mean_active_time'] = read_data_from_file(stk_data_path, study_name, "MeanActive")
    data_set['mean_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MeanBlackout")
    data_set['max_stored_power_time'] = read_data_from_file(stk_data_path, study_name, 'MaxStoredPowerEvent')
    data_set['total_stored_power_time'] = read_data_from_file(stk_data_path, study_name, 'TotalStoredPowerEvent')
    ####################################################################################################################

    # ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS
    ####################################################################################################################
    # Orbital perturbations on argument of perigee [0], eccentricity [1], and inclination [2]
    perturbations = calculate_orbital_perturbations(study['semi-maj-axis'], study['eccentricity'], study['inclination'], study['arg_perigee'])
    # Calculate skew in argument of perigee in degrees per year
    data_set['arg_perigee_skew'] = [abs(i * (365.0 * 24.0 * 3600.0) * 180.0 / np.pi) for i in perturbations[0]]
    data_set['eccentricity_dt'] = perturbations[1]
    ####################################################################################################################

    # CALCULATE LINK EFFICIENCY AND POWER/ENERGY DELIVERED, APPLY POINTING CONSTRAINT
    ####################################################################################################################
    # Initialize lists
    data_set['min_link_efficiency'] = []
    data_set['min_power_received'] = []
    data_set['mean_link_efficiency'] = []
    data_set['mean_power_received'] = []
    # Cycle through all access events
    for i in range(0, len(data_set['max_range'])):
        # If range is "np.nan", no access periods exist for that orbit, therefore treat orbit as infeasible design
        if data_set['max_range'] == np.nan:
            data_set['min_link_efficiency'].append(np.nan)
            data_set['mean_link_efficiency'].append(np.nan)
            data_set['min_power_received'].append(np.nan)
            data_set['mean_power_received'].append(np.nan)
        # Otherwise calculate efficiency and power, apply pointing constraints
        else:
            # Minimum beam radius as defined by pointing error and minimum range of access period
            min_beam_radius = rover['rec_radius'] + (constraints['point_error'] * data_set['min_range'][i] * 1000.0)

            # Actual beam radius as defined by Gaussian beam divergence
            max_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * (data_set['max_range'][i] * 1000.0) / (np.pi * transmitter['radius'] ** 2)) ** 2)
            mean_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * (data_set['mean_range'][i] * 1000.0) / (np.pi * transmitter['radius'] ** 2)) ** 2)

            # IF constraint is active, remove data points which violate pointing constraint
            if active_constraints['point_error'] == 1:
                if max_surf_beam_radius < min_beam_radius:
                    data_set['min_link_efficiency'].append(np.nan)
                    data_set['mean_link_efficiency'].append(np.nan)
                    data_set['min_power_received'].append(np.nan)
                    data_set['mean_power_received'].append(np.nan)
                    # Apply constraint all other data lists
                    for j in data_set:
                        data_set[j][i] = np.nan
                else:
                    data_set['min_link_efficiency'].append((rover['rec_radius'] / max_surf_beam_radius) ** 2)
                    data_set['mean_link_efficiency'].append((rover['rec_radius'] / mean_surf_beam_radius) ** 2)
                    data_set['min_power_received'].append(data_set['min_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])
                    data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])

            # ELSE calculate as normal
            else:
                # Calculate min link efficiency based on maximum surface beam size
                if max_surf_beam_radius <= rover['rec_radius']:
                    data_set['min_link_efficiency'].append(1.0)
                    data_set['min_power_received'].append(rover['rec_efficiency'] * transmitter['power'])
                else:
                    data_set['min_link_efficiency'].append((rover['rec_radius'] / max_surf_beam_radius) ** 2)
                    data_set['min_power_received'].append(data_set['min_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])
                # Calculate mean efficiency based on mean surface beam size
                if mean_surf_beam_radius <= rover['rec_radius']:
                    data_set['mean_link_efficiency'].append(1.0)
                    data_set['mean_power_received'].append(rover['rec_efficiency'] * transmitter['power'])
                else:
                    data_set['mean_link_efficiency'].append((rover['rec_radius'] / mean_surf_beam_radius) ** 2)
                    data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * transmitter['power'])

    # Calculate total energy delivered to receiver based on mean power
    data_set['total_energy'] = []
    data_set['total_energy'] = [i * j for i, j in zip(data_set['mean_power_received'], data_set['total_active_time'])]
    ####################################################################################################################

    # ENFORCE CONSTRAINTS
    ####################################################################################################################
    # Convert data to appropriate units for applying constraints
    data_set['max_blackout_time'] = [i / 3600.0 for i in data_set['max_blackout_time']]
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
    # Remove data points for which the estimated skew in argument of perigee is too great
    if active_constraints['max_arg_perigee_skew'] == 1:
        data_set = enforce_constraints(data_set, 'arg_perigee_skew', constraints, 'max_arg_perigee_skew', 'max')
    else:
        pass
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
    ####################################################################################################################

    # SELECT SOLUTION WITH HIGHEST LINK EFFICIENCY
    ####################################################################################################################
    # Find best orbit according to weighted objective function
    from numpy import unravel_index
    best_orbit_idx = unravel_index(np.nanargmax(sorted_data_set['mean_link_efficiency']), sorted_data_set['mean_link_efficiency'].shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    eccentricity = (best_apogee - best_perigee)/(best_apogee + best_perigee)
    semi_maj_axis = best_perigee / (1 - eccentricity)
    mu_moon = 6.674e-11 * 7.347673e22
    orbit_period = 2 * np.pi * np.sqrt((semi_maj_axis * 1000.0) ** 3 / mu_moon)
    ####################################################################################################################

    # SELECT SMALLEST TRANSMITTER WHICH ACCOMMODATES POWER
    ####################################################################################################################
    max_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['max_range'][best_orbit_idx] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)
    min_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['min_range'][best_orbit_idx] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)
    mean_surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * sorted_data_set['mean_range'][best_orbit_idx] * 1000.0 / (np.pi * transmitter['radius'] ** 2)) ** 2)

    # Determine minimum allowable transmitter power
    min_allowable_transmitter_pwr = constraints['min_power'] * max_surf_beam_radius ** 2 / (rover['rec_efficiency'] * rover['rec_radius'] ** 2)
    # Adjust mean power delivered according to the minimum allowable transmitter power
    sorted_data_set['mean_power_received'] = sorted_data_set['mean_power_received'] * (min_allowable_transmitter_pwr / transmitter['power'])
    # Switch to smallest transmitter which accommodates power
    if min_allowable_transmitter_pwr < 4000.0:
        transmitter = trans_metrics('4kW')
    elif 4000.0 < min_allowable_transmitter_pwr < 15000.0:
        transmitter = trans_metrics('15kW')
    elif min_allowable_transmitter_pwr > 15000.0:
        transmitter = trans_metrics('100kW')
    transmitter['radius'] = optimum.x

    # Assume transmitter operates at minimum allowable power
    transmitter['power'] = min_allowable_transmitter_pwr

    # Calculate flux variations
    mean_surf_flux = min_allowable_transmitter_pwr / (np.pi * mean_surf_beam_radius ** 2)
    min_surf_flux = min_allowable_transmitter_pwr / (np.pi * max_surf_beam_radius ** 2)
    max_surf_flux = min_allowable_transmitter_pwr / (np.pi * min_surf_beam_radius ** 2)
    ####################################################################################################################

    # EVALUATE REQUIRED POWER GENERATOR METRICS
    ####################################################################################################################
    # Power generator parameters - Stretched lens array SquareRigger platform
    solar_array_eff = 0.3
    solar_array_spec_pwr = 300.0
    # Power required from generator
    generator_pwr = transmitter['power'] / transmitter['efficiency']
    # Mass of generator
    generator_mass = generator_pwr / solar_array_spec_pwr
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
    lipo_specific_power = 270.0
    # Fuel cell generator
    fuel_cell_specific_pwr = 500.0
    sps_battery_mass = sps_battery_capacity / lipo_specific_power
    sps_fuel_cell_mass = sps_battery_capacity / fuel_cell_specific_pwr
    ####################################################################################################################

    # EVALUATE FLUX AND HEAT LOAD AT RECEIVER
    ####################################################################################################################
    surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * (sorted_data_set['mean_range'][best_orbit_idx] * 1000.0) / (np.pi * transmitter['radius'] ** 2)) ** 2)
    target_flux = transmitter['power'] / (np.pi * surf_beam_radius ** 2)
    target_heat_load = target_flux * (1 - rover['rec_efficiency']) * np.pi * rover['rec_radius'] ** 2
    ####################################################################################################################

    # ESTIMATE STATION KEEPING
    ####################################################################################################################
    # Find events when SPS could perform station keeping
    sps_lighting_raw = '{}\DVP_{}_{}perigee{}apogee_lighting.csv'.format(stk_data_path, study_name, best_perigee, best_apogee)
    sps_lighting = parse_csv_to_array(sps_lighting_raw, study['start'])
    sps_access_raw = '{}\DVP_{}_{}perigee{}apogee_access.csv'.format(stk_data_path, study_name, best_perigee, best_apogee)
    sps_access = parse_csv_to_array(sps_access_raw, study['start'])
    sps_station_keeping_events = determine_battery_chargeup_events(sps_lighting, sps_access, study['duration'])

    # Available delta v integrated across all events, assuming 1N thrust for 20kW power
    delta_v = 1.0 * np.sum(sps_station_keeping_events[2]) / (generator_mass + transmitter['mass'])

    # delta v required to adjust argument of perigee to negate slew rate
    delta_v_to_maintain = 2 * eccentricity * np.sqrt(mu_moon / (semi_maj_axis * (1 - eccentricity ** 2))) * np.sin((sorted_data_set['arg_perigee_skew'][best_orbit_idx] / 2) * (np.pi / 180.0))
    ####################################################################################################################

    # DISPLAY RESULTS
    ####################################################################################################################
    # Print out performance results for optimal orbit
    print('-----------------------------------------------------------------------------------------------------------')
    print('Optimal orbit altitudes --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Orbital period --> {} minutes'.format(round(orbit_period / 60.0, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Minimum allowable transmitter power --> {} kW'.format(round(transmitter['power'] / 1000.0, 2)))
    print('Mean flux at receiver --> {} AM0'.format(round(mean_surf_flux / solar_irradiance, 2)))
    print('Maximum flux at receiver --> {} AM0'.format(round(max_surf_flux / solar_irradiance, 2)))
    print('Minimum flux at receiver --> {} AM0'.format(round(min_surf_flux / solar_irradiance, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Mean link efficiency --> {} %'.format(round(sorted_data_set['mean_link_efficiency'][best_orbit_idx] * 100.0, 5)))
    print('Mean power received --> {} W'.format(round(sorted_data_set['mean_power_received'][best_orbit_idx], 2)))
    print('Mean heat load on receiver --> {} W'.format(round(target_heat_load, 2)))
    print('Total energy transferred --> {} MJ per year'.format(round(sorted_data_set['total_energy'][best_orbit_idx] / 2e6, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Transmitter aperture radius: {} cm'.format(round(transmitter['radius'] * 100.0, 2)))
    print('Combined mass of generator and transmitter --> {} kg'.format(round(generator_mass + transmitter['mass'], 2)))
    print('Combined heat load of generator and transmitter --> {} kW'.format(round(total_heat / 1000.0, 2)))
    print('Steady state temperature --> {} Celsius'.format(round(steady_state_temp - 273.15, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Estimated argument of perigee drift rate --> {} deg/yr'.format(round(sorted_data_set['arg_perigee_skew'][best_orbit_idx], 2)))
    print('Estimated delta v available to correct drift with EP --> {} km/s per year'.format(round(delta_v / 1000.0, 2)))
    print('Estimated delta v required to negate drift --> {} km/s per year'.format(round(delta_v_to_maintain / 1000.0, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Total active time (blackout reduction) --> {} %'.format(round(sorted_data_set['total_active_time'][best_orbit_idx], 2)))
    print('Total blackout time --> {} %'.format(round(100.0 * sorted_data_set['total_blackout_time'][best_orbit_idx] / study['duration'], 2)))
    print('Max active period duration --> {} hours'.format(round(sorted_data_set['max_active_time'][best_orbit_idx] / 3600.0, 2)))
    print('Max blackout period duration --> {} hours'.format(round(sorted_data_set['max_blackout_time'][best_orbit_idx], 2)))
    print('-----------------------------------------------------------------------------------------------------------')
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
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['arg_perigee_skew'], 500)
    plt.title('Argument Perigee Skew Rate [deg/yr]')
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

