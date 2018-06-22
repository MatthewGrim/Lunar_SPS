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


def main():

    # INITIALIZATION
    ####################################################################################################################
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Set bounds on parametric scan
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Get orbit data set
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Name study
    study_name = 'SouthPole_IncrementedRes_Inertial'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    ####################################################################################################################

    # SELECT RECEIVER
    ####################################################################################################################
    rover = rover_metrics('sorato')
    ####################################################################################################################

    # SELECT TRANSMITTER
    ####################################################################################################################
    # Optimize transmitter radius to maximize link efficiency
    optimum = optimize_link_efficiency()

    # IPG YLS10000 (high power)
    # wavelength = 1070e-9
    # trans_radius = optimum.x
    # trans_power = 100e3
    # trans_mass = 3600.0
    # trans_eff = 0.35

    # IPG YLS-CUT (low power)
    wavelength = 1070e-9
    trans_radius = optimum.x
    trans_power = 15e3
    trans_mass = 440.0
    trans_eff = 0.35
    ####################################################################################################################

    # SPECIFY DESIGN CONSTRAINTS
    ####################################################################################################################
    constraints, active_constraints = sps_constraints(rover)
    ####################################################################################################################

    # READ IN DATA FILES
    ####################################################################################################################
    data_set = {"total_active_time": [], 'total_blackout_time': [], 'max_active_time': [],
                'max_blackout_time': [], 'mean_range': [], 'mean_active_time': [], 'mean_blackout_time': []}

    data_set['total_active_time'] = read_data_from_file(stk_data_path, study_name, "TotalActive_Inertial_Extended")
    data_set['total_blackout_time'] = read_data_from_file(stk_data_path, study_name, "TotalBlackout_Inertial_Extended")
    data_set['max_active_time'] = read_data_from_file(stk_data_path, study_name, "MaxActive_Inertial_Extended")
    data_set['max_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MaxBlackout_Inertial_Extended")
    data_set['mean_range'] = read_data_from_file(stk_data_path, study_name, "MeanRange_Inertial_Extended")
    data_set['mean_active_time'] = read_data_from_file(stk_data_path, study_name, "MeanActive_Inertial_Extended")
    data_set['mean_blackout_time'] = read_data_from_file(stk_data_path, study_name, "MeanBlackout_Inertial_Extended")
    ####################################################################################################################

    # ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS
    ####################################################################################################################
    # Orbital perturbations on argument of perigee [0], eccentricity [1], and inclination [2]
    perturbations = calculate_orbital_perturbations(semi_maj_axis, eccentricity)
    # Calculate skew in argument of perigee in degrees per year
    data_set['arg_perigee_skew'] = [i * (365.0 * 24.0 * 3600.0) * 180.0 / np.pi for i in perturbations[0]]
    ####################################################################################################################

    # CALCULATE LINK EFFICIENCY, MEAN POWER, and TOTAL ENERGY DELIVERED
    # POINTING ERROR CONSTRAINT ALSO APPLIED
    ####################################################################################################################
    data_set['mean_link_efficiency'] = []
    data_set['mean_power_received'] = []
    for i in range(len(data_set['mean_range'])):
        # Minimum beam radius as defined by pointing error
        min_beam_radius = rover['rec_radius'] + (constraints['point_error'] * data_set['mean_range'][i] * 1000.0)
        # Actual beam radius as defined by Gaussian beam divergence
        surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (data_set['mean_range'][i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
        # Calculate link efficiency
        if surf_beam_radius < min_beam_radius:
            if active_constraints['point_error'] == 1:
                data_set['mean_link_efficiency'].append(np.nan)
                data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * trans_power)
                # apply constraint all other data lists
                for j in data_set:
                    data_set[j][i] = np.nan
            else:
                data_set['mean_link_efficiency'].append(1.0)
                data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * trans_power)
        else:
            data_set['mean_link_efficiency'].append((rover['rec_radius'] / surf_beam_radius) ** 2)
            data_set['mean_power_received'].append(data_set['mean_link_efficiency'][i] * rover['rec_efficiency'] * trans_power)
    # Calculate total energy delivered to receiver
    data_set['total_energy'] = []
    data_set['total_energy'] = [i * j for i, j in zip(data_set['mean_power_received'], data_set['total_active_time'])]
    ####################################################################################################################

    # ENFORCE CONSTRAINTS
    ####################################################################################################################
    # Convert data to appropriate units for applying constraints
    data_set['max_blackout_time'] = [i / 3600.0 for i in data_set['max_blackout_time']]
    data_set['total_active_time'] = [100.0 * i / total_duration for i in data_set['total_active_time']]
    # Remove data points for which not enough power is delivered on average
    if active_constraints['min_power'] == 1:
        data_set = enforce_constraints(data_set, 'mean_power_received', constraints, 'min_power', 'min')
    else:
        pass
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
    # Remove data points for which the overall blackout time is not sufficiently reduced
    if active_constraints['max_arg_perigee_skew'] == 1:
        data_set = enforce_constraints(data_set, 'arg_perigee_skew', constraints, 'max_perigee_skew', 'max')
    else:
        pass
    ####################################################################################################################

    # EVALUATE REQUIRED POWER GENERATOR METRICS
    ####################################################################################################################
    # Power generator parameters - Stretched lens array SquareRigger platform
    generator_eff = 0.4
    generator_spec_pwr = 300
    # Power required from generator
    generator_pwr = trans_power / trans_eff
    # Mass of generator
    generator_mass = generator_pwr / generator_spec_pwr
    ####################################################################################################################

    # EVALUATE HEAT GENERATED AND STEADY STATE TEMPERATURE
    ####################################################################################################################
    # Heat absorbed by generator
    generator_heat = (1 - generator_eff) * generator_pwr
    # Heat generated by transmitter
    transmitter_heat = (1 - trans_eff) * trans_power
    # Total heat in system
    total_heat = transmitter_heat + generator_heat
    # Steady state temperature
    solar_irradiance = 1367.0
    emissivity = 0.9
    steady_state_temp = (solar_irradiance * (1 - generator_eff * trans_eff) / (2 * emissivity * 5.67e-8)) ** 0.25
    ####################################################################################################################

    # PARSING DATA
    ####################################################################################################################
    # Reorganize the data lists into 2D arrays
    sorted_data_set = sort_data_lists(data_set, orbit_data)

    # Extract unique perigees and apogees tested for plotting
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    r_moon = 1737.0
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

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
    ####################################################################################################################

    # ESTIMATE SPS BATTERY MASS
    ####################################################################################################################
    sps_battery_capacity = trans_power * sorted_data_set['mean_blackout_time'][best_orbit_idx] / 3600.0 * trans_eff
    lipo_specific_power = 270.0
    sps_battery_mass = sps_battery_capacity / lipo_specific_power
    ####################################################################################################################

    # EVALUATE FLUX AT RECEIVER
    ####################################################################################################################
    surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (sorted_data_set['mean_range'][best_orbit_idx] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
    target_flux = trans_power / (np.pi * surf_beam_radius ** 2)
    ####################################################################################################################

    # DISPLAY RESULTS
    ####################################################################################################################
    # Print out performance results for optimal orbit
    print('-----------------------------------------------------------------------------------------------------------')
    print('Combined mass of generator and transmitter --> {} kg'.format(round(generator_mass + trans_mass, 2)))
    print('Combined heat load of generator and transmitter --> {} kW'.format(round(total_heat / 1000.0, 2)))
    print('Steady state temperature --> {} Celsius'.format(round(steady_state_temp - 273.15, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Optimal orbit altitudes --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Estimated argument of perigee slew rate --> {} deg/yr'.format(round(sorted_data_set['arg_perigee_skew'][best_orbit_idx], 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Total active time (blackout reduction) --> {}%'.format(round(sorted_data_set['total_active_time'][best_orbit_idx], 2)))
    print('Total blackout time --> {} %'.format(round(100.0 * sorted_data_set['total_blackout_time'][best_orbit_idx] / total_duration, 2)))
    print('Max active period duration --> {} hours'.format(round(sorted_data_set['max_active_time'][best_orbit_idx] / 3600.0, 2)))
    print('Max blackout period duration --> {} hours'.format(round(sorted_data_set['max_blackout_time'][best_orbit_idx], 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Mean blackout period duration --> {} hours'.format(round(sorted_data_set['mean_blackout_time'][best_orbit_idx] / 3600.0, 2)))
    print('Approximate battery mass required to eliminate mean blackout --> {} kg'.format(round(sps_battery_mass, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Link efficiency --> {}%'.format(round(sorted_data_set['mean_link_efficiency'][best_orbit_idx] * 100.0, 5)))
    print('Power received per access period --> {} W'.format(round(sorted_data_set['mean_power_received'][best_orbit_idx], 2)))
    print('Mean flux at receiver --> {} W/m2'.format(round(target_flux, 2)))
    print('Total energy transferred --> {} MJ per year'.format(round(sorted_data_set['total_energy'][best_orbit_idx] / 2e6, 2)))

    # Plot weighted objective function
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
    plt.contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['mean_link_efficiency'] * 100.0, 500)
    plt.title('Mean Link Efficiency [%]')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.scatter(best_apogee - r_moon, best_perigee - r_moon, marker='x')
    plt.show()


main()
