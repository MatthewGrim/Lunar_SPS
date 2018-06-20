""""
01/06/2018
Author: Darian van Paridon

This script takes the processed data for the Lunar SPS for south pole target orbit parametric scan, and applies
constraints relevant to the SPS and target systems to remove data points which do no the meet the requirements. Then
a weighted objective function is applied to determine the highest performing orbit within the feasible design space.
The resulting system performance metrics for the highest performing design are then printed.

The user can select from a small variety of parametrized transmitters and receivers to analyze the design space for
different combinations. The user can also toggle the weightings of the objective function to preference.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *


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
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name study
    study_name = 'SouthPole_IncrementedRes_Inertial'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    ####################################################################################################################

    # SELECT TRANSMITTER
    ####################################################################################################################
    # IPG YLS10000 (high power)
    wavelength = 1070e-9
    trans_radius = 0.25
    trans_power = 100e3
    trans_mass = 3600.0
    trans_eff = 0.35

    # IPG YLS-CUT (low power)
    # wavelength = 1070e-9
    # trans_radius = 0.3
    # trans_power = 15e3
    # trans_mass = 440.0
    # trans_eff = 0.35
    ####################################################################################################################

    # SELECT RECEIVER
    ####################################################################################################################
    # Team ITALIA AMALIA (intermediate)
    # rec_radius = 0.5
    # rec_efficiency = 0.40
    # operation_pwr = 100.0
    # hibernation_pwr = 7.0
    # battery_capacity = 100.0

    # ispace Sorato (miniature)
    rec_radius = 0.1
    rec_efficiency = 0.40
    operation_pwr = 21.5
    hibernation_pwr = 4.5
    battery_capacity = 38

    # NASA Curiosity (large)
    # rec_radius = 1
    # rec_efficiency = 0.40
    # operation_pwr = 270
    # hibernation_pwr = 23.5
    # battery_capacity = 1600.0
    ####################################################################################################################

    # SPECIFY DESIGN CONSTRAINTS
    ####################################################################################################################
    # Minimum pointing error of SPS system in radians
    point_error = 1e-6
    # Minimum reduction in overall blackout time in percent
    min_active_time = 0.0
    # Minimum power requirement at target in Watts
    min_power = operation_pwr
    # Maximum time rover can survive without recharging in hours
    rover_blackout_limit = battery_capacity / hibernation_pwr
    ####################################################################################################################

    # SELECT DESIGN OBJECTIVE WEIGHTINGS
    ####################################################################################################################
    # Define weights for generalized objective function - must sum to one
    # Relative importance of maximizing total active time / minimizing total blackout time
    weight_active_time = 0.4
    # Relative importance of maximizing energy delivered to target
    weight_energy = 0.2
    # Relative importance of minimizing duration of mean duration of individual blackout events
    weight_blackout_time = 0.2
    # Relative importance of minimizing orbital perturbations of initial orbit
    weight_perturbations = 0.2
    ####################################################################################################################

    # READ IN DATA FILES
    ####################################################################################################################
    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive_Inertial_Extended")
    total_blackout_time = read_data_from_file(stk_data_path, study_name, "TotalBlackout_Inertial_Extended")
    max_active_time = read_data_from_file(stk_data_path, study_name, "MaxActive_Inertial_Extended")
    max_blackout_time = read_data_from_file(stk_data_path, study_name, "MaxBlackout_Inertial_Extended")
    mean_range = read_data_from_file(stk_data_path, study_name, "MeanRange_Inertial_Extended")
    mean_active_time = read_data_from_file(stk_data_path, study_name, "MeanActive_Inertial_Extended")
    mean_blackout_time = read_data_from_file(stk_data_path, study_name, "MeanBlackout_Inertial_Extended")
    ####################################################################################################################

    # ENFORCE POINTING CONSTRAINTS
    ####################################################################################################################
    mean_link_efficiency = []
    for i in range(len(mean_range)):
        # Minimum beam radius as defined by pointing error
        min_beam_radius = rec_radius + (point_error * mean_range[i] * 1000.0)
        # Actual beam radius as defined by Gaussian beam divergence
        surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (mean_range[i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
        # Calculate mean link efficiency
        link_eff_temp = (rec_radius / surf_beam_radius) ** 2
        # Check pointing error constraint
        if surf_beam_radius < min_beam_radius:
            mean_link_efficiency.append(np.nan)
            # Apply effect of constraint to other data sets
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            max_active_time[i] = np.nan
            max_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan

        else:
            mean_link_efficiency.append(link_eff_temp)
    ####################################################################################################################

    # ENFORCE POWER CONSTRAINT
    ####################################################################################################################
    mean_power_received = []
    for i in range(len(mean_link_efficiency)):
        mean_power_temp = mean_link_efficiency[i] * rec_efficiency * trans_power
        if mean_power_temp < min_power:
            mean_power_received.append(np.nan)
            # Apply effect of constraint to other data sets
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            max_active_time[i] = np.nan
            max_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_link_efficiency[i] = np.nan
        else:
            mean_power_received.append(mean_power_temp)

    # Calculate total energy delivered to receiver
    total_energy = [i * j for i, j in zip(mean_power_received, total_active_time)]
    ####################################################################################################################

    # ENFORCING BLACKOUT DURATION CONSTRAINTS
    ####################################################################################################################
    # Remove data points for which blackout durations exceed the limit
    for i in range(len(max_blackout_time)):
        if max_blackout_time[i] > rover_blackout_limit * 3600.0:
            max_blackout_time[i] = np.nan

            # Apply constraint to all other data sets
            max_active_time[i] = np.nan
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_link_efficiency[i] = np.nan
            mean_power_received[i] = np.nan
            total_energy[i] = np.nan
    ####################################################################################################################

    # ENFORCE ACTIVE DURATION CONSTRAINT
    ####################################################################################################################
    # Remove data points for which the overall blackout time is not sufficiently reduced
    for j in range(len(total_active_time)):
        if (100.0 * total_active_time[j] / total_duration) < min_active_time:
            total_active_time[j] = np.nan

            # Apply constraint to all other data sets
            max_active_time[i] = np.nan
            total_active_time[i] = np.nan
            total_blackout_time[i] = np.nan
            mean_active_time[i] = np.nan
            mean_blackout_time[i] = np.nan
            mean_link_efficiency[i] = np.nan
            mean_power_received[i] = np.nan
            total_energy[i] = np.nan
    ####################################################################################################################

    # EVALUATE REQUIRED POWER GENERATOR METRICS
    ####################################################################################################################
    # Power generator parameter - Stretched lens array SquareRigger platform
    generator_eff = 0.4
    generator_spec_pwr = 300
    # Power required from generator
    generator_pwr = trans_power / (generator_eff * trans_eff)
    # Mass of generator
    generator_mass = generator_pwr / generator_spec_pwr
    ####################################################################################################################

    # EVALUATE HEAT GENERATED
    ####################################################################################################################
    # Heat absorbed by generator
    generator_heat = (1 - generator_eff) * generator_pwr
    # Heat generated by transmitter
    transmitter_heat = (1 - trans_eff) * trans_power
    total_heat = transmitter_heat + generator_heat
    ####################################################################################################################

    # ESTIMATE MAGNITUDE OF ORBITAL PERTURBATIONS
    ####################################################################################################################
    # Orbital perturbations on argument of perigee and eccentricity due to Earth's gravity
    perturbations = calculate_orbital_perturbations(semi_maj_axis, eccentricity)

    arg_perigee_pert_norm = abs(perturbations[0]) / np.max(abs(perturbations[0]))
    ecc_pert_norm = perturbations[1] / np.max(perturbations[1])
    combined_pert = (arg_perigee_pert_norm + ecc_pert_norm) / max(arg_perigee_pert_norm + ecc_pert_norm)
    ####################################################################################################################

    # PARSING DATA
    ####################################################################################################################
    # Reorganize the data lists into 2D arrays
    total_active_times_sorted = sort_incremented_resolution_data(orbit_data, total_active_time)
    total_blackout_times_sorted = sort_incremented_resolution_data(orbit_data, total_blackout_time)
    max_active_times_sorted = sort_incremented_resolution_data(orbit_data, max_active_time)
    max_blackout_times_sorted = sort_incremented_resolution_data(orbit_data, max_blackout_time)
    mean_blackout_times_sorted = sort_incremented_resolution_data(orbit_data, mean_blackout_time)
    link_efficiency_sorted = sort_incremented_resolution_data(orbit_data, mean_link_efficiency)
    power_received_sorted = sort_incremented_resolution_data(orbit_data, mean_power_received)
    total_energy_sorted = sort_incremented_resolution_data(orbit_data, total_energy)
    combined_pert_sorted = sort_incremented_resolution_data(orbit_data, combined_pert)

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

    # DEFINE OBJECTIVE FUNCTION and SELECT BEST SOLUTION
    ####################################################################################################################
    # Normalize design variables for direct comparison in objective function
    active_times_normalized = total_active_times_sorted / np.nanmax(total_active_times_sorted)
    total_energy_normalized = total_energy_sorted / np.nanmax(total_energy_sorted)
    blackout_times_normalized = mean_blackout_times_sorted / np.nanmax(mean_blackout_times_sorted)

    # Evaluate objective function
    objective = np.array([weight_active_time * a + weight_energy * b + weight_blackout_time * c + weight_perturbations * d
                          for a, b, c, d in zip(active_times_normalized,
                                             total_energy_normalized,
                                             (1 - blackout_times_normalized),
                                             (1 - combined_pert_sorted))])

    # Find best orbit according to weighted objective function
    from numpy import unravel_index
    best_orbit_idx = unravel_index(np.nanargmax(objective), objective.shape)
    best_perigee = unique_perigees[best_orbit_idx[0]]
    best_apogee = unique_apogees[best_orbit_idx[1]]
    ####################################################################################################################

    # DISPLAY RESULTS
    ####################################################################################################################
    # Print out performance results for optimal orbit
    print('Combined mass of generator and transmitter --> {} kg'.format(round(generator_mass + trans_mass, 2)))
    print('Combined heat load of generator and transmitter --> {} kW'.format(round(total_heat / 1000.0, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Optimal orbit altitudes --> Perigee: {} km, Apogee: {} km'.format(round(best_perigee - r_moon, 2), round(best_apogee - r_moon, 2)))
    print('Total active time (blackout reduction) --> {} hours, or {}%'.format(round(total_active_times_sorted[best_orbit_idx] / 3600.0, 2), round(
        100.0 * total_active_times_sorted[best_orbit_idx] / total_duration, 2)))
    print('Max active period duration --> {} hours'.format(round(max_active_times_sorted[best_orbit_idx] / 3600.0, 2)))
    print('Total blackout time --> {} hours, or {} %'.format(round(total_blackout_times_sorted[best_orbit_idx] / 3600.0, 2), round(
        100.0 * total_blackout_times_sorted[best_orbit_idx] / total_duration, 2)))
    print('Max blackout period duration --> {} hours'.format(round(max_blackout_times_sorted[best_orbit_idx] / 3600.0, 2)))
    print('-----------------------------------------------------------------------------------------------------------')
    print('Link efficiency --> {}%'.format(round(link_efficiency_sorted[best_orbit_idx] * 100.0, 5)))
    print('Power received per access period --> {} W'.format(round(power_received_sorted[best_orbit_idx], 2)))
    print('Total energy transferred --> {} MJ per year'.format(round(total_energy_sorted[best_orbit_idx] / 2e6, 2)))

    # Plot weighted objective function
    plt.figure(1)
    plt.contourf(apogee_altitudes, perigee_altitudes, objective, 500)
    plt.colorbar()
    plt.scatter(best_apogee - r_moon, best_perigee - r_moon, marker='x')
    textstr = 'Mean Power Received: {} W\nTotal Active Time Weighting: {}%\nMean Blackout Time Weighting {}%\nTotal Energy Weighting: {}%\nOrbital Perturbation Weighting: {}%'\
        .format(round(power_received_sorted[best_orbit_idx], 2), weight_active_time * 100.0, weight_blackout_time * 100.0, weight_energy * 100.0, weight_perturbations * 100.0)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(100, 4000, textstr, fontsize=7, verticalalignment='top', bbox=props)
    plt.title("Weighted Objective Function")
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.show()


main()
