""""
27/06/2018
Author: Darian van Paridon

This script is a simple analysis used to approximate the number of SPS in a constellation that would be required
to eliminate blackouts at a lunar south pole target.

"""
from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import *


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

    # READ IN DATA FILES
    ####################################################################################################################
    total_active_time = read_data_from_file(stk_data_path, study_name, "TotalActive_Inertial_Extended")
    # Import target illumination events
    target_lighting_raw = '{}\DVP_{}_Target_Lighting.csv'.format(stk_data_path, study_name)
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)
    total_target_eclipse = np.sum(target_eclipse[2])
    ####################################################################################################################

    # Calculate SPS constellation size required
    number_of_sps = [int(math.ceil(total_target_eclipse / i)) for i in total_active_time]
    # Determine the number of simulations necessary to append current data set
    number_of_sims = np.sum(number_of_sps) - np.size(number_of_sps)
    print('Number of simulations: {}'.format(number_of_sims))

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

    # Sort data
    number_of_sps = sort_incremented_resolution_data(orbit_data, number_of_sps)

    # Make plot
    plt.figure(1)
    plt.contourf(apogee_altitudes, perigee_altitudes, number_of_sps, 32)
    plt.colorbar()
    plt.title("Number of SPS to Eliminate Blackouts")
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.show()


main()
