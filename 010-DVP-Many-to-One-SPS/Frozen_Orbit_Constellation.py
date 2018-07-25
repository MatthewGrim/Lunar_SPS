""""
25/06/2018
Author: Darian van Paridon

This script is a first attempt addressing the many-to-one problem for solar power satellites. It contains analysis for
determining the total reduction in blackout time at the target for a particular constellation of SPS.

"""

from DVP_general_SPS_functions import *
from DVP_Programmatic_Functions import *
from SPS_Constrained_DesignFunctions import *


def calculate_mean_link_eff(trans_radius, args):

    mean_range = args[0]
    transmitter = args[1]
    rover = args[2]

    mean_link_eff = []
    for i in range(len(mean_range)):
        surf_beam_radius = trans_radius * np.sqrt(1 + (transmitter['wavelength'] * (mean_range[i] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
        mean_link_eff.append((rover['rec_radius'] / surf_beam_radius) ** 2)

    return 1 - np.mean(mean_link_eff)


def main():

    # INITIALIZATION
    ####################################################################################################################
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Name study
    study_name = 'FrozenOrbit'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    true_anomalies = [0, 120, 180, 240]
    ####################################################################################################################

    # IMPORT AND PROCESS DATA
    ####################################################################################################################
    sps_lighting = {}
    sps_access = {}
    sps_active = {}
    sps_range = {}
    target_blackout = {}

    target_lighting = parse_csv_to_array('{}/FrozenOrbit_Target_Lighting.csv'.format(stk_data_path), start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    for i in true_anomalies:
        sps_lighting['{}deg'.format(i)] = parse_csv_to_array('{}/FrozenOrbit_4400sma_lighting_{}deg.csv'.format(stk_data_path, i), start)
        sps_range['{}deg'.format(i)] = import_range_data_statistics('FrozenOrbit_4400sma_range_{}deg'.format(i), stk_data_path)
        sps_access['{}deg'.format(i)] = parse_csv_to_array('{}/FrozenOrbit_4400sma_access_{}deg.csv'.format(stk_data_path, i), start)
        sps_active['{}deg'.format(i)] = determine_SPS_active_time(sps_lighting['{}deg'.format(i)], target_eclipse, sps_access['{}deg'.format(i)])
        target_blackout['{}deg'.format(i)] = determine_blackout_data(sps_active['{}deg'.format(i)], target_eclipse, total_duration)
    ####################################################################################################################

    # DETERMINE CONSTELLATION COVERAGE PERFORMANCE
    ####################################################################################################################
    double_sps = combine_events(sps_active['0deg'], sps_active['180deg'])

    triple_sps_temp = combine_events(sps_active['0deg'], sps_active['120deg'])
    triple_sps = combine_events(triple_sps_temp, sps_active['240deg'])

    single_sps_blackout = determine_blackout_data(sps_active['0deg'], target_eclipse, total_duration)
    double_sps_blackout = determine_blackout_data(double_sps, target_eclipse, total_duration)
    triple_sps_blackout = determine_blackout_data(triple_sps, target_eclipse, total_duration)
    print('\n')
    print('Single SPS active time: {} %'.format(round(100.0 * np.sum(sps_active['0deg'][2]) / total_duration, 4)))
    print('Single SPS blackout time: {} %'.format(round(100.0 * np.sum(single_sps_blackout[2]) / total_duration, 4)))
    print('\n')
    print('Double SPS active time: {} %'.format(round(100.0 * np.sum(double_sps[2]) / total_duration, 4)))
    print('Double SPS blackout time: {} %'.format(round(100.0 * np.sum(double_sps_blackout[2]) / total_duration, 4)))
    print('\n')
    print('Triple SPS active time: {} %'.format(round(100.0 * np.sum(triple_sps[2]) / total_duration, 4)))
    print('Triple SPS blackout time: {} %'.format(round(100.0 * np.sum(triple_sps_blackout[2]) / total_duration, 4)))
    ####################################################################################################################

    # DETERMINE ORBIT LINK PERFORMANCE
    ####################################################################################################################
    # Select transmitter and rover
    transmitter = trans_metrics('100kW')
    rover = rover_metrics('curiosity')

    # Optimize transmitter size
    from scipy.optimize import minimize_scalar
    args = [sps_range['0deg'][2], transmitter, rover]
    optimum = minimize_scalar(calculate_mean_link_eff, bounds=(0, 2.5), method='bounded', args=args)
    transmitter['radius'] = optimum.x
    print('Optimal aperture radius: {} cm'.format(round(optimum.x * 100.0, 2)))

    # Calculate mean link efficiency and power
    mean_link_eff = []
    for i in range(len(sps_range['0deg'][2])):
        surf_beam_radius = transmitter['radius'] * np.sqrt(1 + (transmitter['wavelength'] * (sps_range['0deg'][2][i] * 1000.0) / (np.pi * transmitter['radius'] ** 2)) ** 2)
        mean_link_eff.append((rover['rec_radius'] / surf_beam_radius) ** 2)

    mean_power = [i * rover['rec_efficiency'] * transmitter['power'] for i in mean_link_eff]

    # Calculate surface beam size and assess required pointing accuracy
    # Minimum beam radius as defined by pointing error
    # Actual beam radius as defined by Gaussian beam divergence
    surf_beam_radius = [transmitter['radius'] * np.sqrt(1 + (
                transmitter['wavelength'] * (i * 1000.0) / (
                    np.pi * transmitter['radius'] ** 2)) ** 2) for i in sps_range['0deg'][1]]

    pointing_error = [(i - rover['rec_radius']) / (j * 1000.0) for i,j in zip(surf_beam_radius, sps_range['0deg'][1])]
    print('Mean acceptable pointing error: {} rad'.format(round(np.mean(pointing_error), 10)))
    plt.plot([i * 1e6 for i in pointing_error])
    plt.title('Acceptable pointing error for mean range of every access period.')
    plt.ylabel('1e-6 rad')
    plt.xlabel('Access event')
    plt.show()

    print('\n')
    print('Optimum transmitter radius: {} cm'.format(round(100.0 * transmitter['radius'])))
    print('Mean link efficiency: {} %'.format(round(np.mean(mean_link_eff) * 100.0, 2)))
    print('Mean power delivered: {} W'.format(round(np.mean(mean_power), 2)))
    ####################################################################################################################

    plt.subplot(311)
    plt.bar([i / 86400.0 for i in single_sps_blackout[0]], [j / 3600.0 for j in single_sps_blackout[2]])
    plt.title('One SPS')
    plt.ylabel('Blackout Duration [h]')
    plt.subplot(312)
    plt.bar([i / 86400.0 for i in double_sps_blackout[0]], [j / 3600.0 for j in double_sps_blackout[2]])
    plt.title('Two SPS')
    plt.ylabel('Blackout Duration [h]')
    plt.subplot(313)
    plt.bar([i / 86400.0 for i in triple_sps_blackout[0]], [j / 3600.0 for j in triple_sps_blackout[2]])
    plt.title('Three SPS')
    plt.ylabel('Blackout Duration [h]')
    plt.xlabel('Days')
    plt.show()


main()
