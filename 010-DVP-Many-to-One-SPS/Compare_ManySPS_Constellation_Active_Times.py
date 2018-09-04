""""
25/06/2018
Author: Darian van Paridon

This script is a first attempt addressing the many-to-one problem for solar power satellites. It contains analysis for
determining the total reduction in blackout time at the target for a particular constellation of SPS.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import *


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
    study_name = 'SouthPole_ManyToOne\Test2'

    # Set file path for data
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    true_anomalies = [0, 60, 72, 90, 120, 144, 180, 216, 240, 270, 288, 300]
    ####################################################################################################################

    # IMPORT AND PROCESS DATA
    ####################################################################################################################
    sps_lighting = {}
    sps_access = {}
    sps_active = {}
    target_blackout = {}

    target_lighting = parse_csv_to_array('{}/Target_lighting.csv'.format(stk_data_path), start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    for i in true_anomalies:
        sps_lighting['{}deg'.format(i)] = parse_csv_to_array('{}/SPS_{}deg_lighting.csv'.format(stk_data_path, i), start)
        sps_access['{}deg'.format(i)] = parse_csv_to_array('{}/SPS_{}deg_access.csv'.format(stk_data_path, i), start)
        sps_active['{}deg'.format(i)] = determine_SPS_active_time(sps_lighting['{}deg'.format(i)], target_eclipse, sps_access['{}deg'.format(i)])
        target_blackout['{}deg'.format(i)] = determine_blackout_data(sps_active['{}deg'.format(i)], target_eclipse, total_duration)
    ####################################################################################################################

    # DETERMINE CONSTELLATION PERFORMANCE
    ####################################################################################################################
    number_of_sps = 2
    for i in range(number_of_sps):
        true_anomalies[i] = 0.0 + (i * 360.0 / number_of_sps)

    double_sps = combine_events(sps_active['0deg'], sps_active['180deg'])

    triple_sps = combine_events(sps_active['0deg'], sps_active['120deg'])
    triple_sps = combine_events(triple_sps, sps_active['240deg'])

    quad_sps1 = combine_events(sps_active['0deg'], sps_active['90deg'])
    quad_sps2 = combine_events(sps_active['180deg'], sps_active['270deg'])
    quad_sps = combine_events(quad_sps1, quad_sps2)

    quint_sps1 = combine_events(sps_active['0deg'], sps_active['72deg'])
    quint_sps2 = combine_events(sps_active['144deg'], sps_active['216deg'])
    quint_sps3 = combine_events(quint_sps1, quint_sps2)
    quint_sps = combine_events(quint_sps3, sps_active['288deg'])

    ses_sps1 = combine_events(sps_active['0deg'], sps_active['60deg'])
    ses_sps2 = combine_events(sps_active['120deg'], sps_active['180deg'])
    ses_sps3 = combine_events(sps_active['240deg'], sps_active['300deg'])
    ses_sps4 = combine_events(ses_sps1, ses_sps2)
    ses_sps = combine_events(ses_sps3, ses_sps4)

    single_sps_blackout = determine_blackout_data(sps_active['0deg'], target_eclipse, total_duration)
    double_sps_blackout = determine_blackout_data(double_sps, target_eclipse, total_duration)
    triple_sps_blackout = determine_blackout_data(triple_sps, target_eclipse, total_duration)
    quad_sps_blackout = determine_blackout_data(quad_sps, target_eclipse, total_duration)
    quint_sps_blackout = determine_blackout_data(quint_sps, target_eclipse, total_duration)
    ses_sps_blackout = determine_blackout_data(ses_sps, target_eclipse, total_duration)

    print('Single SPS active time: {} %'.format(round(100.0 * np.sum(sps_active['0deg'][2]) / total_duration, 4)))
    print('Single SPS blackout time: {} %'.format(round(100.0 * np.sum(single_sps_blackout[2]) / total_duration, 4)))

    print('Double SPS active time: {} %'.format(round(100.0 * np.sum(double_sps[2]) / total_duration, 4)))
    print('Double SPS blackout time: {} %'.format(round(100.0 * np.sum(double_sps_blackout[2]) / total_duration, 4)))

    print('Triple SPS active time: {} %'.format(round(100.0 * np.sum(triple_sps[2]) / total_duration, 4)))
    print('Triple SPS blackout time: {} %'.format(round(100.0 * np.sum(triple_sps_blackout[2]) / total_duration, 4)))

    print('Quadruple SPS active time: {} %'.format(round(100.0 * np.sum(quad_sps[2]) / total_duration, 4)))
    print('Quadruple SPS blackout time: {} %'.format(round(100.0 * np.sum(quad_sps_blackout[2]) / total_duration, 4)))

    print('Quintuple SPS active time: {} %'.format(round(100.0 * np.sum(quint_sps[2]) / total_duration, 4)))
    print('Quintuple SPS blackout time: {} %'.format(round(100.0 * np.sum(quint_sps_blackout[2]) / total_duration, 4)))

    print('Sestuple SPS active time: {} %'.format(round(100.0 * np.sum(ses_sps[2]) / total_duration, 4)))
    print('Sestuple SPS blackout time: {} %'.format(round(100.0 * np.sum(ses_sps_blackout[2]) / total_duration, 4)))
    ####################################################################################################################

    plt.subplot(321)
    plt.bar([i / 86400.0 for i in single_sps_blackout[0]], [j / 3600.0 for j in single_sps_blackout[2]])
    plt.title('One SPS')
    plt.subplot(322)
    plt.bar([i / 86400.0 for i in double_sps_blackout[0]], [j / 3600.0 for j in double_sps_blackout[2]])
    plt.title('Two SPS')
    plt.subplot(323)
    plt.bar([i / 86400.0 for i in triple_sps_blackout[0]], [j / 3600.0 for j in triple_sps_blackout[2]])
    plt.title('Three SPS')
    plt.subplot(324)
    plt.bar([i / 86400.0 for i in quad_sps_blackout[0]], [j / 3600.0 for j in quad_sps_blackout[2]])
    plt.title('Four SPS')
    plt.subplot(325)
    plt.bar([i / 86400.0 for i in quint_sps_blackout[0]], [j / 3600.0 for j in quint_sps_blackout[2]])
    plt.title('Five SPS')
    plt.subplot(326)
    plt.bar([i / 86400.0 for i in ses_sps_blackout[0]], [j / 3600.0 for j in ses_sps_blackout[2]])
    plt.title('Six SPS')
    plt.xlabel('Days')
    plt.show()


main()
