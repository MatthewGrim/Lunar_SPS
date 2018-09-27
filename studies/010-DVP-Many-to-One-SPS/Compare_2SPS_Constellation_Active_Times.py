""""
25/06/2018
Author: Darian van Paridon

In this script the performance of a dual-SPS constellation, targeting the lunar south pole, is assessed as a function
of the relative initial true anomaly between the satellites.

"""

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


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
    constellation_total_active = []
    for i in range(1, len(true_anomalies)):
        constellation_active = combine_events(sps_active['0deg'], sps_active['{}deg'.format(true_anomalies[i])])
        constellation_total_active.append(np.sum(constellation_active[2]))
    ####################################################################################################################

    plt.plot(true_anomalies[1:], [100.0 * i / total_duration for i in constellation_total_active])
    plt.xlabel('Relative True Anomaly of Second SPS [deg]')
    plt.ylabel('Total Active Time [%]')
    plt.title("Comparison of Dual-SPS Constellation Performance")
    plt.show()


main()
