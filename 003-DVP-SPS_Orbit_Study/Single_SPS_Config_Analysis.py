""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the active events for a lunar SPS in any given orbit, where the access to the target,
lighting events, and target eclipse events data reports have been extracted from the STK simulation.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def main():
    # Initialize start and end times of simulation
    start = convert_string_to_datetime(['2008', '07', '1', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Select data reports for simulation/configuration
    sps_access_raw = 'SPS(3k)_Circular_45Ntarget_Access.csv'
    sps_lighting_raw = 'SPS(3k)_Circular_45Ntarget_Lighting.csv'
    target_eclipse_raw = 'Target1-Eclipse-Edited.csv'

    # Parse data reports and import data
    sps_access = parse_csv_to_array(sps_access_raw, start)
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    eclipse_target = parse_csv_to_array(target_eclipse_raw, start)

    # Calculate SPS active events
    sps_active = determine_SPS_active_time(sps_lighting, eclipse_target, sps_access)

    # Plot active events
    plt.figure(1)
    plt.bar([i / 86400.0 for i in sps_active[0]], [j / 3600.0 for j in sps_active[2]])
    plt.xlabel('Days')
    plt.ylabel('Active Event Duration [h]')
    plt.title('SPS Active Events')
    plt.show()


main()
