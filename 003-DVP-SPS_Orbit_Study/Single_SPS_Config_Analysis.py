""""
17/05/2018
Author: Darian van Paridon

This script is for analyzing the access time for a lunar SPS in polar orbit,
targeting the south pole. Satellite flies at 90 degree inclination, with 90
degree argument of perigee

"""

from DVP_general_SPS_functions import *
import matplotlib.pyplot as plt


def main():
    start = convert_string_to_datetime(['2008', '07', '1', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    sps_access_raw = 'SPS(3k)_Circular_45Ntarget_Access.csv'
    sps_lighting_raw = 'SPS(3k)_Circular_45Ntarget_Lighting.csv'
    target_eclipse_raw = 'Target1-Eclipse-Edited.csv'
    sps_access = parse_csv_to_array(sps_access_raw, start)
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    target_eclipse = parse_csv_to_array(target_eclipse_raw, start)

    SPS_active_time = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)

main()
