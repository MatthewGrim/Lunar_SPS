""""
11/05/2018
Author: Darian van Paridon

This script is for verifying the optimal apogee altitude for a lunar SPS mission, which Brandhorst claims
is 30000 km in the paper "A solar electric propulsion mission for lunar power beaming"

"""

import matplotlib.pyplot as plt
from DVP_general_SPS_functions import *
from general_functions import convert_string_to_datetime


def scan_apogee_altitudes(apogees, start):
    # This function analyzes the availability of the SPS to a lunar target at 45 N as a function of the
    # altitude of the orbital apogee

    total_active_time = np.zeros(len(apogees))

    raw_eclipse_target = 'Target1-Eclipse-Edited.csv'
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)

    for i in range(len(apogees)):
        raw_sunlight = "Lighting({}k).csv".format(apogees[i])
        raw_access = "Access({}k).csv".format(apogees[i])

        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)
        print('\n')
        print('Apogee altitude {}000 km'.format(apogees[i]))
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_target, access_times)

    return total_active_time


def main():

    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    apogees = (1, 2, 5, 10, 15, 20, 25, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55)
    active_times = scan_apogee_altitudes(apogees, start)
    apogees_scaled = [i * 1000 for i in apogees]
    active_times_scaled = [i / 3600.0 for i in active_times]
    plt.plot(apogees_scaled, active_times_scaled)
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Hours of Access (Over 2 Years)')
    plt.show()


main()
