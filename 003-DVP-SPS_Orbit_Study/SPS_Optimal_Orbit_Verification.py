""""
11/05/2018
Author: Darian van Paridon

This script is for verifying Brandhorst's claim that the optimal apogee altitude for the proposed SPS constellation
(in " A solar electric propulsion mission for lunar power beaming ") is 30000 km. In that paper, the perigee altitude
is 500 km. This hypothesis is tested as well.

The verification is addressed by fixing the above mentioned perigee and apogee altitudes, and varying the other altitude.
The highest performing orbit is determined by calculating the total active time and maximum blackout event duration. The
optimal orbit would have the highest total active, and lowest maximum blackout times.

"""

from DVP_general_SPS_functions import *
from general_functions import convert_string_to_datetime


def scan_apogee_altitudes(apogees, start, total_duration):
    # This function analyzes the availability of the SPS to a lunar target at 45 N as a function of the
    # altitude of the orbital apogee

    total_active_time = []
    max_blackout_time = []

    raw_eclipse_target = 'Target1-Eclipse-Edited.csv'
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)

    for i in range(len(apogees)):
        raw_sunlight = "Lighting({}k).csv".format(apogees[i])
        raw_access = "Access({}k).csv".format(apogees[i])
        sps_access = parse_csv_to_array(raw_access, start)
        sps_lighting = parse_csv_to_array(raw_sunlight, start)
        print('\n')
        print('Apogee {}000 km, Perigee 500 km'.format(apogees[i]))
        sps_active = determine_SPS_active_time(sps_lighting, eclipse_target, sps_access)
        total_active_time.append(np.sum(sps_active[2]))
        print('Maximum eclipse without SPS: {} hrs'.format(round(max(eclipse_target[2]) / 3600.0, 2)))
        sps_blackout = determine_blackout_data(sps_access, eclipse_target, total_duration)
        max_blackout_time.append(max(sps_blackout[2]))

    return total_active_time, max_blackout_time


def scan_perigee_altitudes(perigees, start, total_duration):
    # This function analyzes the availability of the SPS to a lunar target at 45 N as a function of the
    # altitude of the orbital apogee

    total_active_time = []
    max_blackout_time = []

    raw_eclipse_target = 'Target1-Eclipse-Edited.csv'
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)

    for i in range(len(perigees)):
        raw_sunlight = "Lighting({},30k).csv".format(perigees[i])
        raw_access = "Access({},30k).csv".format(perigees[i])
        sps_access = parse_csv_to_array(raw_access, start)
        sps_lighting = parse_csv_to_array(raw_sunlight, start)

        print('\n')
        print('Perigee {} km, Apogee 30000 km'.format(perigees[i]))
        sps_active = determine_SPS_active_time(sps_lighting, eclipse_target, sps_access)
        total_active_time.append(np.sum(sps_active[2]))
        print('Maximum eclipse without SPS: {} hrs'.format(round(max(eclipse_target[2]) / 3600.0, 2)))
        sps_blackout = determine_blackout_data(sps_access, eclipse_target, total_duration)
        max_blackout_time.append(max(sps_blackout[2]))

    return total_active_time, max_blackout_time


def main():

    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()
    apogees = (1, 2, 5, 10, 15, 20, 25, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 52.5, 55)
    perigees = (100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000)

    active_times_varying_apogee, blackout_times_varying_apogee = scan_apogee_altitudes(apogees, start, total_duration)
    active_times_varying_perigee, blackout_times_varying_perigee = scan_perigee_altitudes(perigees, start, total_duration)

    apogees_scaled = [i * 1000 for i in apogees]
    active_apogee_scaled = [i / 3600.0 for i in active_times_varying_apogee]
    active_perigee_scaled = [i / 3600.0 for i in active_times_varying_perigee]
    blackout_apogee_scaled = [i / 3600.0 for i in blackout_times_varying_apogee]
    blackout_perigee_scaled = [i / 3600.0 for i in blackout_times_varying_perigee]

    plt.figure(1)
    plt.subplot(221)
    plt.plot(perigees, active_perigee_scaled)
    plt.title('Apogee Fixed 30000 km')
    plt.xlabel('Perigee Altitude [km]')
    plt.ylabel('Total Active Time [h]')
    plt.subplot(222)
    plt.plot(perigees, blackout_perigee_scaled)
    plt.title('Apogee Fixed 30000 km')
    plt.xlabel('Perigee Altitude [km]')
    plt.ylabel('Max Blackout Duration [h]')
    plt.subplot(223)
    plt.plot(apogees_scaled, active_apogee_scaled)
    plt.title('Perigee Fixed 500 km')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Total Active Time [h]')
    plt.subplot(224)
    plt.plot(apogees_scaled, blackout_apogee_scaled)
    plt.title('Perigee Fixed 500 km')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Max Blackout Duration [h]')
    plt.show()


main()
