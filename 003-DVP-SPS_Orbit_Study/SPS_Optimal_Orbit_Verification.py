""""
11/05/2018
Author: Darian van Paridon

This script is for verifying the optimal apogee altitude for a lunar SPS mission, which Brandhorst claims
is 30000 km in the paper "A solar electric propulsion mission for lunar power beaming"

"""

import matplotlib.pyplot as plt
from DVP_general_SPS_functions import *
from general_functions import convert_string_to_datetime


def scan_apogee_altitudes(apogees, start, total_duration):
    # This function analyzes the availability of the SPS to a lunar target at 45 N as a function of the
    # altitude of the orbital apogee

    total_active_time = np.zeros(len(apogees))
    total_blackout_time = np.zeros(len(apogees))

    raw_eclipse_target = 'Target1-Eclipse-Edited.csv'
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)

    for i in range(len(apogees)):
        raw_sunlight = "Lighting({}k).csv".format(apogees[i])
        raw_access = "Access({}k).csv".format(apogees[i])
        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)
        print('\n')
        print('Apogee {}000 km, Perigee 500 km'.format(apogees[i]))
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_target, access_times)
        total_blackout_time[i] = determine_blackout_data(access_times, eclipse_target, total_duration)

    return total_active_time, total_blackout_time


def scan_perigee_altitudes(perigees, start, total_duration):
    # This function analyzes the availability of the SPS to a lunar target at 45 N as a function of the
    # altitude of the orbital apogee

    total_active_time = np.zeros(len(perigees))
    total_blackout_time = np.zeros(len(perigees))

    raw_eclipse_target = 'Target1-Eclipse-Edited.csv'
    eclipse_target = parse_csv_to_array(raw_eclipse_target, start)

    for i in range(len(perigees)):
        raw_sunlight = "Lighting({},30k).csv".format(perigees[i])
        raw_access = "Access({},30k).csv".format(perigees[i])
        access_times = parse_csv_to_array(raw_access, start)
        sunlight_times = parse_csv_to_array(raw_sunlight, start)

        print('\n')
        print('Perigee {} km, Apogee 30000 km'.format(perigees[i]))
        total_active_time[i] = determine_SPS_active_time(sunlight_times, eclipse_target, access_times)
        total_blackout_time[i] = determine_blackout_data(access_times, eclipse_target, total_duration)

    return total_active_time, total_blackout_time


def get_times_without_coverage(access_times):
    length = len(access_times[0]) - 1
    no_coverage = np.zeros((3, length))
    # Check to if first access period begins at start of simulation
    # If not, no_coverage begins at start of simulation
    if access_times[0][0] != 0:
        no_coverage[0][0] = 0
        no_coverage[1][0] = access_times[0][0]
        no_coverage[2][0] = no_coverage[1][0] - no_coverage[0][0]
        for j in range(1, length):
            no_coverage[1][j] = access_times[0][j]
            no_coverage[0][j] = access_times[1][j-1]
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]
    # If access period begins at start of simulation
    elif access_times[0][0] == 0:
        for j in range(length - 2):
            no_coverage[0][j] = access_times[1][j]
            no_coverage[1][j] = access_times[0][j+1]
            no_coverage[2][j] = no_coverage[1][j] - no_coverage[0][j]

    return no_coverage


def determine_blackout_data(access_times, eclipse_target,  duration):

    # This function replicates the results of Figure 2 in Brandhorst's paper
    # regarding black-out periods for the lunar target (ie. times when there is
    # no SPS coverage, and the target is in eclipse.

    # Calculates the times when the SPS cannot access the target
    SPS_no_coverage = get_times_without_coverage(access_times)
    check_sum(access_times, SPS_no_coverage, duration)
    # Calculates overlap of no coverage with target eclipses
    dark_periods = get_event_overlaps(SPS_no_coverage, eclipse_target)
    print('Total blackout time: {} hours'.format(round(sum(dark_periods[2]) / 3600.0, 2)))
    return sum(dark_periods[2])

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

    sum_apogee = np.zeros(len(apogees))
    sum_perigee = np.zeros(len(perigees))
    for i in range(len(apogees)):
        sum_apogee[i] = (active_times_varying_apogee[i] + blackout_times_varying_apogee[i]) * 100.0 / total_duration
    for i in range(len(perigees)):
        sum_perigee[i] = (active_times_varying_perigee[i] + blackout_times_varying_perigee[i]) * 100.0 / total_duration

    plt.figure(1)
    plt.subplot(231)
    plt.plot(perigees, active_perigee_scaled)
    plt.title('Apogee Fixed 30000 km')
    plt.xlabel('Perigee Altitude [km]')
    plt.ylabel('Hours of Access (Over 2 Years)')
    plt.subplot(232)
    plt.plot(perigees, blackout_perigee_scaled)
    plt.title('Apogee Fixed 30000 km')
    plt.xlabel('Perigee Altitude [km]')
    plt.ylabel('Hours of Blackout (Over 2 Years)')
    plt.subplot(233)
    plt.plot(perigees, sum_perigee)
    plt.title('Total Active SPS/Blackout Duration')
    plt.xlabel('Perigee Altitude [km]')
    plt.ylabel('Percentage of Two Years')
    plt.subplot(234)
    plt.plot(apogees_scaled, active_apogee_scaled)
    plt.title('Perigee Fixed 500 km')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Hours of Access (Over 2 Years)')
    plt.subplot(235)
    plt.plot(apogees_scaled, blackout_apogee_scaled)
    plt.title('Perigee Fixed 500 km')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Hours of Blackout (Over 2 Years)')
    plt.subplot(236)
    plt.plot(apogees_scaled, sum_apogee)
    plt.title('Total Active SPS/Blackout Duration')
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Percentage of Two Years')
    plt.show()


main()
