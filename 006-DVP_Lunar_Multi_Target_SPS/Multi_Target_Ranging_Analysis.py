""""
24/05/2018
Author: Darian van Paridon

This script is meant to analyze the beaming capabilities for a single SPS with
multiple targets which are distributed at various latitudes relative to the orbital plane of the SPS.

The data currently associated with this script are for an elliptic (300 km perigee to 15000 km apogee altitude) and a
circular orbit (3000 km altitude) in the equatorial plane of the moon. The apogee of the elliptic orbits is directly
overhead the longitude of the targets.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
import os


def eliminate_event_periods(range_data, event_data):

    # Remove range data points which are time stamped during an
    # conditional event in which the SPS would not be active
    for i in range(0, len(range_data[0])):
        for j in range(0, len(event_data[0])):
            if event_data[0][j] <= range_data[0][i] <= event_data[1][j]:
                range_data[1][i] = 'nan'

    return range_data


def scan_latitudes(start, total_duration, scan_array, orbit):
    # This function is used to scan through a series of .csv files from STK to analyze
    # the coverage availability for various targets

    # For accessing the file in a folder contained in the current folder
    fileDir = os.path.dirname(os.path.realpath('__file__'))

    # Import SPS illumination events
    sps_lighting_raw = os.path.join(fileDir, '{}/SPS_Lighting.csv'.format(orbit))
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    sps_eclipse = invert_events_list(sps_lighting, total_duration)

    times = []
    ranges = []
    blackout_durations = []
    surf_flux = []

    # Cycle through available target positions, as function of latitude
    for i in scan_array:
        print("For target at {} degree latitude".format(round(i / 10.0, 2)))

        # Import target lighting data, and invert to find eclipse events
        raw_target_lighting = os.path.join(fileDir, "{}/Target_Lighting_{}.csv".format(orbit, i))
        target_lighting = parse_csv_to_array(raw_target_lighting, start)
        target_eclipse = invert_events_list(target_lighting, total_duration)

        # Import range data for SPS (time, range to target)
        raw_range = os.path.join(fileDir, "{}/SPS_Range_{}.csv".format(orbit, i))
        sps_range = import_range_data(raw_range, start)

        # Import access data for sps
        raw_access = os.path.join(fileDir, "{}/SPS_Access_{}.csv".format(orbit, i))
        sps_access = parse_csv_to_array(raw_access, start)

        # Remove range samples which occur during inactive SPS periods
        sps_range_wo_eclipse = eliminate_event_periods(sps_range, sps_eclipse)
        corrected_range_data = eliminate_event_periods(sps_range_wo_eclipse, target_lighting)

        # Calculate surface beam intensity
        surface_intensity = determine_surface_beam_data(corrected_range_data)
        surf_flux.append(surface_intensity)

        # Collect remaining time and range samples
        times.append(corrected_range_data[0])
        ranges.append(corrected_range_data[1])

        # Determine blackout periods
        sps_active = determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)
        blackout_data = determine_blackout_data(sps_active, target_eclipse, total_duration)

        # Collect blackout durations
        blackout_durations.append(blackout_data[2])

    return times, ranges, blackout_durations, surf_flux


def determine_surface_beam_data(range_data):

    # Transmitter parameters used in Brandhorst study
    wavelength = 1070e-9
    trans_diameter = 0.2
    trans_power = 100e3

    # Gaussian beam divergence to determine size of beam at the target
    trans_radius = 0.5 * trans_diameter
    surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (range_data[1] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)

    # Determine surface flux relative to solar flux for comparison to Brandhorst
    surf_flux = (trans_power / (np.pi * surf_beam_radius ** 2))

    return surf_flux


def main():

    # Initialize simulation times
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2019', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Select which orbit to use. Options are Circle and Elliptic
    orbit = 'Elliptic'

    # Latitude of list of targets (multiplied by ten for processing the file names which don't contain decimals)
    if orbit == 'Elliptic':
        latitudes = np.asarray([0, 125, 250, 375, 500])
    elif orbit == 'Circle':
        latitudes = np.asarray([0, 125, 250, 375, 500, 625])
    times, ranges, blackouts, surf_flux = scan_latitudes(start, total_duration, latitudes, orbit)

    max_blackout = []
    for i in range(len(latitudes)):
        max_blackout.append(max(blackouts[i]))

    plt.bar([i / 10.0 for i in latitudes], [i / 3600.0 for i in max_blackout])
    plt.xlabel('Degrees Latitude')
    plt.ylabel('Duration [h]')
    plt.title('Maximum Blackout Duration vs Target Latitude')
    plt.show()

    for i in range(len(times)):
        plt.subplot(211)
        plt.plot(times[i] / 86400.0, ranges[i], label='{} Degree Latitude'.format(latitudes[i]/10.0))
        plt.xlabel("Days")
        plt.ylabel('Range [km]')
        plt.title('Elliptical SPS, Multiple Targets')
        plt.legend()
        plt.subplot(212)
        plt.plot(times[i] / 86400.0, surf_flux[i], label='{} Degree Latitude'.format(latitudes[i]/10.0))
        plt.xlabel("Days")
        plt.ylabel('Surface Beam Flux [W/m2]')
        plt.title('Surface Beam Flux (100 kW transmitter)')
        plt.legend()
    plt.show()


main()
