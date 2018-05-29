""""
Date: 14/05/2018
Author: Darian van Paridon

This file is for verifying the results of Brandhorst (shown in Figures 8 and 9 in his lunar SPS paper ___).

"""

from DVP_general_SPS_functions import *
import matplotlib.pylab as plt


def import_range_data(file_name, sim_start):
    the_file = open(file_name, "r")
    size = file_len(file_name)
    sps_range = np.zeros(size - 1)
    start_time_sec_from_simstart = np.zeros(size - 1)
    parsed_data = np.array((3, size))

    # For loop parses data into three categories: start, end and duration of event
    # This function can be used for target to satellite access times, or target/satellite illumination events
    # .csv file must be three columns, with start time, end time, duration in that order
    # The function outputs the start and end times of the event in seconds since the beginning of the simulation
    for i, line in enumerate(the_file):
        if i == 0:
            continue
        # Split line into table components, assuming comma delimited
        components = line.split(",")

        # Work out sunlit times and set new time t
        start = components[0].split(":")
        start_time = convert_string_to_datetime(start)
        start_time_sec_from_simstart[i - 1] = (start_time - sim_start).total_seconds()

        sps_range[i - 1] = components[1]

        parsed_data = [start_time_sec_from_simstart, sps_range]

    return parsed_data


def determine_surface_beam_data(range_data):

    # Transmitter parameters used in Brandhorst study
    wavelength = 850e-9
    trans_area = 4
    trans_power = 40e3

    # Receiver parameters used in Brandhorst study
    pwr_area_metric = 300
    total_rec_pwr = 60e3
    rec_area = total_rec_pwr / pwr_area_metric
    # Efficiency assumed in Brandhorst study
    solar_cell_eff = 0.45

    # General AM0 solar irradiance
    solar_flux = 1367

    # Gaussian beam divergence to determine size of beam at the target
    trans_radius = np.sqrt(trans_area / np.pi)
    surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (range_data[1] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
    surf_beam_area = np.pi * surf_beam_radius ** 2

    # Determine surface flux relative to solar flux for comparison to Brandhorst
    surf_flux = (trans_power / (np.pi * surf_beam_radius ** 2))
    # Determine percentage of the total receiver area which is covered for comparison to Brandhorst
    percent_rec_covered = (surf_beam_area / rec_area) * 100
    # Determine total power delivered to receiver
    power_delivered = surf_beam_area * surf_flux * solar_cell_eff

    print(max(percent_rec_covered) * 2.0)
    print(min(percent_rec_covered) * 2.0)

    return (surf_flux / solar_flux), percent_rec_covered, power_delivered


def eliminate_event_periods(range_data, event_data):

    # Remove range data points which are time stamped during an
    # conditional event in which the SPS would not be active
    for i in range(0, len(range_data[0])):
        for j in range(0, len(event_data[0])):
            if event_data[0][j] <= range_data[0][i] <= event_data[1][j]:
                range_data[1][i] = 'nan'

    return range_data


def main():

    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    raw_range_data = 'SPS-to-Target-Range.csv'
    range_data = import_range_data(raw_range_data, start)

    raw_sps_eclipse = 'SPS1-Eclipse(0)-Edited.csv'
    raw_target_sunlit = 'Target1-Lighting-Edited.csv'
    sps_eclipse = parse_csv_to_array(raw_sps_eclipse, start)
    target_sunlit = parse_csv_to_array(raw_target_sunlit, start)

    range_data_wo_eclipses = eliminate_event_periods(range_data, sps_eclipse)
    range_data_wo_eclipse_and_target_sunlight = eliminate_event_periods(range_data_wo_eclipses, target_sunlit)

    surf_flux, percent_rec_covered, power_delivered = determine_surface_beam_data(range_data_wo_eclipse_and_target_sunlight)

    range_days = range_data[0] / 86400.0

    # Brandhorst figures begin on the 414th (7:55 am) day of the simulation, and end on the 416th day (7:55 am)

    plt.figure(1)
    plt.subplot(211)
    plt.scatter(range_days, range_data[1], s=1.0)
    plt.title('SPS Beaming Range and Surface Flux')
    plt.ylabel('Distance to Target [km]')
    plt.xlim([206, 217])
    plt.subplot(212)
    plt.scatter(range_days, surf_flux, s=1.0)
    plt.ylabel('Surface Beam Flux [AM0]')
    plt.xlim([206, 217])
    plt.xlabel('Days')

    plt.figure(2)
    plt.subplot(211)
    plt.scatter(range_days, percent_rec_covered, s=1.0)
    plt.xlabel('Days Since Simulation Start')
    plt.ylabel('Percentage of Receiver Covered')
    plt.xlim([206, 217])
    plt.subplot(212)
    plt.scatter(range_days, power_delivered / 1000.0, s=1.0)
    plt.xlabel('Days Since Simulation Start')
    plt.ylabel('Power Delivered at Receiver [kW]')
    plt.xlim([206, 217])
    plt.show()


main()
