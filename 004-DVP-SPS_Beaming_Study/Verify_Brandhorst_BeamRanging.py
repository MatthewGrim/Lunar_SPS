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

    return (surf_flux / solar_flux), percent_rec_covered, power_delivered


def main():

    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
    raw_range_data = 'SPS-to-Target-Range.csv'
    range_data = import_range_data(raw_range_data, start)
    surf_flux, percent_rec_covered, power_delivered = determine_surface_beam_data(range_data)

    range_days = range_data[0] / 86400.0

    # Brandhorst figures begin on the 414th (7:55 am) day of the simulation, and end on the 416th day (7:55 am)
    brandhorst_range = range_days[(range_days >= 414.0) & (range_days <= 416.0)]

    print(range_days)
    print(brandhorst_range)
    plt.figure(1)
    plt.subplot(221)
    plt.scatter(range_days, range_data[1], s=0.1)
    plt.ylabel('Distance to Target [km]')
    plt.subplot(222)
    plt.scatter(range_days, surf_flux,s=0.1)
    plt.ylabel('Surface Beam Flux [AM0]')
    plt.subplot(223)
    plt.scatter(range_days, percent_rec_covered,s=0.1)
    plt.xlabel('Days Since Simulation Start')
    plt.ylabel('Percentage of Receiver Covered')
    plt.subplot(224)
    plt.scatter(range_days, power_delivered / 1000.0, s=0.1)
    plt.xlabel('Days Since Simulation Start')
    plt.ylabel('Power Delivered at Receiver [kW]')
    plt.show()


main()
