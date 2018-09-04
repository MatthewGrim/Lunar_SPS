""""
Date: 14/05/2018
Author: Darian van Paridon

This file is for verifying the SPS wireless power transmission link results from Brandhorst (shown in Figures 8 and 9
in the lunar SPS paper " A solar electric propulsion mission for lunar power beaming "). The assumption is that the
SPS which the beam sizing study in Brandhorst's paper is the SPS with an argument of perigee of 0 degrees.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
import matplotlib.pylab as plt


def determine_surface_beam_data(sps_range):

    # Transmitter parameters used in Brandhorst study
    wavelength = 850e-9
    trans_area = 4.0
    trans_power = 40e3

    # Receiver parameters used in Brandhorst study
    pwr_area_metric = 300.0
    total_rec_pwr = 60e3
    rec_area = total_rec_pwr / pwr_area_metric

    # Efficiency assumed in Brandhorst study
    solar_cell_eff = 0.45

    # General AM0 solar irradiance
    solar_flux = 1367.0

    # Gaussian beam divergence to determine size of beam at the target
    trans_radius = np.sqrt(trans_area / np.pi)
    surf_beam_radius = trans_radius * np.sqrt(1 + (wavelength * (sps_range[1] * 1000.0) / (np.pi * trans_radius ** 2)) ** 2)
    surf_beam_area = np.pi * surf_beam_radius ** 2

    # Determine surface flux relative to solar flux for comparison to Brandhorst
    surf_flux = (trans_power / (np.pi * surf_beam_radius ** 2))
    # Determine percentage of the total receiver area which is covered for comparison to Brandhorst
    percent_rec_covered = (surf_beam_area / rec_area) * 100.0
    # Determine total power delivered to receiver
    power_delivered = surf_beam_area * surf_flux * solar_cell_eff

    return (surf_flux / solar_flux), percent_rec_covered, power_delivered


def eliminate_event_periods(sps_range, event_data):

    # Remove range data points which are time stamped during an
    # conditional event in which the SPS would not be active
    for i in range(0, len(sps_range[0])):
        for j in range(0, len(event_data[0])):
            if event_data[0][j] <= sps_range[0][i] <= event_data[1][j]:
                sps_range[1][i] = 'nan'

    return sps_range


def main():

    # Initialize simulation time
    start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])

    # Import sps to target range data
    raw_sps_range = 'SPS-to-Target-Range.csv'
    sps_range = import_range_data(raw_sps_range, start)

    # Import lighting data
    raw_sps_eclipse = 'SPS1-Eclipse(0)-Edited.csv'
    raw_target_lighting = 'Target1-Lighting-Edited.csv'
    sps_eclipse = parse_csv_to_array(raw_sps_eclipse, start)
    target_lighting = parse_csv_to_array(raw_target_lighting, start)

    # Remove sps range data points which are outside active events
    sps_range_wo_eclipses = eliminate_event_periods(sps_range, sps_eclipse)
    sps_range_wo_eclipse_and_target_sunlight = eliminate_event_periods(sps_range_wo_eclipses, target_lighting)

    surf_flux, percent_rec_covered, power_delivered = determine_surface_beam_data(sps_range_wo_eclipse_and_target_sunlight)

    # Brandhorst figures begin on the 414th (7:55 am) day of the simulation, and end on the 416th day (7:55 am)

    plt.figure(1)
    plt.subplot(311)
    plt.scatter(sps_range[0] / 86400.0, sps_range[1], s=1.0)
    plt.title('SPS Beaming Range and Surface Flux')
    plt.ylabel('Distance to Target [km]')
    plt.xlim([206, 217])
    plt.subplot(312)
    plt.scatter(sps_range[0] / 86400.0, surf_flux, s=1.0)
    plt.ylabel('Surface Beam Flux [AM0]')
    plt.xlim([206, 217])
    plt.xlabel('Days')
    plt.subplot(313)
    plt.scatter(sps_range[0] / 86400.0, percent_rec_covered, s=1.0)
    plt.xlabel('Days Since Simulation Start')
    plt.ylabel('Percentage of Receiver Covered')
    plt.xlim([206, 217])
    plt.show()


main()
