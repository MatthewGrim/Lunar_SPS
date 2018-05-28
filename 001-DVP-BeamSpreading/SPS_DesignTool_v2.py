"""
Author: Darian van Paridon
Date: 03/05/2018

This script is used to generate feasible design configurations for a solar power satellite.

Inputs: Surface beam size, surface beam power, transmitter wavelength

Outputs: Required transmitter power, compatible altitudes and corresponding transmitter aperture size.

"""

import numpy as np
import matplotlib.pyplot as plt


def get_altitude_from_wavelength_surfbeamsize(wave, diameter):
    #############################################################################
    # This function determines altitudes (within the L1 point) for a solar power satellite
    # based on the available transmitter wavelength and desired surface beam diameter.

    # Inputs: Wavelength and desired surface beam diameter
    # Variables: Transmitter aperture size
    # Determine: Compatible altitudes within L1 point, and corresponding transmitter size
    #############################################################################
    wavelength = wave
    rec_radius = 0.5 * diameter

    trans_radius = np.linspace(1.0e-3, diameter, 1e3)

    # Calculate altitude and surface radius
    altitude = (np.pi * trans_radius ** 2) * np.sqrt((rec_radius / trans_radius) ** 2 - 1) / wavelength

    # Remove non-physical values, by setting altitude to zero
    L1point = 5.82e7
    flag_alt_too_high = altitude > L1point  # Where altitude is beyond L1
    flag_alt_too_low = altitude < 50e3  # Where altitude is below 50 km
    altitude[flag_alt_too_high] = None  # Set values to zero
    altitude[flag_alt_too_low] = None  # Set values to zero

    return altitude, trans_radius


def get_transpower_from_recieverarea_surfaceflux(flux, diameter):
    #############################################################################
    # This function determines the required transmitter power based on the specified
    # receiver size and surface flux requirement

    # Inputs: Desired surface flux and receiver area
    # Determine: transmitter power and aperture size
    #############################################################################

    surf_radius = 0.5 * diameter
    surf_flux = flux

    # Calculate the necessary transmitter power based on desired surface flux and surface beam diameter
    trans_power = surf_flux * np.pi * surf_radius ** 2

    return trans_power


def take_closest(array, number):
    closest = 0
    for i in range(1, len(array)):
        if abs(array[i] - number) < abs(array[closest] - number):
            closest = i
    return closest


def get_solar_array_size(trans_power, transmitter_eff):

    solar_cell_efficiency = 0.45
    power_management_efficiency = 0.95
    solar_flux = 1367
    solar_array_size = trans_power / (solar_cell_efficiency * transmitter_eff * power_management_efficiency * solar_flux)

    return solar_array_size


def main():

    # This example calculation describes a 100 kW laser (industrial welding/cutting applications) from IPG (YLS-100000).
    trans_eff = 100.0 / 290.0
    trans_wavelength = 1070e-9
    # This surface beam corresponds to the 100 kW transmitter
    surf_flux = 200
    rec_diameter = 25

    # wavelength = input('Enter the chosen transmitter wavelength (in m): ')
    # surf_flux = input('Enter the desired flux at the surface (in W/m2): ')
    altitude, trans_radius = get_altitude_from_wavelength_surfbeamsize(trans_wavelength, rec_diameter)

    plt.plot(trans_radius * 2, altitude / 1000.0)
    plt.xlabel('Transmitter Diameter [m]')
    plt.ylabel('SPS Altitude [km]')
    plt.title('Compatible Transmitter Sizes and Altitudes')
    plt.show()

    # rec_area = input('Enter the desired surface beam diameter (in m): ')
    # surf_flux = input('Enter the desired flux at the surface (in W/m2): ')
    req_trans_power = get_transpower_from_recieverarea_surfaceflux(surf_flux, rec_diameter)
    print('The minimum required transmitter power (assuming no loss) is {} kW.'.format(round(req_trans_power / 1000.0, 4)))
    solar_array_size = get_solar_array_size(req_trans_power, trans_eff)
    print('Which can be gathered with a solar array of size {} m across.'.format(round(np.sqrt(solar_array_size), 2)))
    print('\n')

    i = np.nanargmin(altitude)
    print('The minimum compatible altitude is {} km.'.format(round(altitude[i] / 1000.0, 2)))
    print('Which pairs with a transmitter diameter of {} m.'.format(round(trans_radius[i] * 2.0, 6)))
    print('\n')
    j = np.nanargmax(altitude)
    print('The maximum compatible altitude is {} km.'.format(round(altitude[j] / 1000.0, 2)))
    print('Which pairs with a transmitter diameter of {} m.'.format(round(trans_radius[j] * 2.0, 6)))


main()
