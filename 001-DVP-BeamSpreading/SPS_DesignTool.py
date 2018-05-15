"""
Author: Darian van Paridon
Date: 03/05/2018
This file contains functions which are used to analyze the design space for a
solar power satellite architecture. The functions output feasible
solutions for an SPS configuration based on certain desired inputs. The relationship
between the transmitter and surface beam sizes is based on Gaussian beam divergence.
"""

import numpy as np
import matplotlib.pyplot as plt

# Relevant parameters for solar power satellite configuration:

# transmitter wavelength
# transmitter power
# transmitter beam size
# altitude
# surface beam size
# surface flux


def get_altitude_from_wavelength_surfaceflux(wave, surf):
    #############################################################################
    # This function determines altitudes (within the L1 point) for a solar power satellite
    # based on the available transmitter wavelength and desired surface flux.

    # Inputs: Wavelength and desired surface flux
    # Variables: Transmitter power and aperture size
    # Determine: Compatible altitudes within L1 point
    #############################################################################
    wavelength = wave
    surf_flux = surf

    trans_radius = np.logspace(-3, 0, 1000)
    trans_power = np.logspace(1, 7, 1000)

    # Calculate altitude and surface radius
    R, P = np.meshgrid(trans_radius, trans_power, indexing="ij")
    altitude = (np.pi * R ** 2) * np.sqrt((P / (surf_flux * np.pi * R ** 2)) - 1) / wavelength

    # Remove non-physical values, by setting altitude to zero
    L1point = 5.82e7
    flag_alt_too_high = altitude > L1point # Where altitude is beyond L1
    flag_alt_too_low = altitude < 100e3 # Where altitude is below 100 km
    altitude[flag_alt_too_high] = None  # Set values to zero
    altitude[flag_alt_too_low] = None  # Set values to zero

    # Calculate the resulting surface beam size
    surf_diameter = 2 * R * np.sqrt(1 + ((wavelength * altitude) / (np.pi * R ** 2)) ** 2)

    # This plot shows compatible altitudes for orbits about the Moon within the L1 point.
    # All non-compatible altitudes (greater than L1 distance) are set to nan.
    plt.figure(1)
    plt.subplot(211)
    plt.contourf(np.log10(trans_radius), np.log10(trans_power), altitude/1000, 100)
    cbar = plt.colorbar()
    cbar.set_label = 'Altitude [km]'
    plt.ylabel('Logarithm Transmitter Power')
    plt.title('Compatible Altitudes Within L1 [km]')
    # This plot shows the surface beam diameters associated with the solutions from
    # the previous plot (for altitude of a range of transmitter designs)
    plt.subplot(212)
    plt.contourf(np.log10(trans_radius), np.log10(trans_power), surf_diameter, 100)
    cbar = plt.colorbar()
    cbar.set_label = 'Beam Diameter [m]'
    plt.xlabel('Logarithm Transmitter Radius')
    plt.ylabel('Logarithm Transmitter Power')
    plt.title('Resulting Surface Beam Diameter [m]')
    plt.show()

    return altitude


def get_transpower_from_recieverarea_surfaceflux(flux, diameter, trans_wavelength):
    #############################################################################
    # This function determines the required transmitter power and possible aperture
    # sizes based on the specified receiver area and surface.

    # Inputs: Desired surface flux and receiver area
    # Variables: wavelength, altitude
    # Determine: transmitter power and aperture size
    #############################################################################

    surf_radius = 0.5 * diameter
    surf_flux = flux

    L1point = 5.82e7

    wavelength = np.logspace(-8, -6, 1000)
    altitude = np.linspace(1e5, L1point, 1000)

    # Calculate the necessary transmitter power based on desired surface flux and surface beam diameter
    trans_power = surf_flux * np.pi * surf_radius ** 2

    print('The minimum required transmitter power (assuming no loss) is {} kW.'.format(round(trans_power / 1000.0, 2)))

    W, A = np.meshgrid(wavelength, altitude, indexing="ij")

    # Solve quadratic equation for transmitter diameter for range of wavelengths and altitudes
    # Based on chosen surface beam size
    # Solution to quadratic equation results in "folded" solution
    trans_diameter_p = 2 * np.sqrt(0.5 * surf_radius ** 2 + np.sqrt(0.25 * surf_radius ** 4 - (W * A / np.pi) ** 2))
    trans_diameter_m = 2 * np.sqrt(0.5 * surf_radius ** 2 - np.sqrt(0.25 * surf_radius ** 4 - (W * A / np.pi) ** 2))

    i = take_closest(wavelength, trans_wavelength)
    trans_diameter_close_m = trans_diameter_m[i][:]
    trans_diameter_close_p = trans_diameter_p[i][:]

    plt.plot(altitude / 1000.0, trans_diameter_close_m)
    # plt.plot(altitude / 1000.0, trans_diameter_close_p)
    plt.show()


    # This plot shows compatible transmitter sizes for the specified
    # surface beam area and flux, as a function of wavelength and altitude
    plt.figure(1)
    plt.subplot(211)
    plt.contourf(np.log10(wavelength), altitude / 1e3, trans_diameter_p, 100)
    cbar = plt.colorbar()
    cbar.set_label = ('Altitude [km]')
    plt.ylabel('Altitude [km]')
    plt.title('Compatible Transmitter Diameter (Plus Solution) [m]')
    plt.subplot(212)
    plt.contourf(np.log10(wavelength), altitude / 1e3, trans_diameter_m, 100)
    cbar = plt.colorbar()
    cbar.set_label = ('Altitude [km]')
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Altitude [km]')
    plt.title('Compatible Transmitter Diameter (Minus Solution) [m]')
    plt.show()

    return trans_power, trans_diameter_close_m, altitude


def get_apogee_from_altitude_constraints(altitude, desired_trans_power, desired_trans_radius):

    # Copy-paste of transmitter design space from altitude function,
    # Consider taking as input to function instead of redefining
    trans_radius = np.logspace(-3, 0, 1000)
    trans_power = np.logspace(1, 7, 1000)

    i = take_closest(trans_power, desired_trans_power)
    j = take_closest(trans_radius, desired_trans_radius)
    apogee = altitude[i][j]
    print("The altitude for this transmitter is {} km.".format(round(apogee / 1000.0, 2)))

    return apogee

def take_closest(array, number):
    closest = 0
    for i in range(1, len(array)):
        if abs(array[i] - number) < abs(array[closest] - number):
            closest = i
    return closest


def main():

    trans_wavelength = 1064e-9
    surf_flux = 400
    rec_diameter = 200

    print('This first function will determine compatible altitudes for a SPS which delivers a specified surface '
          'power flux with a specified wavelength transmitter, as a function of transmitter size and power.')
    # wavelength = input('Enter the chosen transmitter wavelength (in m): ')
    # surf_flux = input('Enter the desired flux at the surface (in W/m2): ')
    altitude = get_altitude_from_wavelength_surfaceflux(trans_wavelength, surf_flux)

    print('This second function will determine minimum necessary transmitter power and size for a '
         'desired receiver area and the specified flux at the target, as a function of wavelength and altitude.')
    # rec_area = input('Enter the desired surface beam diameter (in m): ')
    # surf_flux = input('Enter the desired flux at the surface (in W/m2): ')
    req_trans_power, trans_dia_m, altitude_fortransdia = get_transpower_from_recieverarea_surfaceflux(surf_flux, rec_diameter, trans_wavelength)

    print('This function takes in the altitudes compatible with the desired radiation wavelength and surface flux, '
          'then finds the upper altitude limits for a given transmitter power and size.')
    # trans_power = input('Enter maximum transmitter power (in W): ')
    # trans_diameter = input('Enter diameter of transmitter (in m): ')
    apogee = get_apogee_from_altitude_constraints(altitude, req_trans_power, 0.1)

    i = take_closest(altitude_fortransdia, apogee)
    print('The compatible transmitter diameter is {} m'.format(round(trans_dia_m[i], 2)))


main()



