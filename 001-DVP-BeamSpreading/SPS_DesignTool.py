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
    altitude = (np.pi * R ** 2) * np.sqrt(P / (np.pi * surf_flux * R ** 2) - 1) / wavelength

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


def get_transpower_from_recieverarea_surfaceflux(flux, diameter):
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

    wavelength = np.logspace(-7, -4, 1000)
    altitude = np.linspace(100e3, L1point, 1000)

    # Calculate the necessary transmitter power based on desired surface flux and surface beam diameter
    trans_power = surf_flux * np.pi * surf_radius ** 2

    print('Required transmitter power (assuming no loss):')
    print(trans_power)

    W, A = np.meshgrid(wavelength, altitude, indexing="ij")

    # Solve quadratic equation for transmitter diameter for range of wavelengths and altitudes
    # Based on chosen surface beam size
    # Solution to quadratic equation results in "folded" solution
    trans_diameter_p = 2 * np.sqrt(0.5 * surf_radius ** 2 + np.sqrt(0.25 * surf_radius ** 4 - (W * A / np.pi) ** 2))
    trans_diameter_m = 2 * np.sqrt(0.5 * surf_radius ** 2 - np.sqrt(0.25 * surf_radius ** 4 - (W * A / np.pi) ** 2))

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


def main():
    print('This first function will determine possible altitudes for a SPS with a specified wavelength transmitter '
          'and surface flux quantity, as a function of transmitter size and power.')
    wavelength = input('Enter the chosen transmitter wavelength (in m): ')
    surf_flux = input('Enter the desired flux at the surface (in W/m2): ')
    get_altitude_from_wavelength_surfaceflux(wavelength, surf_flux)

    print('This second function will determine available transmitter power and size configurations for a '
         'specific receiver area and flux at the target, as a function of wavelength and altitude.')
    rec_area = input('Enter the desired surface beam diameter (in m): ')
    surf_flux = input('Enter the desired flux at the surface (in W/m2): ')
    get_transpower_from_recieverarea_surfaceflux(rec_area, surf_flux)


main()



