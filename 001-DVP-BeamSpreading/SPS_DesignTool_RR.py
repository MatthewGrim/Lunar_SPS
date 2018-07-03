"""
Author: Rohan Ramasamy
Date: 28/05/2018

This script carries out some simple analysis for the laser architecture required for a particular
altitude, subject to received power and pointing error constraints.  
"""

import numpy as np
import matplotlib.pyplot as plt


def get_surfaceflux_from_wavelength_and_laser_power(wavelength, laser_powers, receiver_area, 
                                                    power_req, pointing_error=1e-6):
    """
    Plot the surface flux, and beam radius of a design at a specified wavelength, laser power, 
    and receiver size.

    wavelength: in m
    laser_power: in watts
    receiver_area: in m2
    power_req: in watts
    pointing_error: error in pointing in radians
    """
    # Set the parameter space
    trans_radius = np.logspace(-3, 1, 1000)
    altitudes = np.logspace(4, 6.699, 1001)
    R, Z = np.meshgrid(trans_radius, altitudes, indexing="ij")

    fig, ax = plt.subplots(1, len(laser_powers), sharey=True)
    for i, laser_power in enumerate(laser_powers):
        # Get the beam radius
        beam_radius = R * np.sqrt(1.0 + (Z * wavelength / (np.pi * R ** 2)) ** 2)
        receiver_radius = np.sqrt(receiver_area / np.pi)
        radius_constraint = pointing_error * Z + receiver_radius
        beam_radius[beam_radius < radius_constraint] = np.nan

        # Calculate the resulting surface flux
        receiver_power = laser_power / (np.pi * beam_radius ** 2) * receiver_area
        receiver_power[np.pi * beam_radius ** 2  < receiver_area] = laser_power
        receiver_power[receiver_power < power_req] = np.nan

        # Normalise result by input power to get total efficiency
        receiver_power /= laser_power
        receiver_power[receiver_power < 0.01] = np.nan

        im = ax[i].contourf(np.log10(R), Z / 1e3, receiver_power, 100)
        fig.colorbar(im, ax=ax[i])
        ax[i].set_title('Laser Power: {}kW'.format(laser_power / 1e3))
        ax[i].set_xlabel('Logarithm Radius [m]')
    ax[0].set_ylabel('Altitude [km]')
    fig.suptitle("Total efficiency for possible architectures to achieve greater than {}W receiver power using a {} wavelength laser.".format(power_req, wavelength))
    plt.show()

    return beam_radius, receiver_power


def main():
    trans_wavelength = 1070e-9
    laser_power = [4e3, 15e3, 100e3]
    rec_area = 1.0
    power_req = 270.0
    _, _ = get_surfaceflux_from_wavelength_and_laser_power(trans_wavelength, laser_power, rec_area, power_req)


if __name__ == '__main__':
    main()

