"""
Author: Rohan Ramasamy
Date: 01/08/2018

This script contains a script to generate the 0D design plot in the conference paper being generated for the IAC. The
functions in this scripts are taken from 001-DVP-BeamSpreading
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def get_surfaceflux_from_wavelength_and_laser_power(wavelength, rover_specs, laser_powers, receiver_areas,
                                                    power_reqs, pointing_error=1e-6):
    """
    Plot the surface flux, and beam radius of a design at a specified wavelength, laser power,
    and receiver size.

    wavelength: in m
    rover_specs: name of the rover being analysed
    laser_power: in watts
    receiver_area: in m2
    power_req: in watts
    pointing_error: error in pointing in radians
    """
    assert len(power_reqs) == len(receiver_areas)
    assert len(power_reqs) == len(rover_specs)

    # Set the parameter space
    trans_radius = np.logspace(-3, 1, 1000)
    altitudes = np.logspace(4, 6.699, 1001)
    R, Z = np.meshgrid(trans_radius, altitudes, indexing="ij")

    fig, ax = plt.subplots(len(power_reqs), len(laser_powers), sharey=True, sharex=True, figsize=(15, 8))
    for i, laser_power in enumerate(laser_powers):
        for j in range(len(power_reqs)):
            rover_spec = rover_specs[j]
            receiver_area = receiver_areas[j]
            power_req = power_reqs[j]

            # Get the beam radius
            beam_radius = R * np.sqrt(1.0 + (Z * wavelength / (np.pi * R ** 2)) ** 2)
            receiver_radius = np.sqrt(receiver_area / np.pi)
            radius_constraint_one = pointing_error * Z + receiver_radius
            radius_constraint_two = pointing_error * Z + beam_radius
            mask_one = beam_radius < radius_constraint_one
            mask_two = receiver_radius > radius_constraint_two
            final_mask = np.logical_and(mask_one, np.logical_not(mask_two))
            beam_radius[final_mask] = np.nan

            # Calculate the resulting surface flux
            receiver_power = laser_power/ (np.pi * beam_radius ** 2) * receiver_area
            receiver_power[np.pi * beam_radius ** 2 < receiver_area] = laser_power
            receiver_power[receiver_power < power_req] = np.nan

            # Normalise result by input power to get total efficiency
            receiver_power /= laser_power
            receiver_power[receiver_power < 0.001] = np.nan

            log_power = np.log10(receiver_power * 100)
            ax[j, i].contourf(np.log10(R), Z / 1e3, log_power, 100)
            m = cm.ScalarMappable()
            m.set_array(log_power)
            m.set_clim(-1.0, 2.0)
            fig.colorbar(m, ax=ax[j, i])
            ax[j, 0].set_ylabel('{} \n Transmission distance [km]'.format(rover_spec))
        ax[0, i].set_title('Laser Power: {}kW'.format(laser_power / 1e3))
        ax[1, i].set_xlabel('Logarithm of Transmitter Radius [m]')
    plt.tight_layout()
    plt.show()

    return beam_radius, receiver_power


def main():
    trans_wavelength = 850e-9
    laser_power = [4e3, 15e3, 100e3]
    rover_specs = ["AMALIA in Hibernation", "AMALIA in Operation"]
    rec_area = [0.366, 0.366]
    power_req = [21.0, 300.0]
    _, _ = get_surfaceflux_from_wavelength_and_laser_power(trans_wavelength, rover_specs, laser_power, rec_area, power_req)


if __name__ == '__main__':
    main()
