"""
Author: Rohan Ramasamy
Date: 22/11/2018

This script contains a script to generate the 0D design plots that show the minimum power necessary to enable a link
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def get_surfaceflux_from_wavelength_and_laser_power(wavelength, rover_specs, max_laser_power, receiver_areas,
                                                    power_reqs, pointing_error=[1e-6, 1e-7]):
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
    altitudes = np.logspace(4, 6.67, 1001)
    R, Z = np.meshgrid(trans_radius, altitudes, indexing="ij")

    fig, ax = plt.subplots(len(power_reqs), sharey=True, sharex=True, figsize=(12, 7))
    for j in range(len(power_reqs)):
        rover_spec = rover_specs[j]
        receiver_area = receiver_areas[j]
        power_req = power_reqs[j]

        # Get the beam radius
        beam_radius = R * np.sqrt(1.0 + (Z * wavelength / (np.pi * R ** 2)) ** 2)
        receiver_radius = np.sqrt(receiver_area / np.pi)
        radius_constraint_one = pointing_error[j] * Z + receiver_radius
        radius_constraint_two = pointing_error[j] * Z + beam_radius
        mask_one = beam_radius < radius_constraint_one
        mask_two = receiver_radius > radius_constraint_two
        final_mask = np.logical_and(mask_one, np.logical_not(mask_two))
        final_mask = mask_one
        beam_radius[final_mask] = np.nan

        # Calculate the resulting surface flux
        transmitter_power = power_reqs[j] * (np.pi * beam_radius ** 2) / receiver_area
        transmitter_power[np.pi * beam_radius ** 2 < receiver_area] = power_req
        
        # Remove powers above threshold
        transmitter_power[transmitter_power > max_laser_power] = np.nan

        im = ax[j].contourf(np.log10(R), Z * 1e-3, transmitter_power * 1e-3, 100)
        fig.colorbar(im, ax=ax[j])
        ax[j].set_ylabel('{} \n Transmission distance [km]'.format(rover_spec))
    ax[1].set_xlabel('Logarithm of Transmitter Radius [m]')
    plt.tight_layout()
    plt.show()

    return beam_radius, transmitter_power


def main():
    max_laser_power = 10e3
    trans_wavelength = 850e-9
    # rover_specs = ["ISRU Excavator $\sigma_p = 0.1 \mu Rad$", "ISRU Plant $\sigma_p = 0.1 \mu Rad$"]
    # rec_area = [2.779809802487198, 19.568397951719092]
    # power_req = [1520, 10700]
    # pointing_error = [1e-7, 1e-7]
    # rover_specs = ["PHILIP Rover $\sigma_p = 1 \mu Rad$", "PHILIP Rover $\sigma_p = 0.1 \mu Rad$"]
    # rec_area = [0.7088, 0.7088]
    # power_req = [140.0, 140.0]
    # pointing_error = [1e-6, 1e-7]
    # rover_specs = ["Sorato $\sigma_p = 1 \mu Rad$", "Sorato $\sigma_p = 0.1 \mu Rad$"]
    # rec_area = [0.079, 0.079]
    # power_req = [43.0, 43.0]
    # pointing_error = [1e-6, 1e-7]
    rover_specs = ["AMALIA $\sigma_p = 1 \mu Rad$", "AMALIA $\sigma_p = 0.1 \mu Rad$"]
    rec_area = [0.366, 0.366]
    power_req = [200.0, 200.0]
    pointing_error = [1e-6, 1e-7]
    _, _ = get_surfaceflux_from_wavelength_and_laser_power(trans_wavelength, rover_specs, max_laser_power,rec_area, power_req, pointing_error=pointing_error)

if __name__ == '__main__':
    main()

