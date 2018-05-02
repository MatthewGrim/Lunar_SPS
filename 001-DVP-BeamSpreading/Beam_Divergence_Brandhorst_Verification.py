"""
Author: Darian van Paridon
Date: 01/05/2018
This file contains analysis of beam waist size for solar power satellite
for verification of the laser beaming results in the Brandhorst paper
"""

import numpy as np

# Altitudes
perigee = 5e5
apogee = 3e7
altitude = np.asarray([500000.0, 30000000.0])

# Transmitter parameters, based on Brandhorst
trans_area = 4.0
trans_power = 40000.0
wavelength = 850e-9
trans_beam_radius = np.sqrt(trans_area/np.pi)

# Surface beam parameters, based on Gaussian beam divergence
surf_beam_radius = trans_beam_radius * np.sqrt(1 + (altitude * wavelength / trans_area)**2)
surf_beam_area = np.pi * (surf_beam_radius ** 2)

print(surf_beam_area)

# Brandhorst claims that at maximum surface beam size (assumed to occur at apogee)
# the beam covers only 60% of the receiver area
rec_area = (1/0.6)*surf_beam_area

print(rec_area.item(1))

# The flux at the surface of the moon
surf_flux = np.divide(trans_power, surf_beam_area)

# AM0 solar irradiance
sun = 1377  # W/m2

# Comparing laser beam flux to solar irradiance to verify "laser beaming results" in Brandhorst
relative_flux = np.divide(surf_flux, sun)

print(relative_flux)




