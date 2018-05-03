"""
Author: Darian van Paridon
Date: 02/05/2018
This file contains analysis of beam waist size for solar power satellite
to analyze configurations which shine a kW laser onto the lunar surface to power 10W rover(s)
PART 1 concerns size of beam waist for converting small (mm - cm) kW transmitter beam into 10 W/m2 surface beam
PART 2 concerns separation of small (mm-cm) transmitter for large surface beam (km)
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# PART 1
# Define a range of transmitter powers and determine the magnitude of relative
# beam divergence (surface beam radius divided by transmitter radius)
# based on a requirement of 10 W/m2 surface flux

# Altitudes
altitude = np.arange(5e5, 3e7, 5e4)  # 500 to 30000 km altitude

# Transmitter parameters
trans_power = np.arange(1e3, 1e5, 1e2)  # 10 to 100 kW transmitter power
laser_wavelength = 850e-9  # For GaAs cells
microwave_wavelength = 122e-6  # For 2.45 GHz microwave
trans_radius = np.arange(1e-3, 1e-1, 1e-4)
trans_area = np.pi * trans_radius ** 2

# Surface beam parameters
surf_power = 10.0
# surface beam radius relative to transmitter beam radius, according to conservation of flux
relative_beam_radius = np.sqrt(trans_power / surf_power)
relative_beam_area = trans_power / surf_power

plt.figure(1)
plt.subplot(211)
plt.plot(trans_power / 1000.0, relative_beam_radius)
plt.ylabel('Rs / Rt')
plt.title('Beam Divergence for 10 W/m2 Surface Flux')
plt.subplot(212)
plt.plot(trans_power / 1000.0, relative_beam_area)
plt.ylabel('As / At')
plt.xlabel('Transmitter Power [kW]')
plt.show()
plt.show()

# PART 2
# Specify a range of transmitter sizes (mm - cm range), and determine the distance to the surface target which
# is required for a increase in beam radius into the km range, for laser and microwave radiation

# Satellite orbit
# For a desired surface beam size relative to transmitter size, what altitude is required?
# Chosen beam radius ratio is 1e6 --> 1 mm goes to 1 km, 10 cm goes to 100 km
desired_beam_div = 1e6
z_laser = trans_area * np.sqrt(desired_beam_div ** 2 - 1) / laser_wavelength
z_microwave = trans_area * np.sqrt(desired_beam_div ** 2 - 1) / microwave_wavelength

# What is the resulting surface beam flux, based on a 10 MW transmitter?
# 10 kW transmitter results in a flux << 1W/m^2
surf_radius = trans_radius * desired_beam_div
surf_flux = 1e8 / (np.pi * surf_radius ** 2)

plt.figure(2)
laser, = plt.plot(np.log10(trans_radius), z_laser / 1000.0, 'r--', label='Laser 850 nm')
microwave, = plt.plot(np.log10(trans_radius), z_microwave / 1000.0, 'b--', label='Microwave 122 um')
plt.ylabel('Required Altitude [km]')
plt.xlabel('Logarithm of Transmitter Radius')
plt.title('Altitude Required for Surface Beam Radius 1000000X Larger Than Transmitter')
plt.legend(handles=[laser, microwave])
# set upper limit at distance from Earth to Moon: 384,400 km
plt.ylim(-100, 384400)
plt.xlim(-3, -1)
plt.show()

plt.figure(3)
flux, = plt.plot(np.log10(trans_radius), surf_flux, label='Surface Flux [W/m2]')
plt.ylabel('Surface Flux [W/m^2]')
plt.xlabel('Logarithm of Transmitter Radius')
plt.title('Surface Flux for 1000000X Increase in Beam Radius, for 100 MW laser')
plt.xlim(-3, -1)
plt.show()



