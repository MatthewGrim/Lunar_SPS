"""
Author: Darian van Paridon
Date: 01/05/2018
This file contains analysis of beam waist size for solar power satellite
to analyze configurations which shine a kW laser onto the lunar surface to power 10W rover(s)
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# Altitudes
altitude = np.arange(500000.0, 30000000.0, 50000.0)  # 500 to 30000 km altitude


# Define a range of transmitter powers to determine the magnitude of relative
# beam divergence (surface beam radius divided by transmitter radius)
# based on a requirement of 10 W/m2 surface flux

# Transmitter parameters
trans_power = np.arange(1000.0, 1000000.0, 1000.0)  # 10 to 40 kW transmitter power
laser_wavelength = 850.0*10**(-9)  # For GaAs cells
microwave_wavelength = 122.0*10**(-6)  # For 2.45 GHz microwave
trans_radius = np.array([0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10])
trans_area = np.pi*(trans_radius**2)

print(trans_area)

# Surface beam parameters
surf_power = 10.0
# surface beam waist relative to transmitter beam radius, according to conservation of flux
relative_beam_radius = np.sqrt(trans_power/surf_power)

# Satellite orbit
# For a desired surface beam size relative to transmitter size, what altitude is required?
desired_beam_div = 1000
z_laser = trans_area*np.sqrt(desired_beam_div**2 - 1)/laser_wavelength
z_microwave = trans_area*np.sqrt(desired_beam_div**2 - 1)/microwave_wavelength

plt.plot(trans_power/1000, relative_beam_radius)
plt.ylabel('Beam Radius on Surface (Normalized by Transmitter Radius)')
plt.xlabel('Transmitter Power [kW]')
plt.title('Beam Divergence for 10 W Surface Power')
plt.show()

plt.figure(1)
plt.subplot(211)
plt.plot(np.log10(trans_area), z_laser/1000, 'r--')
plt.ylabel('Required Altitude [km]')
plt.xlabel('Transmitter Area [m2]')
plt.title('Altitude Required for 1000x Beam Size Increase (Laser)')
plt.subplot(212)
plt.plot(np.log10(trans_area), z_microwave/1000, 'b--')
plt.ylabel('Required Altitude [km]')
plt.xlabel('Transmitter Area [m2]')
plt.title('Altitude Required for 1000x Beam Size Increase (Microwave)')
plt.show()




