""""
01/06/2018
Author: Darian van Paridon

This script calculates the required battery mass for an SPS system as a function of the transmitter power
and maximum beaming duration (while eclipsed)

"""

import numpy as np
import matplotlib.pyplot as plt

trans_power = np.linspace(1e3, 100e3, 1000)  # transmitter power in watts
duration = np.linspace(0.5, 4.0, 1001)  # beaming duration in hours
P, T = np.meshgrid(trans_power, duration, indexing="ij")


battery_capacity = P * T  # battery capacity in watt hours
li_ion_specific_power = 270  # specific power in Wh/kg
battery_mass = battery_capacity / li_ion_specific_power  # estimate battery mass

plt.contourf(P / 1000.0, T, battery_mass, 100)
plt.xlabel('Transmitter Power [kW]')
plt.ylabel('Eclipsed Beaming Duration [h]')
plt.colorbar()
plt.show()
