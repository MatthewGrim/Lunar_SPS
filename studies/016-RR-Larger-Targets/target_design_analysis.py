"""
Author: Rohan Ramasamy
Date: 05/12/2018

This script contains code to define possible SPS target designs
"""

import numpy as np
import matplotlib.pyplot as plt

from Lunar_SPS.pysrc.utils.unit_conversions import UnitConversions
from Lunar_SPS.pysrc.utils.physical_constants import PhysicalConstants


def altitudes_for_pointing(eta_link, solar_efficiency=0.2, laser_efficiency=0.5, P_hib_ratio=0.0,
	                       pointing_requirement=1e-7, laser_requirement=4e3):
	a = np.linspace(8e5, 2e6, 100)
	r_moon = 1737.0e3
	r = a + r_moon

	P_rec = np.linspace(20.0, 300.0, 200)
	
	# Ger receiver size and hibernation power
	A, P_REC = np.meshgrid(a, P_rec, indexing='ij')
	R = A + r_moon
	Z = np.sqrt(r_moon ** 2 + R ** 2 - np.sqrt(2) * r_moon * R)
	Z_MAX = np.sqrt(R ** 2 - r_moon ** 2)
	Z = Z_MAX
	W_R = P_REC / (solar_efficiency * PhysicalConstants.solar_irradiance)
	P_REC_HIB = P_rec * P_hib_ratio

	# Get approximate percentage active time
	percentage_active_best_case = np.arccos(np.sqrt(2) * r_moon / R) / np.pi

	# Get required battery capacity
	mu_moon = PhysicalConstants.mu_moon
	v_orbit = np.sqrt(mu_moon / R)
	omega_orbit = v_orbit / R
	T_ORBIT = 2 * np.pi / omega_orbit
	T_hibernation = P_hib_ratio * 0.5 * T_ORBIT
	E_CAPACITY = T_hibernation * P_REC_HIB * UnitConversions.J_to_Whr

	# Get required pointing
	eta_root = np.sqrt(eta_link)
	SIGMA_POINTING = W_R * (1 - eta_root) / (Z * eta_root) * UnitConversions.radian_to_microradian

	# Get laser power
	P_LAS = P_REC / (laser_efficiency * eta_link) * UnitConversions.W_to_kW
	P_LAS /= (1 - np.exp(-2))

	# Apply constraints - laser exists, pointing is manageable, power is low, active time is achievable at altitude
	active_time_mask = percentage_active_best_case < P_hib_ratio
	wavelength = 1070e-9
	W_B = np.sqrt(eta_link * W_R ** 2)
	transmitter_mask = W_B ** 4 > 4 * (Z * wavelength / np.pi) ** 2
	pointing_mask = SIGMA_POINTING < pointing_requirement * UnitConversions.radian_to_microradian
	power_mask = P_LAS > laser_requirement * UnitConversions.W_to_kW
	for res_array in [W_R, SIGMA_POINTING, P_LAS, E_CAPACITY]:
		res_array[pointing_mask] = 0.0
		res_array[power_mask] = 0.0
		res_array[transmitter_mask] = 0.0
		res_array[active_time_mask] = 0.0

	fig, ax = plt.subplots(4, figsize=(12, 7), sharex=True)

	num_contours = 100
	A *= UnitConversions.m_to_km
	im = ax[0].contourf(A, P_REC, W_R, num_contours)
	ax[0].set_title("Receiver radius [$m$]")
	ax[0].set_ylabel("Power required [W]")
	
	fig.colorbar(im, ax=ax[0])
	im = ax[1].contourf(A, P_REC, SIGMA_POINTING, num_contours)
	ax[1].set_title("Pointing accuracy [$\mu Rad$]")
	ax[1].set_ylabel("Power required [W]")
	fig.colorbar(im, ax=ax[1])	

	im = ax[2].contourf(A, P_REC, P_LAS, num_contours)
	ax[2].set_title("Laser Power [kW]")
	ax[2].set_ylabel("Power required [W]")
	fig.colorbar(im, ax=ax[2])

	im = ax[3].contourf(A, P_REC, E_CAPACITY, num_contours)
	ax[3].set_title("Rover Battery Capacity [Whr]")
	ax[3].set_ylabel("Power required [W]")
	ax[3].set_xlabel("Altitude [km]")
	fig.colorbar(im, ax=ax[3])

	plt.show()


if __name__ == '__main__':
	eta_link = 0.08

	altitudes_for_pointing(eta_link)

	