"""
Author: Rohan Ramasamy
Date: 05/12/2018

This script contains code to define possible SPS target designs
"""

import numpy as np
import matplotlib.pyplot as plt

from Lunar_SPS.pysrc.utils.unit_conversions import UnitConversions
from Lunar_SPS.pysrc.utils.physical_constants import PhysicalConstants


def altitudes_for_pointing(pointing_accuracy, solar_efficiency=0.1, laser_efficiency=0.5, P_hib_ratio=0.1,
	                       panel_heating=12.31, eta_link_requirement=0.01, laser_requirement=7.5e3, DoD=0.3,
	                       a=None):
	"""
	pointing_accuracy: Assumed pointing accuracy of beam
	solar_efficiency: day time efficiency of rover receiver cells
	laser_efficiency: night time efficiency of rover receiver cells
	P_hib_radio: assumed ratio of hibernation and charging power
	panel_heating: approxomimate heating required per m2
	eta_link_requirement: minimum acceptable link efficiency
	laser_requirement: maximum acceptable laser requirement
	DoD: depth of discharge for battery systems
	"""
	if a is None:
		a = np.linspace(8e5, 2e6, 500)
	r_moon = 1737.0e3
	r = a + r_moon

	P_op = np.linspace(20.0, 200.0, 1000)
	
	# Get receiver size and hibernation power
	if isinstance(a, np.ndarray):
		A, P_OP = np.meshgrid(a, P_op, indexing='ij')
	else:
		A = a
		P_OP = P_op

	R = A + r_moon
	Z = np.sqrt(r_moon ** 2 + R ** 2 - np.sqrt(2) * r_moon * R)
	Z_MAX = np.sqrt(R ** 2 - r_moon ** 2)
	Z = Z_MAX
	W_R = np.sqrt(P_OP / (np.pi * solar_efficiency * PhysicalConstants.solar_irradiance))
	A_R = np.pi * W_R ** 2
	# Multiplying by 2 to add a buffer
	P_HIB = A_R * panel_heating * 2
	P_CHARGE = P_HIB / P_hib_ratio

	# Get approximate percentage active time
	percentage_active_best_case = np.arccos(np.sqrt(2) * r_moon / R) / np.pi

	# Get required battery capacity
	mu_moon = PhysicalConstants.mu_moon
	v_orbit = np.sqrt(mu_moon / R)
	omega_orbit = v_orbit / R
	T_ORBIT = 2 * np.pi / omega_orbit
	T_hibernation = P_hib_ratio * 0.5 * T_ORBIT
	E_CAPACITY_ROVER = T_hibernation * P_HIB * UnitConversions.J_to_Whr / DoD

	# Get required pointing
	W_P = Z * pointing_accuracy + W_R
	ETA_LINK = W_R ** 2 / W_P ** 2

	# Get laser power
	P_LAS = P_CHARGE / (laser_efficiency * ETA_LINK) * UnitConversions.W_to_kW
	P_LAS /= (1 - np.exp(-2))

	# Apply constraints - laser exists, pointing is manageable, power is low, active time is achievable at altitude
	active_time_mask = percentage_active_best_case < P_hib_ratio
	wavelength = 1070e-9
	W_B = np.sqrt(ETA_LINK * W_R ** 2)
	transmitter_mask = W_B ** 4 > 4 * (Z * wavelength / np.pi) ** 2
	efficiency_mask = ETA_LINK < eta_link_requirement
	power_mask = P_LAS > laser_requirement * UnitConversions.W_to_kW
	charge_power_mask = P_CHARGE > 2 * P_OP
	for res_array in [W_R, A_R, P_HIB, ETA_LINK, P_LAS, E_CAPACITY_ROVER, P_CHARGE]:
		res_array[efficiency_mask] = 0.0
		res_array[power_mask] = 0.0
		res_array[transmitter_mask] = 0.0
		res_array[active_time_mask] = 0.0
		res_array[charge_power_mask] = 0.0

	fig, ax = plt.subplots(3, 2, figsize=(14, 7), sharex=True)
	if isinstance(a, np.ndarray):
		num_contours = 100
		A *= UnitConversions.m_to_km
		im = ax[0, 0].contourf(A, P_OP, A_R, num_contours)
		ax[0, 0].set_title("Receiver area [$m^2$]")
		ax[0, 0].set_ylabel("$P_{op}$ [W]")
		fig.colorbar(im, ax=ax[0, 0])

		im = ax[1, 0].contourf(A, P_OP, P_HIB, num_contours)
		ax[1, 0].set_title("Hibernation power [W]")
		ax[1, 0].set_ylabel("$P_{op}$[W]")
		fig.colorbar(im, ax=ax[1, 0])
		
		im = ax[2, 0].contourf(A, P_OP, ETA_LINK * 100, num_contours)
		ax[2, 0].set_title("Link Efficiency")
		ax[2, 0].set_ylabel("$P_{op}$ [W]")
		ax[2, 0].set_xlabel("Altitude [km]")
		fig.colorbar(im, ax=ax[2, 0])	

		im = ax[0, 1].contourf(A, P_OP, P_LAS, num_contours)
		ax[0, 1].set_title("Laser Power [kW]")
		ax[0, 1].set_ylabel("$P_{op}$ [W]")
		fig.colorbar(im, ax=ax[0, 1])

		im = ax[1, 1].contourf(A, P_OP, E_CAPACITY_ROVER, num_contours)
		ax[1, 1].set_title("Rover Battery Capacity [Whr]")
		ax[1, 1].set_ylabel("$P_{op}$ [W]")
		fig.colorbar(im, ax=ax[1, 1])

		im = ax[2, 1].contourf(A, P_OP, P_CHARGE, num_contours)
		ax[2, 1].set_title("Charge Power")
		ax[2, 1].set_ylabel("Power required [W]")
		ax[2, 1].set_xlabel("Altitude [km]")
		fig.colorbar(im, ax=ax[2, 1])
	else:
		im = ax[0, 0].plot(P_OP, A_R)
		ax[0, 0].set_title("Receiver area [$m^2$]")

		im = ax[1, 0].plot(P_OP, P_HIB)
		ax[1, 0].set_title("Hibernation power [W]")
		
		im = ax[2, 0].plot(P_OP, ETA_LINK * 100)
		ax[2, 0].set_title("Link Efficiency [%]")
		ax[2, 0].set_xlabel("Operating power [W]")

		im = ax[0, 1].plot(P_OP, P_LAS)
		ax[0, 1].set_title("Laser Power [kW]")

		im = ax[1, 1].plot(P_OP, E_CAPACITY_ROVER)
		ax[1, 1].set_title("Rover Battery Capacity [Whr]")

		im = ax[2, 1].plot(P_OP, P_CHARGE)
		ax[2, 1].set_title("Charge Power [W]")
		ax[2, 1].set_xlabel("Operating power [W]")

	plt.show()


if __name__ == '__main__':
	pointing_accuracy = 5e-7
	a=None

	altitudes_for_pointing(pointing_accuracy, a=a)

	