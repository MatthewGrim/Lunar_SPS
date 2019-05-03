"""
Author: Rohan Ramasamy
Date: 05/12/2018

This script contains code to define possible SPS target designs
"""

import numpy as np
import matplotlib.pyplot as plt


def altitudes_for_pointing(pointing_accuracy, solar_efficiency=0.2, laser_efficiency=0.5, P_hib_ratio=0.04,
	                       eta_link_requirement=0.01, laser_requirement=10e3):
	a = np.linspace(8e5, 2e6, 500)
	r_moon = 1737.0e3
	r = a + r_moon

	P_rec = np.linspace(20.0, 500.0, 1000)
	
	# Get receiver size and hibernation power
	A, P_REC = np.meshgrid(a, P_rec, indexing='ij')
	R = A + r_moon
	Z = np.sqrt(r_moon ** 2 + R ** 2 - np.sqrt(2) * r_moon * R)
	Z_MAX = np.sqrt(R ** 2 - r_moon ** 2)
	Z = Z_MAX
	W_R = np.sqrt(P_REC / (np.pi * solar_efficiency * 1367.0))
	P_REC_HIB = P_rec * P_hib_ratio

	# Get approximate percentage active time
	percentage_active_best_case = np.arccos(np.sqrt(2) * r_moon / R) / np.pi

	# Get required battery capacity
	mu_moon = 6.674e-11 * 7.347673e22
	v_orbit = np.sqrt(mu_moon / R)
	omega_orbit = v_orbit / R
	T_ORBIT = 2 * np.pi / omega_orbit
	T_hibernation = P_hib_ratio * 0.5 * T_ORBIT
	E_CAPACITY = T_hibernation * P_REC_HIB / 3600.0

	# Get required pointing
	W_P = Z * pointing_accuracy + W_R
	ETA_LINK = W_R ** 2 / W_P ** 2

	# Get laser power
	P_LAS = P_REC / (laser_efficiency * ETA_LINK) * 1e-3
	P_LAS /= (1 - np.exp(-2))

	# Apply constraints - laser exists, pointing is manageable, power is low, active time is achievable at altitude
	active_time_mask = percentage_active_best_case < P_hib_ratio
	wavelength = 1070e-9
	W_B = np.sqrt(ETA_LINK * W_R ** 2)
	transmitter_mask = W_B ** 4 > 4 * (Z * wavelength / np.pi) ** 2
	efficiency_mask = ETA_LINK < eta_link_requirement
	power_mask = P_LAS > laser_requirement * 1e-3
	for res_array in [W_R, ETA_LINK, P_LAS, E_CAPACITY]:
		res_array[efficiency_mask] = 0.0
		res_array[power_mask] = 0.0
		res_array[transmitter_mask] = 0.0
		res_array[active_time_mask] = 0.0

	fig, ax = plt.subplots(4, figsize=(12, 7), sharex=True)

	num_contours = 100
	A *= 1e-3
	im = ax[0].contourf(A, P_REC, W_R, num_contours)
	ax[0].set_title("Receiver radius [$m$]")
	ax[0].set_ylabel("Power required [W]")
	
	fig.colorbar(im, ax=ax[0])
	im = ax[1].contourf(A, P_REC, ETA_LINK, num_contours)
	ax[1].set_title("Link Efficiency")
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
	pointing_accuracy = 8e-7

	altitudes_for_pointing(pointing_accuracy)

	