"""
Author: Rohan Ramasamy
Date: 02/12/2018

This script is used to assess the potential of using low altitude orbits to power target.
"""

import numpy as np
import matplotlib.pyplot as plt

from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignFunctions import rover_metrics


def get_battery_options(P_hib, P_op, phi_req):
	a = np.linspace(8e5, 3e6, 100)
	r_moon = 1737.0e3
	r = a + r_moon
	mu_moon = 4.9048695e12
	v_orbit = np.sqrt(mu_moon / r)
	omega_orbit = v_orbit / r
	T_orbit = 2 * np.pi / omega_orbit
	T_charge = T_orbit / 2
	theta = 2 * (np.pi / 2 - np.arccos(r_moon / r))
	charge_ratio = theta / np.pi

	theta_eclipse = (2 * np.arcsin(r_moon / r))
	T_eclipse = theta_eclipse / (2 * np.pi) * T_orbit
	plt.figure()
	plt.plot(a * 1e-3, T_eclipse / 3600)
	plt.title("Eclipse time")
	plt.show()

	percentage_active_best_case = np.arccos(np.sqrt(2) * r_moon / r) / np.pi
	percentage_active_worst_case = (np.arccos(np.sqrt(2) * r_moon / r) - np.arcsin(r_moon / r)) / np.pi
	plt.figure()
	plt.plot(a * 1e-3, 0.5 * (percentage_active_best_case + percentage_active_worst_case))
	plt.axvline(1300, linestyle='--')
	plt.title("% Active")
	plt.show()

	plt.figure()
	plt.plot(a * 1e-3, T_orbit / 3600)
	plt.title("Orbit Period [hrs]")
	plt.show()

	plt.figure()
	plt.plot(a * 1e-3, charge_ratio)
	plt.title("Ratio of battery use to charge time")
	plt.show()

	n_cycles = 24 / T_orbit * 3600 * 365 * 10
	reduction_factor = 2 * np.arcsin(r_moon / r) / (2 * np.arccos(np.sqrt(2) * r_moon / r))
	reduction_factor[reduction_factor > 1] = 1
	n_cycles *= reduction_factor
	plt.figure()
	plt.plot(a * 1e-3, n_cycles)
	plt.title("Number of battery cycles")
	plt.show()

	wavelength = 1070e-9
	z = np.sqrt(r_moon ** 2 + r ** 2 - np.sqrt(2) * r_moon * r)

	plt.figure()
	plt.plot(a * 1e-3, z * 1e-3)
	plt.title("Approximate range [km]")
	plt.show()

	w_t = np.sqrt(z * wavelength / np.pi)
	plt.figure()
	plt.plot(a * 1e-3, w_t)
	plt.title("Transmitter radius [m]")
	plt.show()

	w_b = np.sqrt(2) * w_t
	P_las = phi_req * np.pi * w_b ** 2 / 0.875

	specific_energy = 140.0
	T_req = T_orbit * 0.5 * P_hib / P_op
	E_req = P_las * T_req / 3600.0
	battery_mass = E_req / 140.0

	plt.figure()
	plt.plot(a * 1e-3, P_las * 1e-3)
	plt.title("Laser power [kW]")
	plt.show()

	plt.figure()
	plt.plot(a * 1e-3, E_req)
	plt.title("Battery capacity [Whr]")
	plt.show()

	plt.figure()
	plt.plot(a * 1e-3, battery_mass)
	plt.title("Battery mass [kg]")
	plt.show()

	i = np.arcsin(r_moon / r)
	plt.figure()
	plt.plot(a * 1e-3, i / np.pi * 180)
	plt.title("Inclination to avoid eclipse")
	plt.show()

	theta = np.arcsin(r_moon / (np.sqrt(2) * z))
	misalignment_loss = np.cos(theta)
	plt.figure()
	plt.plot(a * 1e-3, misalignment_loss)
	plt.title("Loss in transmission power due to misalignement")
	plt.show()


if __name__ == '__main__':
	rover_name = "sorato"
	rover = rover_metrics(rover_name)
	P_hib = rover["hibernation_pwr"]
	P_op = rover["operation_pwr"]
	phi_req = 550.0

	get_battery_options(P_hib, P_op, phi_req)

