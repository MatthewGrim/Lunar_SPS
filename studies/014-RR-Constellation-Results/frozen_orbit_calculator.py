"""
Author: Rohan Ramasamy
Date: 13/11/2018

This script is used to generate an equivalent frozen orbit to the data generated by polar orbits.
"""

from matplotlib import pyplot as plt
import numpy as np 
from scipy import interpolate
import os


def read_frozen_orbits():
	main_dir = "frozen_orbit_data"
	data_file = "frozen_orbit.csv"
	data_path = os.path.join(main_dir, data_file)
	apogee_inclination_data = np.loadtxt(data_path, delimiter=",")
	print(apogee_inclination_data)
	apogee_inclination_interpolator = interpolate.interp1d(apogee_inclination_data[:, 0], apogee_inclination_data[:, 1])

	return apogee_inclination_interpolator

def get_inclination_transform(i_op):
	axis_transform = 6.8 / 180.0 * np.pi
	i_ep = np.arccos(np.cos(axis_transform) * np.cos(i_op) - np.sin(axis_transform) * np.sin(i_op))

	return i_ep


def get_inertial_inclination():
	i_op = np.linspace(0.0, np.pi / 2, 100)
	i_op = 56.2 / 180.0 * np.pi
	i_ep = get_inclination_transform(i_op)

	if isinstance(i_ep, np.ndarray):
		plt.figure()
		plt.plot(i_op * 180 / np.pi, i_ep * 180 / np.pi)
		plt.show()
	else:
		print(i_ep * 180 / np.pi)

def get_eccentricities(apogee_inclination_data, perigee):
	eccentricity = (apogee_inclination_data[:, 0] - perigee) / (apogee_inclination_data[:, 0] + perigee)
	return eccentricity

if __name__ == '__main__':
	get_inertial_inclination()

	apogee_inclination_interpolator = read_frozen_orbits()
	r_moon = 1737.4
	r_a_equivalent_polar = 5e3
	r_p = 1000.0
	r_a = r_a_equivalent_polar
	r_a_old = 0.0
	while not np.isclose(r_a, r_a_old):
		r_a_old = r_a
		inclination = apogee_inclination_interpolator(r_a)
		r_a = r_a_equivalent_polar / np.sin(inclination / 180.0 * np.pi)
		ecc = (r_a - r_p) / (r_a + r_p + 2 * r_moon)

	print("Apogee radius: ", r_a + r_moon)
	print("Inclination: ", inclination)
	print("Eccentricity", ecc)
	print("Perigee radius: ",  (1 - ecc) / (1 + ecc) * (r_moon + r_a))
	print("Equivalent polar orbit:", r_a_equivalent_polar, (1 - ecc) / (1 + ecc) * (r_moon + r_a_equivalent_polar) - r_moon)

