"""
Author: Rohan Ramasamy
Date: 29/11/2018

This file is used to compute the gaussian beam divergence intensity for a particular beam radius.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


def evaluate_w_z(w_0, wavelength, z):
	z_R = np.pi * w_0 ** 2 / wavelength
	w_z = w_0 * np.sqrt(1 + (z / z_R) ** 2)

	return w_z

def gaussian_beam_intensity(z, r, w_0, wavelength, P_in):
	w_z = evaluate_w_z(w_0, wavelength, z)
	I_0 = 2 * P_in / (np.pi * w_0 ** 2)
	I = I_0 * (w_0 / w_z) ** 2 * np.exp(-2 * r ** 2 / w_z ** 2)

	return I

def get_flux_transmitter_radius_dependence(w_t, z, wavelength, P_las, P_rec, A_rec):
	phi_rec = P_rec / A_rec
	phi_beam = P_las / (np.pi * w_t ** 2 * (1 + (z * wavelength / (np.pi * w_t ** 2)) ** 2))

	plt.figure()
	plt.plot(w_t, phi_beam)
	plt.axhline(phi_rec, linestyle='--')
	plt.show()

def get_transmitter_options():
	z = 2531e3
	wavelength = 859e-9
	P_rec = 43.0
	A_rec = 0.079
	P_las = 3.44e3
	w_t = np.linspace(0.05, 1.12, 100)
	get_flux_transmitter_radius_dependence(w_t, z, wavelength, P_las, P_rec, A_rec)

def get_intensity_profile():
	params = "AMALIA_1300_submicro"
	if params is "AMALIA_1300_submicro":
		ranges = [2183.26e3, 2531.95e3]
		r_radius = 2.5
		r = np.linspace(0.0, r_radius, 500)
		wavelength = 1070e-9
		w_0 = 0.8844
		P_in = 2.98e3
		P_rec = 200.0
		target_radius = 0.341
	elif params is "Sorato_2300_submicro":
		ranges = [3065150.0, 3644010.0]
		r_radius = 2.5
		r = np.linspace(0.0, r_radius, 500)
		wavelength = 1070e-9
		w_0 = 1.0562
		P_in = 4.29e3
		P_rec = 43.0
		target_radius = 0.158
	else:
		raise ValueError("Invalid parameter selection")
	r_moon = 1737e3
	sigma = 1e-7
	phi_req = P_rec / (np.pi * target_radius ** 2)

	fig, ax = plt.subplots(2, sharex=True, figsize=(12, 7))

	color_i = np.linspace(0, 1, len(ranges))
	for i, z in enumerate(ranges):
		c = plt.cm.viridis(color_i[i])

		I = gaussian_beam_intensity(z, r, w_0, wavelength, P_in)
		r_interp = interpolate.interp1d(I, r)
		I_interp = interpolate.interp1d(r, I)
		I_max = I_interp(0.0)

		r_e2 = r_interp(I_max / np.exp(2))
		ax[0].plot(r, I, label="{}$km$".format(z * 1e-3), color=c)
		ax[0].axvline(r_e2, linestyle='--', color=c)
		ax[0].axvline(sigma * z, linestyle=':', color=c)
		if i == 0:
			ax[0].axvline(target_radius, linestyle='-', color="k", label="Target Radius")
			ax[0].axhline(phi_req, linestyle='--', color="grey", label="$\phi_{req}$")
	
		I_grad = np.diff(I) / np.diff(r)

		ax[1].plot((r[0:-1] + r[1:]) / 2, I_grad, label="{}$km$".format(z * 1e-3), color=c)
		ax[1].axvline(r_e2, linestyle='--', color=c)
		ax[1].axvline(sigma * z, linestyle=':', color=c)
		if i == 0:
			ax[1].axvline(target_radius, linestyle='-', color="k", label="Target Radius")

	ax[0].set_xlim([0, r_radius])
	ax[0].set_ylabel("Normalised Intensity [$Wm^{-2}]$")
	ax[1].set_ylabel("Intensity gradient [$Wm^{-3}$]")
	ax[1].set_xlabel("Offset from beam centre [m]")
	ax[0].legend()
	ax[1].legend()
	fig.suptitle("Parameter selection: {}\nTarget Radius: {}$m$".format(params, target_radius))

	plt.savefig(params)
	plt.show()


if __name__ == '__main__':
	get_intensity_profile()
	# get_transmitter_options()

