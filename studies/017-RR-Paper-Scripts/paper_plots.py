"""
Author: Rohan Ramasamy
Date: 23/03/2019

Functions to generate the plots for the paper. Data is digitised from original plots from studies and so is smoothed with a 
cubic spline interpolation.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d


r_moon = 1743e3
r_orbit =1300e3 + r_moon

def get_misalignment_efficiency(interp_points, elevation, azimuth, range_0, lat):
    interpolator = CubicSpline(range_0[:, 0], range_0[:, 1])
    interp_r = interpolator(interp_points)
    interpolator = CubicSpline(azimuth[:, 0], azimuth[:, 1]) if lat != 0.0 else interp1d(azimuth[:, 0], azimuth[:, 1])
    interp_a = interpolator(interp_points) 
    if lat == 0.0:
        interpolator = interp1d(elevation[:, 0], elevation[:, 1])
    else:
        interpolator = CubicSpline(elevation[:, 0], elevation[:, 1])
    interp_e = interpolator(interp_points)

    interp_e *= np.pi / 180
    interp_a *= np.pi / 180

    # Convert coordinates into cartesian frame with axes normal to surface
    z = interp_r * np.sin(interp_e)
    y = interp_r * np.cos(interp_e) * np.cos(interp_a)
    x = interp_r * np.cos(interp_e) * np.sin(interp_a)

    # Rotate coordinates so that z axis is within the lunar equatorial plane
    rotation_angles = [-lat * np.pi / 180.0, 0.0, 0.0]
    c_x = np.cos(rotation_angles[0])
    s_x = np.sin(rotation_angles[0])
    c_y = np.cos(rotation_angles[1])
    s_y = np.sin(rotation_angles[1])
    c_z = np.cos(rotation_angles[2])
    s_z = np.sin(rotation_angles[2])

    # Arbitrary 3D rotation transform
    R = np.asarray([[c_y * c_z, c_z * s_x * s_y - c_x * s_z, c_x * c_z * s_y + s_x * s_z],
                   [c_y * s_z, c_x * c_z + s_x * s_y * s_z, -c_z * s_x + c_x * s_y * s_z],
                   [-s_y, c_y * s_x, c_x * c_y]])

    x_rot = R[0, 0] * x + R[0, 1] * y + R[0, 2] * z
    y_rot = R[1, 0] * x + R[1, 1] * y + R[1, 2] * z
    z_rot = R[2, 0] * x + R[2, 1] * y + R[2, 2] * z

    # Recalculate elevation, azimuth
    elevation_rot = np.arcsin(z_rot / interp_r)
    azimuth_rot = np.arctan(x_rot / y_rot)

    # Get values in degrees
    interp_e *= 180.0 / np.pi
    elevation_rot *= 180.0 / np.pi
    interp_a *= 180.0 / np.pi
    azimuth_rot *= 180.0 / np.pi

    # Calculate efficiency after change in rotation
    misalignment_efficiency = np.cos((90 - elevation_rot) /180.0 * np.pi)

    return misalignment_efficiency

def get_link_efficiency(interp_points, range_0):
    interpolator = CubicSpline(range_0[:, 0], range_0[:, 1])
    interp_r = interpolator(interp_points)

    transmitter_radius = 0.597
    target_radius = np.sqrt(1.83 / np.pi)
    wavelength = 1070e-9
    pointing_error = 800e-9
    minimum_beam_radius = interp_r * pointing_error
    min_efficiency = target_radius ** 2 / minimum_beam_radius ** 2 * (1 - np.exp(-2))

    return min_efficiency


def get_range_elevation_plot():
	elevation_dir = 'elevation'
	range_dir = 'range'
	azimuth_dir = 'azimuth'

	access_0 = np.loadtxt(os.path.join(elevation_dir, '0_degree_access.csv'), delimiter=',')
	access_15 = np.loadtxt(os.path.join(elevation_dir, '15_degree_access.csv'), delimiter=',')
	access_30 = np.loadtxt(os.path.join(elevation_dir, '30_degree_access.csv'), delimiter=',')

	no_access_0 = np.loadtxt(os.path.join(elevation_dir, '0_degree_no_access.csv'), delimiter=',')
	no_access_15 = np.loadtxt(os.path.join(elevation_dir, '15_degree_no_access.csv'), delimiter=',')
	no_access_30 = np.loadtxt(os.path.join(elevation_dir, '30_degree_no_access.csv'), delimiter=',')
    
    # Fix ordering in access_0 array
	access_0 = access_0[access_0[:, 0].argsort()]
	_, unique_indices = np.unique(access_0[:, 0], return_index=True)
	access_0 = access_0[unique_indices, :]

	fig, ax = plt.subplots(4, figsize=(10, 8), sharex=True)
	colors = ['b', 'r', 'g']
	labels = ['$\\theta_{lat}=0^\circ$', '$\\theta_{lat}=15^\circ$', '$\\theta_{lat}=30^\circ$']
	for i, array in enumerate([access_15, access_30]):
		interp_points = np.linspace(array[0, 0], array[-1, 0], 101)
		interpolator = CubicSpline(array[:, 0], array[:, 1])
		interp_y = interpolator(interp_points)
		
		ax[0].plot(interp_points, interp_y, color=colors[i+1], label=labels[i+1])
	ax[0].plot(access_0[:, 0], access_0[:, 1], color=colors[0], label=labels[0])

	# Get average misalignment loss
	angle_interp = interp1d(access_0[:, 0] * 60, access_0[:, 1])
	times = np.linspace(-24.52 * 30, 24.52 * 30, 200)
	integral = np.trapz(np.cos((90 - angle_interp(times)) / 180.0 * np.pi), dx=times[1]-times[0])
	print('Worst case misalignment loss: {}\%'.format(integral / (24.52 * 60) * 100))

	for i, array in enumerate([no_access_0, no_access_15, no_access_30]):
		array_neg = array[array[:, 0] < 0, :]
		array_pos = array[array[:, 0] > 0, :]
		ax[0].plot(array_neg[:, 0], array_neg[:, 1], linestyle='--', color=colors[i])
		ax[0].plot(array_pos[:, 0], array_pos[:, 1], linestyle='--', color=colors[i])

	range_0 = np.loadtxt(os.path.join(range_dir, '0_degree_access_range.csv'), delimiter=',')
	range_15 = np.loadtxt(os.path.join(range_dir, '15_degree_access_range.csv'), delimiter=',')
	range_30 = np.loadtxt(os.path.join(range_dir, '30_degree_access_range.csv'), delimiter=',')

	no_range_0 = np.loadtxt(os.path.join(range_dir, '0_degree_no_access_range.csv'), delimiter=',')
	no_range_15 = np.loadtxt(os.path.join(range_dir, '15_degree_no_access_range.csv'), delimiter=',')
	no_range_30 = np.loadtxt(os.path.join(range_dir, '30_degree_no_access_range.csv'), delimiter=',')
	for i, array in enumerate([range_0, range_15, range_30]):
		interp_points = np.linspace(array[0, 0], array[-1, 0], 25)
		interpolator = CubicSpline(array[:, 0], array[:, 1])
		interp_y = interpolator(interp_points)
		
		ax[1].plot(interp_points, interp_y, color=colors[i], label=labels[i])
		
	for i, array in enumerate([no_range_0, no_range_15, no_range_30]):
		array_neg = array[array[:, 0] < 0, :]
		array_pos = array[array[:, 0] > 0, :]
		ax[1].plot(array_neg[:, 0], array_neg[:, 1], linestyle='--', color=colors[i])
		ax[1].plot(array_pos[:, 0], array_pos[:, 1], linestyle='--', color=colors[i])

	for i, lat in enumerate([0, 15, 30]):
		print('\nLatitiude: {}'.format(lat))
		elevation = np.loadtxt(os.path.join('elevation', '{}_degree_access.csv'.format(lat)), delimiter=',')
		azimuth = np.loadtxt(os.path.join('azimuth', '{}degree_access_azimuth.csv'.format(lat)), delimiter=',')
		range_0 = np.loadtxt(os.path.join('range', '{}_degree_access_range.csv'.format(lat)), delimiter=',')
		range_0[:, 1] *= 1e3

		# Interpolate data points
		interp_min = max(azimuth[:, 0].min(), elevation[:, 0].min(), range_0[:, 0].min())
		interp_max = min(azimuth[:, 0].max(), elevation[:, 0].max(), range_0[:, 0].max())
		interp_points = np.linspace(interp_min, interp_max, 101)

		# Calculate misalignment efficiency
		misalignment_efficiency = get_misalignment_efficiency(interp_points, elevation, azimuth, range_0, lat) * 100
		ax[2].plot(interp_points, misalignment_efficiency, label=labels[i], c=colors[i])

		# Get average misalignment loss
		angle_interp = interp1d(interp_points * 60, misalignment_efficiency)
		times = np.linspace(-24.52 * 30, 24.52 * 30, 200)
		integral = np.trapz(angle_interp(times), dx=times[1]-times[0])
		average_misalignment_efficiency = integral / (24.52 * 60)
		minimum_misalignment_efficiency = angle_interp(24.52 * 30)
		print(r'Average misalignment efficiency: {} %'.format(average_misalignment_efficiency))
		print(r'Minimum misalignment efficiency: {} %'.format(minimum_misalignment_efficiency))

		min_link_efficiency = get_link_efficiency(interp_points, range_0) * 100
		ax[3].plot(interp_points, min_link_efficiency, c=colors[i])

		efficiency_interp = interp1d(interp_points * 60, min_link_efficiency)
		times = np.linspace(-24.52 * 30, 24.52 * 30, 200)
		integral = np.trapz(efficiency_interp(times), dx=times[1]-times[0])
		average_link_efficiency = integral / (24.52 * 60)
		minimum_link_efficiency = efficiency_interp(24.52 * 30)
		print('Average link efficiency: {} %'.format(average_link_efficiency))
		print('Minimum link efficiency: {} %'.format(minimum_link_efficiency))

		print('Total efficiency: {} %'.format(minimum_link_efficiency * minimum_misalignment_efficiency * 1e-2))

		# Add curves for when there is no access
		elevation = np.loadtxt(os.path.join('elevation', '{}_degree_no_access.csv'.format(lat)), delimiter=',')
		
		# Fix indexing in elevation array
		elevation = elevation[elevation[:, 0].argsort()]
		_, unique_indices = np.unique(elevation[:, 0], return_index=True)
		elevation = elevation[unique_indices, :]

		elevation = elevation[elevation[:, 0].argsort()]
		azimuth = np.loadtxt(os.path.join('azimuth', '{}degree_no_access_azimuth.csv'.format(lat)), delimiter=',')
		range_0 = np.loadtxt(os.path.join('range', '{}_degree_no_access_range.csv'.format(lat)), delimiter=',')
		range_0[:, 1] *= 1e3

		# Interpolate data points
		interp_neg_min = max(azimuth[azimuth[:, 0] < 0][:, 0].min(), elevation[elevation[:, 0] < 0][:, 0].min(), range_0[range_0[:, 0] < 0][:, 0].min())
		interp_neg_max = min(azimuth[azimuth[:, 0] < 0][:, 0].max(), elevation[elevation[:, 0] < 0][:, 0].max(), range_0[range_0[:, 0] < 0][:, 0].max())
		interp_pos_min = max(azimuth[azimuth[:, 0] > 0][:, 0].min(), elevation[elevation[:, 0] > 0][:, 0].min(), range_0[range_0[:, 0] > 0][:, 0].min())
		interp_pos_max = min(azimuth[azimuth[:, 0] > 0][:, 0].max(), elevation[elevation[:, 0] > 0][:, 0].max(), range_0[range_0[:, 0] > 0][:, 0].max())
		
		interp_neg_points = np.linspace(interp_neg_min, interp_neg_max, 101)
		interp_pos_points = np.linspace(interp_pos_min, interp_pos_max, 101)
		interp_points = np.concatenate((interp_neg_points, interp_pos_points))

		# Calculate misalignment efficiency
		misalignment_efficiency = get_misalignment_efficiency(interp_points, elevation, azimuth, range_0, lat) * 100
		neg_mask = interp_points < 0
		pos_mask = interp_points > 0
		ax[2].plot(interp_points[neg_mask], misalignment_efficiency[neg_mask], c=colors[i], linestyle='--')
		ax[2].plot(interp_points[pos_mask], misalignment_efficiency[pos_mask], c=colors[i], linestyle='--')

		min_link_efficiency = get_link_efficiency(interp_points, range_0) * 100
		ax[3].plot(interp_points[neg_mask], min_link_efficiency[neg_mask], c=colors[i], linestyle='--')
		ax[3].plot(interp_points[pos_mask], min_link_efficiency[pos_mask], c=colors[i], linestyle='--')

	for i in range(4):
		ax[i].axvline(-24.52 / 2, color='grey', linestyle='--')
		ax[i].axvline(24.52 / 2, color='grey', linestyle='--')
		ax[i].set_xlim([-42, 42])
		ax[i].grid(linestyle='--')
	ax[0].set_ylim([-10, 100])
	ax[0].axhline(10.0, color='k', linestyle='-')
	ax[0].set_ylabel('Elevation [$^\circ$]', fontsize=font_size)
	ax[0].legend()

	ax[1].axhline(2200.0, color='k', linestyle='-')
	ax[1].set_ylim([1200, 2600])
	ax[1].set_ylabel('Range [$km$]', fontsize=font_size)

	ax[2].set_ylabel('$\eta_{misalignment}$ [%]', fontsize=font_size)
	ax[2].set_ylim([0, 110])
	ax[3].set_ylabel('$\eta_{link}$ [%]', fontsize=font_size)
	ax[3].set_xlabel('Time  [min]', fontsize=font_size)

	fig.tight_layout()
	plt.savefig('link_range_elevation_efficiency')
	plt.show()


def get_rover_temperature():
	rover_dir = 'rover_temp'

	fig, ax = plt.subplots(2, figsize=figsize, sharex=True)

	colors = ['b', 'r', 'g']
	max_power = [250.0, 700.0]
	for i, case in enumerate(['_cold_total', '_hot_total']):
		heater = np.loadtxt(os.path.join(rover_dir, 'heater{}.csv'.format(case)), delimiter=',')
		laser = np.loadtxt(os.path.join(rover_dir, 'laser{}.csv'.format(case)), delimiter=',')
		panel = np.loadtxt(os.path.join(rover_dir, 'panel{}.csv'.format(case)), delimiter=',')

		panel = panel[panel[:, 0].argsort()]

		lns1 = ax[i].plot(heater[:, 0], heater[:, 1], label='heater power', c='b')
		lns2 = ax[i].plot(laser[:, 0], laser[:, 1], label='laser power', c='g')
		ax[i].set_ylabel('Power [$W$]', fontsize=font_size)
		ax[i].set_ylim([0.0, max_power[i]])

		ax2 = ax[i].twinx()
		lns3 = ax2.plot(panel[:, 0], panel[:, 1], label='panel temperature', c='r')
		ax2.set_ylabel('Temperature [$^\circ C$]', fontsize=font_size)
		ax2.set_ylim([-200, 150])
		
		ax[i].set_xlim([9.5, 10.0])
		ax[i].grid(linestyle='--')

		if i == 0:
			lns = lns1+lns2+lns3
			labs = [l.get_label() for l in lns]
			ax[i].legend(lns, labs, loc=0)

	ax[1].set_xlabel("Time [$Earth\ days$]", fontsize=font_size)
	fig.tight_layout()
	plt.savefig('rover_temp')
	plt.show()


def get_laser_temperature():
	rover_dir = 'laser_temp'

	heater = np.loadtxt(os.path.join(rover_dir, 'heater_temp.csv'), delimiter=',')
	laser = np.loadtxt(os.path.join(rover_dir, 'laser_temp.csv'), delimiter=',')

	heater = heater[heater[:, 0].argsort()]
	laser = laser[laser[:, 0].argsort()]

	seconds_to_hrs = 1 / 3600

	x_min = 36000
	x_max = 55000
	interp_pts = np.linspace(x_min, x_max, 240)
	heater_interp = interp1d(heater[:, 0], heater[:, 1])
	laser_interp = interp1d(laser[:, 0], laser[:, 1])
	heater = heater_interp(interp_pts)
	laser = laser_interp(interp_pts)
	average = (heater_interp(interp_pts) + laser_interp(interp_pts)) / 2.0

	fig, ax = plt.subplots(figsize=figsize)

	ax.plot(interp_pts * seconds_to_hrs, heater, label='radiator', c='b')
	ax.plot(interp_pts * seconds_to_hrs, laser, label='laser', c='r')
	ax.plot(interp_pts * seconds_to_hrs, average, label='average', c='grey', linestyle='--')
	ax.scatter(interp_pts * seconds_to_hrs, heater, c='b')
	ax.scatter(interp_pts * seconds_to_hrs, laser, c='r')
	ax.axhline(0.0, c='k')
	ax.axvline(12.58, linestyle='--', c='k')
	ax.axvline(12.94, linestyle='--', c='k')


	ax.grid(linestyle='--')
	ax.set_xlim([x_min  * seconds_to_hrs, x_max  * seconds_to_hrs])
	ax.set_xlabel('Time [$hrs$]', fontsize=font_size)
	ax.set_ylabel('Temperature [$^\circ C$]', fontsize=font_size)
	ax.legend()

	fig.tight_layout()
	plt.savefig('laser_temp')
	plt.show()


if __name__ == '__main__':
	font_size = 14
	figsize = (14, 7)
	get_range_elevation_plot()
	# get_rover_temperature()
	# get_laser_temperature()

