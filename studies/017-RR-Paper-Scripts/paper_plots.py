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


def get_range_elevation_plot():
	elevation_dir = 'elevation'
	range_dir = 'range'

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

	fig, ax = plt.subplots(2, figsize=figsize, sharex=True)
	colors = ['b', 'r', 'g']
	labels = ['Lat $0^\circ$', 'Lat $15^\circ$', 'Lat $30^\circ$']
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

	ax[0].axvline(-24.52 / 2, color='grey', linestyle='--')
	ax[0].axvline(24.52 / 2, color='grey', linestyle='--')
	ax[0].axhline(10.0, color='k', linestyle='-')
	ax[0].set_xlim([-50, 50])
	ax[0].set_ylim([-10, 100])
	ax[0].set_ylabel('Elevation [$^\circ$]')
	ax[0].grid(linestyle='--')
	ax[0].legend()

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
	ax[1].axhline(2200.0, color='k', linestyle='-')
	ax[1].set_xlim([-50, 50])
	ax[1].set_ylim([1200, 2600])
	ax[1].set_ylabel('Range [$km$]')
	ax[1].set_xlabel('Time [$min$]')
	ax[1].grid(linestyle='--')
	# ax[1].legend()

	fig.tight_layout()
	plt.savefig('link_range_elevation')
	plt.show()


def get_rover_temperature():
	rover_dir = 'rover_temp'

	fig, ax = plt.subplots(2, figsize=figsize, sharex=True)

	colors = ['b', 'r', 'g']
	max_power = [60.0, 160.0]
	for i, case in enumerate(['_cold', '_hot']):
		heater = np.loadtxt(os.path.join(rover_dir, 'heater{}.csv'.format(case)), delimiter=',')
		laser = np.loadtxt(os.path.join(rover_dir, 'laser{}.csv'.format(case)), delimiter=',')
		panel = np.loadtxt(os.path.join(rover_dir, 'panel{}ii.csv'.format(case)), delimiter=',')

		panel[:, 1] *= 350.0 / max_power[i]
		panel[:, 1] -= 200
		panel = panel[panel[:, 0].argsort()]

		lns1 = ax[i].plot(heater[:, 0], heater[:, 1], label='heater', c='b')
		lns2 = ax[i].plot(laser[:, 0], laser[:, 1], label='laser', c='g')
		ax[i].set_ylabel('Power [$W$]')
		ax[i].set_ylim([0.0, max_power[i]])

		ax2 = ax[i].twinx()
		lns3 = ax2.plot(panel[:, 0], panel[:, 1], label='panel', c='r')
		ax2.set_ylabel('Temperature [$^\circ C$]')
		ax2.set_ylim([-200, 150])
		
		ax[i].set_xlim([9.5, 10.0])
		ax[i].grid(linestyle='--')

		if i == 0:
			lns = lns1+lns2+lns3
			labs = [l.get_label() for l in lns]
			ax[i].legend(lns, labs, loc=0)

	ax[1].set_xlabel("Time [$Earth\ days$]")
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
	ax.set_xlabel('Time [$hrs$]')
	ax.set_ylabel('Temperature [$^\circ C$]')
	ax.legend()

	fig.tight_layout()
	plt.savefig('laser_temp')
	plt.show()


if __name__ == '__main__':
    figsize = (14, 7)
    get_range_elevation_plot()
    # get_rover_temperature()
    # get_laser_temperature()

