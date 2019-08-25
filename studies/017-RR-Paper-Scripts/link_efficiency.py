import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import CubicSpline, interp1d


range_dir = 'range'
colors = ['b', 'r', 'g']
labels = ['Lat $0^\circ$', 'Lat $15^\circ$', 'Lat $30^\circ$']
transmitter_radius = 0.597
wavelength = 1070e-9
pointing_error = 800e-9
target_radius = np.sqrt(1.83 / np.pi)

range_0 = np.loadtxt(os.path.join(range_dir, '0_degree_access_range.csv'), delimiter=',')
range_15 = np.loadtxt(os.path.join(range_dir, '15_degree_access_range.csv'), delimiter=',')
range_30 = np.loadtxt(os.path.join(range_dir, '30_degree_access_range.csv'), delimiter=',')

no_range_0 = np.loadtxt(os.path.join(range_dir, '0_degree_no_access_range.csv'), delimiter=',')
no_range_15 = np.loadtxt(os.path.join(range_dir, '15_degree_no_access_range.csv'), delimiter=',')
no_range_30 = np.loadtxt(os.path.join(range_dir, '30_degree_no_access_range.csv'), delimiter=',')

fig, ax = plt.subplots(3)
for i, array in enumerate([range_0, range_15, range_30]):
    interp_points = np.linspace(array[0, 0], array[-1, 0], 200)
    interpolator = CubicSpline(array[:, 0], array[:, 1])
    interp_y = interpolator(interp_points)
    
    surface_beam_radius = transmitter_radius * np.sqrt(1 + (wavelength * (interp_y * 1000.0) / (np.pi * transmitter_radius ** 2)) ** 2)
    minimum_beam_radius = interp_y * 1000.0 * pointing_error

    min_efficiency = target_radius ** 2 / minimum_beam_radius ** 2 * (1 - np.exp(-2))
    
    efficiency_interp = interp1d(interp_points * 60, min_efficiency)
    times = np.linspace(-24.52 * 30, 24.52 * 30, 200)
    integral = np.trapz(efficiency_interp(times), dx=times[1]-times[0])
    print('Average efficiency: {}\%'.format(integral / (24.52 * 60) * 100))
    print('Minimum efficiency: {}\%'.format(min(min_efficiency) * 100))

    ax[0].plot(interp_points, interp_y, color=colors[i], label=labels[i])
    ax[1].plot(interp_points, surface_beam_radius, color=colors[i], label=labels[i])
    ax[1].plot(interp_points, minimum_beam_radius, color=colors[i], linestyle='--')
    ax[1].axhline(target_radius)
    ax[2].plot(interp_points, target_radius ** 2 / minimum_beam_radius ** 2 * (1 - np.exp(-2)), linestyle='--')
plt.show()

