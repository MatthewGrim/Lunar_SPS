import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import CubicSpline, interp1d

r_moon = 1743e3
r_orbit =1300e3 + r_moon

def get_misalignment_efficiency(interp_points, elevation, azimuth, range_0):
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

def get_link_efficiency(interp_points):
    interpolator = CubicSpline(range_0[:, 0], range_0[:, 1])
    interp_r = interpolator(interp_points)

    transmitter_radius = 0.597
    target_radius = np.sqrt(1.83 / np.pi)
    wavelength = 1070e-9
    pointing_error = 800e-9
    minimum_beam_radius = interp_r * pointing_error
    min_efficiency = target_radius ** 2 / minimum_beam_radius ** 2 * (1 - np.exp(-2))

    return min_efficiency


fontsize=14
labels = ['Lat $0^\circ$', 'Lat $15^\circ$', 'Lat $30^\circ$']
colors = ['b', 'r', 'g']
fig, ax = plt.subplots(2, figsize=(14, 7), sharex=True)
for i, lat in enumerate([0, 15, 30]):
    elevation = np.loadtxt(os.path.join('elevation', '{}_degree_access.csv'.format(lat)), delimiter=',')
    azimuth = np.loadtxt(os.path.join('azimuth', '{}degree_access_azimuth.csv'.format(lat)), delimiter=',')
    range_0 = np.loadtxt(os.path.join('range', '{}_degree_access_range.csv'.format(lat)), delimiter=',')
    range_0[:, 1] *= 1e3

    # Interpolate data points
    interp_min = max(azimuth[:, 0].min(), elevation[:, 0].min(), range_0[:, 0].min())
    interp_max = min(azimuth[:, 0].max(), elevation[:, 0].max(), range_0[:, 0].max())
    interp_points = np.linspace(interp_min, interp_max, 101)

    # Calculate misalignment efficiency
    misalignment_efficiency = get_misalignment_efficiency(interp_points, elevation, azimuth, range_0)
    ax[0].plot(interp_points, misalignment_efficiency, label=labels[i], c=colors[i])

    # Get average misalignment loss
    angle_interp = interp1d(interp_points * 60, misalignment_efficiency)
    times = np.linspace(-24.52 * 30, 24.52 * 30, 200)
    integral = np.trapz(angle_interp(times), dx=times[1]-times[0])
    average_misalignment_efficiency = integral / (24.52 * 60)
    minimum_misalignment_efficiency = angle_interp(24.52 * 30)
    print(r'Average misalignment efficiency for {} degree latitude: {} %'.format(lat, average_misalignment_efficiency * 100))
    print(r'Minimum misalignment efficiency for {} degree latitude: {} %'.format(lat, minimum_misalignment_efficiency * 100))

    min_link_efficiency = get_link_efficiency(interp_points)
    ax[1].plot(interp_points, min_link_efficiency, c=colors[i])

    efficiency_interp = interp1d(interp_points * 60, min_link_efficiency)
    times = np.linspace(-24.52 * 30, 24.52 * 30, 200)
    integral = np.trapz(efficiency_interp(times), dx=times[1]-times[0])
    average_link_efficiency = integral / (24.52 * 60)
    minimum_link_efficiency = efficiency_interp(24.52 * 30)
    print('Average link efficiency: {} %'.format(average_link_efficiency))
    print('Minimum link efficiency: {} %'.format(minimum_link_efficiency * 100))

    print('Total efficiency: {} %'.format(minimum_link_efficiency * minimum_misalignment_efficiency * 100.0))

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
    misalignment_efficiency = get_misalignment_efficiency(interp_points, elevation, azimuth, range_0)
    neg_mask = interp_points < 0
    pos_mask = interp_points > 0
    ax[0].plot(interp_points[neg_mask], misalignment_efficiency[neg_mask], c=colors[i], linestyle='--')
    ax[0].plot(interp_points[pos_mask], misalignment_efficiency[pos_mask], c=colors[i], linestyle='--')

    min_link_efficiency = get_link_efficiency(interp_points)
    ax[1].plot(interp_points[neg_mask], min_link_efficiency[neg_mask], c=colors[i], linestyle='--')
    ax[1].plot(interp_points[pos_mask], min_link_efficiency[pos_mask], c=colors[i], linestyle='--')


ax[0].grid(linestyle='--')
ax[1].grid(linestyle='--')
ax[0].axvline(-24.52 / 2, color='grey', linestyle='--')
ax[0].axvline(24.52 / 2, color='grey', linestyle='--')
ax[1].axvline(-24.52 / 2, color='grey', linestyle='--')
ax[1].axvline(24.52 / 2, color='grey', linestyle='--')
ax[0].legend()
ax[0].set_ylabel('$\eta_{misalignment}$', fontsize=fontsize)
ax[1].set_ylabel('$\eta_{link}$', fontsize=fontsize)
ax[1].set_xlabel('Time  [min]', fontsize=fontsize)
fig.tight_layout()
plt.savefig('transmission_efficiency')
plt.show()