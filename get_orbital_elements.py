""""
15/05/2018
Author: Darian van Paridon

This file is used to get the orbital elements (semi-major axis and eccentricity)
based on the perigee and apogee of an orbit.

"""

import numpy as np
import matplotlib.pyplot as plt

radius_moon = 1737.0
apogee = 0.0
perigee = 0.0
orbit_data = [0.0, 0.0]

for j in range(0, (10000 / 250) + 1):
    perigee = j * 250.0
    apogee = perigee
    while apogee < 55000.0:
        orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
        apogee += 250.0

eccentricity = np.zeros(len(orbit_data) - 1)
semi_maj_axis = np.zeros(len(orbit_data) - 1)
for i in range(0, len(orbit_data) - 1):
    eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
    semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])

plt.subplot(211)
plt.plot(eccentricity)
plt.subplot(212)
plt.plot(semi_maj_axis)
plt.show()


with open('CC_OrbitGeometry.txt', 'w') as fh:
    for i in range(len(semi_maj_axis)):
        # Sets new orbit for satellite, varying semi major axis and eccentricity
        fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "01 Jul 2008 10:00:00.000" "30 Jun 2010 '
                 '10:00:00.000" 36000 J2000 "1 Jul 2008 10:00:00.000" {} {} '
                 '98.5 0.0 139.7299 360.0\n'.format(semi_maj_axis[i], eccentricity[i]))
        # Generates new report of access time to target
        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File "c:\Users\Darian van Paridan'
                 '\My Documents\STKData\{}perigee{}apogee_access.csv" AccessObject '
                 '*/Target/Target1\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
        # Generates new report of lighting times
        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File "c:\Users\Darian van Paridan'
                 '\My Documents\STKData\{}perigee{}apogee_'
                 'lighting.csv"\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
        # Generates reports of penumbra and umbra times, which together make the eclipses
        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Penumbra_Times" File "c:\Users\Darian van Paridan'
                 '\My Documents\STKData\{}perigee{}apogee_'
                 'penumbra.csv"\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Umbra_Times" File "c:\Users\Darian van Paridan'
                 '\My Documents\STKData\{}perigee{}apogee_'
                 'umbra.csv"\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))

