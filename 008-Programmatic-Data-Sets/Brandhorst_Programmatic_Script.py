""""
15/05/2018
Author: Darian van Paridon

This file is used to get the orbital elements (semi-major axis and eccentricity)
based on the perigee and apogee of an orbit, print the data into a list of connect commands for STK
and then run the commands through STK to generate data. The scenario programmed in this script
is the collection of access and lighting times for an SPS in the Brandhorst lunar SPS configuration
for a variety of orbits.

Note the location of STK scenario, and location of reports being generated, and current python directory.
These will need to be updated for each user/machine.

"""

import numpy as np


def get_orbital_elements():

    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK

    radius_moon = 1737.0
    # Set resolution of data points in km
    resolution = 1000.0
    max_perigee = 10000.0
    max_apogee = 55000.0
    orbit_data = [0.0, 0.0]

    for j in range(1, int(max_perigee / resolution) + 1):
        perigee = j * resolution
        apogee = perigee
        while apogee < max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            apogee += resolution

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])

    return semi_maj_axis, eccentricity, orbit_data


def generate_stk_connect_commands(semi_maj_axis, eccentricity, orbit_data):

    # This function takes in the various data points which are going to be passed to STK,
    # and then it writes a series of Connect Commands for autonomously varying the satellite
    # orbit parameters, and then printing the resulting illumination and access events in reports

    # New file with this name is generated in the current directory
    with open('CC_OrbitGeometry.txt', 'w') as fh:
        for i in range(len(semi_maj_axis)):
            # Sets new orbit for satellite, varying semi major axis and eccentricity
            fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "01 Jul 2008 10:00:00.000" "30 Jun 2010 '
                     '10:00:00.000" 3600 J2000 "1 Jul 2008 10:00:00.000" {} {} '
                     '0 0 0 360\n'.format(semi_maj_axis[i]*1000, eccentricity[i]))
            # Connect commands for generating reports require a save location
            # Generates new report of access time to target
            fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File "c:\Users\Darian van Paridan'
                     '\My Documents\STKData\{}perigee{}apogee_access.csv" AccessObject '
                     '*/Target/Target1\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # Generates new report of lighting times
            fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File "c:\Users'
                     '\Darian van Paridan\My Documents\STKData\{}perigee{}apogee_'
                     'lighting.csv"\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # Generates reports of penumbra and umbra times, which together make the eclipses
            # fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Penumbra_Times" File "c:\Users'
            #          '\Darian van Paridan\My Documents\STKData\{}perigee{}apogee_'
            #          'penumbra.csv"\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Umbra_Times" File "c:\Users'
            #          '\Darian van Paridan\My Documents\STKData\{}perigee{}apogee_'
            #          'umbra.csv"\n'.format(orbit_data[i + 1][0], orbit_data[i + 1][1])


def run_stk_v2():

    # This function opens an instance of STK, loads the desired scenario, and executes the
    # connect commands written by the previous functions

    from win32api import GetSystemMetrics
    # from IPython.display import Image, display, SVG
    import os as os
    import comtypes
    from comtypes.client import CreateObject

    # Open new instance of STK
    app = CreateObject("STK11.Application")
    app.Visible = True
    app.UserControl = True
    app.Top = 0
    app.Left = 0
    app.Width = int(GetSystemMetrics(0) / 2)
    app.Height = int(GetSystemMetrics(1) - 30)

    root = app.Personality2

    comtypes.client.gen_dir
    os.listdir(comtypes.client.gen_dir)

    from comtypes.gen import STKUtil
    from comtypes.gen import STKObjects

    # Load predefined scenario, which contains satellite for which the orbit is varied
    # and the target for which the access time is calculated
    root.LoadScenario(r'C:\Lunar_SouthPole_SPS.sc')
    sc = root.CurrentScenario
    sc2 = sc.QueryInterface(STKObjects.IAgScenario)
    # sc2.SetTimePeriod("1 Jul 2008 10:00:00", "30 Jun 2010 10:00:00")

    # Open file with connect commands, and execute them sequentially
    with open('CC_OrbitGeometry.txt', 'r') as fh:
        commands = fh.readlines()
        size = len(commands)
    for i in range(len(commands)):
        # Print progress update
        print('Progress: {}%'.format(i * 100.0 / size))
        print('\n')
        root.ExecuteCommand(commands[i])


def main():

    sma, ecc, orbit_data = get_orbital_elements()

    generate_stk_connect_commands(sma, ecc, orbit_data)

    run_stk_v2()


main()
