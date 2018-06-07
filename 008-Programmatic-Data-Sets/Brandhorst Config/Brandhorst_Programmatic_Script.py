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
import os


def vary_orbital_elements(resolution):

    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK

    print('Getting orbital data...')

    radius_moon = 1737.0
    # Set resolution of data points in km
    max_perigee = 10000.0
    min_perigee = 0.0
    max_apogee = 55000.0
    orbit_data = [0.0, 0.0]

    for j in range(1, int((max_perigee - min_perigee) / resolution) + 1):
        perigee = min_perigee + (j * resolution)
        apogee = perigee
        while apogee <= max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            apogee += resolution

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])

    return semi_maj_axis, eccentricity, orbit_data


def generate_stk_connect_commands(semi_maj_axis, eccentricity, orbit_data, study_name, file_path):

    # This function takes in the various data points which are going to be passed to STK,
    # and then it writes a series of Connect Commands for autonomously varying the satellite
    # orbit parameters, and then printing the resulting illumination and access events in reports

    print('Writing connect commands...')

    # New file with this name is generated in the current folder
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'w') as fh:
        # Generates new report of lighing events at target
        fh.write('ReportCreate */Target/Target1 Type Export Style "Lighting_Times" File '
                 '"{}\DVP_{}_Target_Lighting.csv"\n'.format(file_path, study_name))
        for i in range(len(semi_maj_axis)):
            # Sets new orbit for satellite, varying semi major axis and eccentricity
            fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "01 Jul 2008 10:00:00.000" "30 Jun 2010 '
                     '10:00:00.000" 3600 J2000 "1 Jul 2008 10:00:00.000" {} {} '
                     '0 0 0 360\n'.format(semi_maj_axis[i]*1000, eccentricity[i]))
            # Connect commands for generating reports require a save location
            # Current location is outside of repository directory to avoid overloading with data files
            # Generates new report of access time to target
            fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File '
                     '"{}\DVP_{}_{}perigee_{}apogee_access.csv" AccessObject '
                     '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # Generates new report of SPS-to-target range during access periods
            fh.write('ReportCreate */Satellite/SPS1 Type Save Style "Access_Range_Stats" File '
                     '"{}\DVP_{}_{}perigee_{}apogee_range.txt" AccessObject '
                     '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # Generates new report of lighting times
            fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File '
                     '"{}\DVP_{}_{}perigee_{}apogee_lighting.csv"\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1]))


def run_stk_v2(scenario_path, study_name):
    # This function opens an instance of STK, loads the desired scenario, and executes the
    # connect commands written by the previous functions

    from win32api import GetSystemMetrics
    # from IPython.display import Image, display, SVG
    import os as os
    import comtypes
    from comtypes.client import CreateObject

    print('Opening STK...')
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

    print('Loading scenario...')
    # Load predefined scenario, which contains satellite for which the orbit is varied
    # and the target for which the access time is calculated
    root.LoadScenario(r'{}'.format(scenario_path))
    sc = root.CurrentScenario
    sc2 = sc.QueryInterface(STKObjects.IAgScenario)
    # sc2.SetTimePeriod("1 Jul 2008 10:00:00", "30 Jun 2010 10:00:00")

    print('Executing connect commands...')
    # Open file with connect commands, and execute them sequentially
    connect_command_file = 'CC_{}_OrbitStudy.txt'.format(study_name)
    with open(connect_command_file, 'r') as fh:
        commands = fh.readlines()
        size = len(commands)
    for i in range(len(commands)):
        # Print progress update
        print('Progress: {}%'.format(round(i * 100.0 / (size-1), 2)))
        root.ExecuteCommand(commands[i])


def main():

    # Step size between data points in km
    resolution = 500.0

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Pathway to STK scenario for this study
    scenario_path = "{}\STK-Scenarios\EquatorialOrbit\Brandhorst_Figure3.sc".format(main_directory)
        
    # Name of study
    study_name = 'Equatorial_{}kmRes_5min_Check'.format(resolution)

    # Create folder inside My Documents for storing data sets
    print('Creating new folder to store data...')
    stk_data_path = r"{}\STK Data".format(main_directory)
    if not os.path.exists(stk_data_path):
        os.makedirs(stk_data_path)

    # Create folder inside the data storage folder which is specific to this study
    new_path = r'{}\STK Data\{}'.format(main_directory, study_name)
    if not os.path.exists(new_path):
        os.makedirs(new_path)

    # Get set of orbit data, varying apogee and perigee
    sma, ecc, orbit_data = vary_orbital_elements(resolution)

    # Generate connect commands for programmatically running STK simulations
    generate_stk_connect_commands(sma, ecc, orbit_data, study_name, new_path)

    # Open STK, load scenario, and execute commands to create data set
    run_stk_v2(scenario_path, study_name)


main()
