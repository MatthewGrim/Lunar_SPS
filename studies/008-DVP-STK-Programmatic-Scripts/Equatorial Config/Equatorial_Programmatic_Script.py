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

import os
import time

import numpy as np
from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *


def generate_stk_connect_commands(semi_maj_axis, eccentricity, orbit_data, time_step, study_name, file_path):

    # This function takes in the various data points which are going to be passed to STK,
    # and then it writes a series of Connect Commands for autonomously varying the satellite
    # orbit parameters, and then printing the resulting illumination and access events in reports

    time_start = time.time()
    # New file with this name is generated in the current folder
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'w') as fh:
        # Generates new report of lighting events at target
        fh.write('ReportCreate */Target/Target1 Type Export Style "Lighting_Times" File '
                 '"{}\DVP_{}_Target_Lighting.csv"\n'.format(file_path, study_name))
        for i in range(len(semi_maj_axis)):
            # Sets new orbit for satellite, varying semi major axis and eccentricity
            fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "01 Jul 2008 10:00:00.000" "30 Jun 2010 '
                     '10:00:00.000" {} J2000 "1 Jul 2008 10:00:00.000" {} {} '
                     '0 0 0 360\n'.format(time_step, semi_maj_axis[i]*1000, eccentricity[i]))
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
    time_end = time.time()
    print('Time required to write connect commands: {} seconds'.format(round(time_end - time_start, 5)))


def run_stk_v2(scenario_path, study_name):
    # This function opens an instance of STK, loads the desired scenario, and executes the
    # connect commands written by the previous functions

    from win32api import GetSystemMetrics
    # from IPython.display import Image, display, SVG
    import os
    import comtypes
    from comtypes.client import CreateObject
    from comtypes.client import GetActiveObject

    print('Opening STK...')
    # Open new instance of STK
    # app = CreateObject("STK11.Application")
    # Pass open instance of STK
    app = GetActiveObject('svchost.Application')
    app.Visible = True
    app.UserControl = True
    app.Top = 0
    app.Left = 0
    app.Width = int(GetSystemMetrics(0) / 2)
    app.Height = int(GetSystemMetrics(1) - 30)

    root = app.Personality2

    comtypes.client.gen_dir
    os.listdir(comtypes.client.gen_dir)

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
    loop_start = time.time()

    duration = np.zeros(size)

    j = 0
    for i in range(size):
        time_start = time.time()
        root.ExecuteCommand(commands[i])
        time_end = time.time()
        if i == 0:
            print('Generating target lighting time report...')
        else:
            if j == 0:
                print('Adjusting Satellite orbit...')
                j += 1
            elif j == 1:
                print('Generating SPS access report...')
                j += 1
            elif j == 2:
                print('Generating SPS range report...')
                j += 1
            elif j == 3:
                print('Generating SPS lighting report...')
                j = 0
        # Print progress update
        print('Progress: {}%, Execution Time: {} seconds'.format(round(i * 100.0 / (size - 1), 2),
                                                                 round(time_end - time_start, 5)))
        duration[i] = time_end - time_start
    loop_end = time.time()

    print('Total time to generate data: {} minutes'.format((loop_end - loop_start) / 60.0))
    print('Average command execution time: {} seconds'.format(np.mean(duration)))


def main():

    # Set resolution of data points in km
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Time step for satellite (not for scenario) in STK
    time_step = 60

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Pathway to STK scenario for this study
    scenario_path = "{}\STK-Scenarios\EquatorialOrbit\Brandhorst_Figure3.sc".format(main_directory)
        
    # Name of study
    study_name = 'Equatorial_IncrementedRes'

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
    print('Getting orbit data...')
    sma, ecc, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # Generate connect commands for programmatically running STK simulations
    print('Writing connect commands...')
    generate_stk_connect_commands(sma, ecc, orbit_data, time_step, study_name, new_path)

    # Open STK, load scenario, and execute commands to create data set
    # run_stk_v2(scenario_path, study_name)


main()
