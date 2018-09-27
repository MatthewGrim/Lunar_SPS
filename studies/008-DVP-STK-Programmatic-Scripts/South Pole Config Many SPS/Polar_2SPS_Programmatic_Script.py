""""
5/06/2018
Author: Darian van Paridon

This file is used to get the orbital elements (semi-major axis and eccentricity)
based on the perigee and apogee of an orbit, print the data into a list of connect commands for STK
and then run the commands through STK to generate data. The scenario programmed in this script
is the collection of access and lighting times for an SPS in a frozen lunar orbit, targeting the south pole.

"""

import numpy as np
import os
import time


def generate_stk_connect_commands(true_anomaly, semi_maj_axis, eccentricity, study_name, file_path):

    # This function takes in the various data points which are going to be passed to STK,
    # and then it writes a series of Connect Commands for autonomously varying the satellite
    # orbit parameters, and then printing the resulting illumination and access events in reports

    # New file with this name is generated in the current directory
    print('Writing connect commands....')
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'w') as fh:
        # Generates new report for illumination events at target
        fh.write('ReportCreate */Target/Target1 Type Export Style "Lighting_Times" File '
                 '"{}\DVP_{}_Target_Lighting.csv"\n'.format(file_path, study_name))
        # Create reports for first SPS, which is not altered
        fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "17 May 2018 10:00:00.000" "17 May 2020 '
                 '10:00:00.000" 3600 Inertial "17 May 2018 10:00:00.000" {} {} '
                 '90.01 90 0 360\n'.format(semi_maj_axis * 1000, eccentricity))
        # Create reports for second SPS, where the relative true anomaly is varied
        for i in range(len(true_anomaly)):
            # Sets new orbit for satellite, varying the true anomaly of satellite 2 relative to 1
            fh.write('SetState */Satellite/SPS2 Classical J4Perturbation "17 May 2018 10:00:00.000" "17 May 2020 '
                     '10:00:00.000" 3600 Inertial "17 May 2018 10:00:00.000" {} {} '
                     '89.99 90 0 {}\n'.format(semi_maj_axis * 1000, eccentricity, true_anomaly[i]))
            # Generates new report of access time to target
            fh.write('ReportCreate */Satellite/SPS2 Type Export Style "Access_Modified" File "{}\DVP_{}_SPS2_{}anomaly_access.csv" AccessObject '
                     '*/Target/Target1\n'.format(file_path, study_name, true_anomaly[i]))
            # Generates new report of lighting times
            fh.write('ReportCreate */Satellite/SPS2 Type Export Style "Lighting_Times" File "{}\DVP_{}_SPS2_{}anomaly_'
                     'lighting.csv"\n'.format(file_path, study_name, true_anomaly[i]))
            # Generates new report of SPS-to-target range during access periods
            fh.write('ReportCreate */Satellite/SPS2 Type Save Style "Access_Range_Stats" File "{}\DVP_{}_SPS2_{}anomaly_range.txt" AccessObject '
                     '*/Target/Target1\n'.format(file_path, study_name, true_anomaly[i]))
        # Generates new report of SPS access to target
        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File "{}\DVP_{}_SPS1_access.csv" AccessObject '
                 '*/Target/Target1\n'.format(file_path, study_name))
        # Generates new report of lighting times
        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File "{}\DVP_{}_SPS1_'
                 'lighting.csv"\n'.format(file_path, study_name))
        # Generates new report of SPS-to-target range during access periods
        fh.write('ReportCreate */Satellite/SPS1 Type Save Style "Access_Range_Stats" File "{}\DVP_{}_SPS1_range.txt" AccessObject '
                 '*/Target/Target1\n'.format(file_path, study_name))


def run_stk_v2(scenario_path, study_name):

    # This function opens an instance of STK, loads the desired scenario, and executes the
    # connect commands written by the previous functions

    from win32api import GetSystemMetrics
    # from IPython.display import Image, display, SVG
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

    duration = []
    print('Executing connect commands...')
    # Open file with connect commands, and execute them sequentially
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'r') as fh:
        commands = fh.readlines()
        size = len(commands)
    for i in range(len(commands)):
        # Print progress update
        time_start = time.time()
        root.ExecuteCommand(commands[i])
        time_end = time.time()
        print('Progress: {}%, Execution Time: {} seconds'.format(round((i + 1) * 100.0 / size, 2), round(time_end - time_start, 5)))
        duration.append(time_end - time_start)
    print('Average command execution time: {} seconds'.format(round(np.mean(duration), 5)))


def main():

    # Get list of relative true anomalies to simulate
    resolution = 3.6
    true_anomaly = np.linspace(0, 180, int(180 / resolution))

    apogee = 5000.0
    perigee = 50.0

    eccentricity = ((apogee / perigee) - 1) / (1 + (apogee / perigee))
    semi_maj_axis = perigee / (1 - eccentricity)

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Pathway to STK scenario for this study
    scenario_path = '{}/STK-Scenarios/PolarOrbit_MultiSPS/Scenario2.sc'.format(main_directory)

    # Name of study
    study_name = 'SouthPole_2SPS_{}degRes'.format(resolution)

    # Create folder inside main directory for storing data sets
    print('Creating new folder to store data...')
    stk_data_path = r'{}/STK Data'.format(main_directory)
    if not os.path.exists(stk_data_path):
        os.makedirs(stk_data_path)

    # Create folder inside the data storage folder which is specific to this study
    new_path = r'{}\{}'.format(stk_data_path, study_name)
    if not os.path.exists(new_path):
        os.makedirs(new_path)

    # Generate connect commands for programmatically running STK simulations
    generate_stk_connect_commands(true_anomaly, semi_maj_axis, eccentricity, study_name, new_path)

    # Open STK, load scenario, and execute commands to create data set
    run_stk_v2(scenario_path, study_name)


main()
