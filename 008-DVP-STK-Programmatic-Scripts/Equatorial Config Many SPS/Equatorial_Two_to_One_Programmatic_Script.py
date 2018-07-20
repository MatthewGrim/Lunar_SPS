""""
15/05/2018
Author: Darian van Paridon

This file is used to get the orbital elements (semi-major axis and eccentricity)
based on the perigee and apogee of an orbit, print the data into a list of connect commands for STK
and then run the commands through STK to generate data.

The scenario programmed in this script is the collection of access and lighting times for an SPS constellations
in a common lunar polar orbit, targeting the south pole.

This particular script assumes that all orbits have been simulated for a single SPS, and thus simulations for
constellations of size two or greater are considered.

Furthermore, this version of the script refers to an STK instance and scenario which is already open, as opposed
to versions of the programmatic scripts which open STK and load the desired scenario
"""

from DVP_Programmatic_Functions import *


def vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee):

    # This function generates a series of orbital data points, in terms of apogee/perigee which
    # is then converted into semi major axis/eccentricity for execution in STK

    radius_moon = 1737.0

    resolution = np.array((10.0, 25.0, 50.0, 100.0))
    thresholds = np.array((100.0, 250.0, 1000.0))
    orbit_data = [0, 0]
    perigee = 0.0
    while perigee <= max_perigee:
        if 0.0 <= perigee < thresholds[0]:
            peri_step = resolution[0]
        elif thresholds[0] <= perigee < thresholds[1]:
            peri_step = resolution[1]
        elif thresholds[1] <= perigee < thresholds[2]:
            peri_step = resolution[2]
        elif thresholds[2] <= perigee:
            peri_step = resolution[3]
        perigee = perigee + peri_step
        apogee = perigee
        while apogee <= max_apogee:
            orbit_data = np.vstack((orbit_data, [perigee + radius_moon, apogee + radius_moon]))
            if 0.0 <= apogee < thresholds[0]:
                apo_step = resolution[0]
            elif thresholds[0] <= apogee < thresholds[1]:
                apo_step = resolution[1]
            elif thresholds[1] <= apogee < thresholds[2]:
                apo_step = resolution[2]
            elif thresholds[2] <= apogee:
                apo_step = resolution[3]
            apogee += apo_step

    eccentricity = np.zeros(len(orbit_data) - 1)
    semi_maj_axis = np.zeros(len(orbit_data) - 1)
    for i in range(0, len(orbit_data) - 1):
        eccentricity[i] = ((orbit_data[i + 1][1] / orbit_data[i + 1][0]) - 1) / (1 + (orbit_data[i + 1][1] / orbit_data[i + 1][0]))
        semi_maj_axis[i] = orbit_data[i + 1][0] / (1 - eccentricity[i])

    return semi_maj_axis, eccentricity, orbit_data


def generate_stk_connect_commands(semi_maj_axis, eccentricity, orbit_data, time_step, study_name, file_path):

    # This function takes in the various data points which are going to be passed to STK,
    # and then it writes a series of Connect Commands for autonomously varying the satellite
    # orbit parameters, and then printing the resulting illumination and access events in reports

    # New file with this name is generated in the current directory
    time_start = time.time()
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'w') as fh:
        # loop through orbits
        for i in range(len(semi_maj_axis)):
            # Connect commands for generating reports require a save location
            # Sets new orbit for satellite, varying semi major axis, eccentricity, and mean anomaly
            fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "17 May 2018 10:00:00.000" "17 May 2020 '
                     '10:00:00.000" {} Inertial "17 May 2018 10:00:00.000" {} {} '
                     '0 {} 0 0 \n'.format(time_step, semi_maj_axis[i] * 1000.0, eccentricity[i], k))
            fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File "{}\DVP_{}_{}perigee{}apogee_180.0argperi_access.csv" AccessObject '
                '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # Generates new report of SPS-to-target range statistics during access periods
            fh.write('ReportCreate */Satellite/SPS1 Type Save Style "Access_Range_Stats" File "{}\DVP_{}_{}perigee{}apogee_180.0argperi_range.txt" AccessObject '
                '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1]))
            # Generates new report of lighting times
            fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File "{}\DVP_{}_{}perigee{}apogee_180.0argperi_lighting.csv"\n'.format(
                    file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1]))

    time_end = time.time()
    print('Time required to write connect commands: {} seconds'.format(time_end - time_start))

    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'r') as fh:
        commands = fh.readlines()
        for i in range(len(commands)):
            print(commands[-i])


def run_stk_v2(scenario_path, study_name, orbit_data, stk_data_path):

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

    from comtypes.gen import STKUtil
    from comtypes.gen import STKObjects

    print('Loading scenario...')
    # Load predefined scenario, which contains satellite for which the orbit is varied
    # and the target for which the access time is calculated
    root.LoadScenario(r'{}'.format(scenario_path))
    sc = root.CurrentScenario
    sc2 = sc.QueryInterface(STKObjects.IAgScenario)

    print('Executing commands...')
    # Open file with connect commands, and execute them sequentially
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'r') as fh:
        commands = fh.readlines()
        size = len(commands)
    loop_start = time.time()

    duration = np.zeros(size)

    j = 0
    for i in range(size):
        time_start = time.time()
        if j == 0:
            print('Adjusting Satellite orbit...')
            root.ExecuteCommand(commands[i])
            j += 1

        elif j == 1:
            print('Generating SPS access report...')
            if not os.path.exists('{}\DVP_{}_{}perigee{}apogee_{}meananom_access.csv'.format(stk_data_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], 180.0)):
                root.ExecuteCommand(commands[i])
            else:
                print('Access report for {} x {} km orbit already exists'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
            j += 1

        elif j == 2:
            print('Generating SPS range report...')
            if not os.path.exists('{}\DVP_{}_{}perigee{}apogee_{}meananom_range.txt'.format(stk_data_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], 180.0)):
                root.ExecuteCommand(commands[i])
            else:
                print('Range report for {} x {} km orbit already exists'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
            j += 1

        elif j == 3:
            print('Generating SPS lighting report...')
            if not os.path.exists('{}\DVP_{}_{}perigee{}apogee_{}meananom_lighting.csv'.format(stk_data_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], 180.0)):
                root.ExecuteCommand(commands[i])
            else:
                print('Lighting for {} x {} km orbit already exists'.format(orbit_data[i + 1][0], orbit_data[i + 1][1]))
            j = 0

        time_end = time.time()
        # Print progress update
        print('Progress: {}%, Execution Time: {} seconds'.format(round(i * 100.0 / (size - 1), 2), round(time_end - time_start, 5)))
        duration[i] = time_end - time_start
    loop_end = time.time()

    print('Total time to generate data: {} minutes'.format((loop_end - loop_start) / 60.0))
    print('Average command execution time: {} seconds'.format(np.mean(duration)))


def main():

    # Set resolution of data points in km
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Set time step size for satellite in STK
    # Note that this does not set the time step for the scenario. Change this by hand.
    time_step = 60

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Pathway to scenario for this study
    scenario_path = '{}/STK-Scenarios/PolarOrbit/SouthPole_ManyToOne.sc'.format(main_directory)

    # Name of study - be descriptive
    study_name = 'SouthPole_IncrementedRes_ManytoOne'

    # Create folder inside main directory for storing data sets
    print('Creating new folder to store reports...')
    stk_data_path = r'{}/STK Data'.format(main_directory)
    if not os.path.exists(stk_data_path):
        os.makedirs(stk_data_path)

    # Create folder inside the data storage folder which is specific to this study
    new_path = r'{}/{}'.format(stk_data_path, study_name)
    if not os.path.exists(new_path):
        os.makedirs(new_path)

    # Get set of orbit data, varying apogee and perigee
    print('Calculating orbit data...')
    sma, ecc, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # Write list of commands for STK
    print('Writing connect commands...')
    generate_stk_connect_commands(sma, ecc, orbit_data, time_step, study_name, new_path)

    # Open STK, load scenario, and execute commands to create data set
    # run_stk_v2(scenario_path, study_name, orbit_data, stk_data_path)


main()
