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


def generate_stk_connect_commands(semi_maj_axis, eccentricity, orbit_data, number_of_sps, arg_perigee, time_step, study_name, file_path):

    # This function takes in the various data points which are going to be passed to STK,
    # and then it writes a series of Connect Commands for autonomously varying the satellite
    # orbit parameters, and then printing the resulting illumination and access events in reports

    # New file with this name is generated in the current directory
    time_start = time.time()
    with open('CC_{}_OrbitStudy.txt'.format(study_name), 'w') as fh:
        # loop through orbits
        for i in range(len(semi_maj_axis)):
            # loop through constellation sizes
            conditions1 = range(2, number_of_sps[i] + 1)
            if number_of_sps[i] == 0 or number_of_sps[i] == 1:
                continue
            elif number_of_sps[i] == 2:
                j = 2
                k = 180.0
                # Connect commands for generating reports require a save location
                # Generates new report of access time to target
                fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File "{}\DVP_{}_{}perigee{}apogee_{}sps_{}argperi_access.csv" AccessObject '
                    '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], j,
                                                k))
                # Generates new report of SPS-to-target range statistics during access periods
                fh.write('ReportCreate */Satellite/SPS1 Type Save Style "Access_Range_Stats" File "{}\DVP_{}_{}perigee{}apogee_{}sps_{}argperi_range.txt" AccessObject '
                    '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], j,
                                                k))
                # Generates new report of lighting times
                fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File "{}\DVP_{}_{}perigee{}apogee_{}sps_{}argperi_lighting.csv"\n'.format(
                        file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], j, k))
                # Sets new orbit for satellite, varying semi major axis, eccentricity, and mean anomaly
                fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "17 May 2018 10:00:00.000" "17 May 2020 '
                         '10:00:00.000" {} Inertial "17 May 2018 10:00:00.000" {} {} '
                         '0 {} 0 0\n'.format(time_step, semi_maj_axis[i] * 1000.0, eccentricity[i], k))
            else:
                for j in range(2, number_of_sps[i] + 1):
                    # loop through mean anomalies of satellites in constellation
                    for k in arg_perigee['{}sps'.format(j)]:
                        # Connect commands for generating reports require a save location
                        # Generates new report of access time to target
                        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Access_Modified" File "{}\DVP_{}_{}perigee{}apogee_{}sps_{}argperi_access.csv" AccessObject '
                                 '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], j, k))
                        # Generates new report of SPS-to-target range statistics during access periods
                        fh.write('ReportCreate */Satellite/SPS1 Type Save Style "Access_Range_Stats" File "{}\DVP_{}_{}perigee{}apogee_{}sps_{}argperi_range.txt" AccessObject '
                                 '*/Target/Target1\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], j, k))
                        # Generates new report of lighting times
                        fh.write('ReportCreate */Satellite/SPS1 Type Export Style "Lighting_Times" File "{}\DVP_{}_{}perigee{}apogee_{}sps_{}argperi_lighting.csv"\n'.format(file_path, study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], j, k))
                        # Sets new orbit for satellite, varying semi major axis, eccentricity, and mean anomaly
                        fh.write('SetState */Satellite/SPS1 Classical J4Perturbation "17 May 2018 10:00:00.000" "17 May 2020 '
                                 '10:00:00.000" {} Inertial "17 May 2018 10:00:00.000" {} {} '
                                 '0 {} 0 0 \n'.format(time_step, semi_maj_axis[i]*1000.0, eccentricity[i], k))
    time_end = time.time()
    print('Time required to write connect commands: {} seconds'.format(time_end - time_start))


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
        root.ExecuteCommand(commands[i])
        time_end = time.time()
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
    scenario_path = '{}/STK-Scenarios/EquatorialOrbit/Equatorial_Orbit_45N_Target.sc'.format(main_directory)

    # Name of study - be descriptive
    study_name = 'Equatorial_IncrementedRes_ManytoOne'

    # Name of reference study on which constellation sizes are based
    reference_study = 'Equatorial_IncrementedRes'

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

    # Get maximum size of constellation required
    print('Calculating constellation sizes...')
    number_of_sps, arg_perigee = determine_constellation_size(ecc, reference_study)

    # Generate connect commands for programmatically running STK simulations
    print('Writing connect commands....')
    generate_stk_connect_commands(sma, ecc, orbit_data, number_of_sps, arg_perigee, time_step, study_name, new_path)

    # Open STK, load scenario, and execute commands to create data set
    # run_stk_v2(scenario_path, study_name)


main()
