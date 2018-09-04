"""
Author: Rohan Ramasamy
Date: 21/08/2018

This script was written to fix mistakes made in the generation of the Equatorial Generic data. Two mistakes were made:

1. Data was ouput for the access times in epoch seconds from 01/07/2018 10:00:00. This means that the results are offset from the
   start of the simulation (17/05/2018 10:00:00) by 45 X 24 X 3600 = 3888000 this offset needs to be fixed so simulations start at
   0. This will allow current functions to work as they did previously.

2. The range data was output not just with the summary for the 35% of simulations. This leads to significantly more data being 
   produced.

"""

from Lunar_SPS.DVP_general_SPS_functions import *
from Lunar_SPS.DVP_Programmatic_Functions import *
import sys

def fix_access_times(sim_file_name, total_duration_offset):
    # Fix access times by adding offset to start and finish times
    access_times = np.loadtxt("{}_access.csv".format(sim_file_name), skiprows=1, delimiter=",")

    access_times[:, 0] += total_duration_offset
    access_times[:, 1] += total_duration_offset

    np.savetxt("{}_access_edit.csv".format(sim_file_name), access_times, delimiter=",",
               header="\"Start Time (EpSec)\", \"Stop Time (EpSec)\", \"Duration (sec)\"")

def fix_range_file(sim_file_name):
    # Fix range files by removing text irrelevant to global statistics
    range_file = "{}_range.txt".format(sim_file_name)
    new_range_file = "{}_range_edit.txt".format(sim_file_name)

    range_lines = open(range_file, "r")
    new_range_lines = open(new_range_file, "w")

    count = 0
    for line in range_lines:
        assert 0 <= count < 5
        if "Global Statistics" in line or "Section Statistics" in line:
            count = 5

        if count > 0:
            new_range_lines.write(line)
            count -= 1


def main(fix_access=True, fix_range=True, rename=False):
    # Get the ammount of time missing from the simulation epoch
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2018', '07', '31', '10', '0', '0.0'])
    total_duration_offset = (end - start).total_seconds()
    print(total_duration_offset)
    
    # Set bounds on parametric scan
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Get orbit data set
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee, 
                                                                         min_perigee=800.0, 
                                                                         resolutions=np.array((50.0, 100.0, 100.0, 250.0)),
                                                                         thresholds=np.array((1000.0, 1500.0, 2500.0)))

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(current_folder)

    # Name study
    study_name = 'Equatorial_IncrementedRes_Generic'
    reference_study = 'Equatorial_IncrementedRes'

    # Set file path for data
    stk_data_path = os.path.join(main_directory, 'STK Data', study_name)
    assert os.path.exists(stk_data_path)
    
    # Get maximum size of constellation required
    print('Calculating constellation sizes...')
    max_constellation_size = 3
    arg_perigee = set_constellation_size(v)
    
    # loop through orbits
    processed_sims = list()
    num_sims = 0
    missing_sims = 0
    for i in range(len(semi_maj_axis)):
        for j in range(1, max_constellation_size + 1):
            # loop through mean anomalies of satellites in constellation
            for k in arg_perigee['{}sps'.format(j)]:
                file_name = "DVP_{}_{}perigee{}apogee_{}argperi".format(study_name, orbit_data[i + 1][0], orbit_data[i + 1][1], k)
                sim_file_name = os.path.join(stk_data_path, file_name)
                
                # If it's already been processed, continue
                if sim_file_name in processed_sims:
                    continue
                print(file_name)
                num_sims += 1

                # Defensive check that all files have been generated
                access_file_name = "{}_access.csv".format(sim_file_name)
                range_file_name = "{}_range.txt".format(sim_file_name)
                lighting_file_name = "{}_lighting.csv".format(sim_file_name)
                assert os.path.exists(access_file_name)
                
                if fix_access:
                    fix_access_times(sim_file_name, total_duration_offset)
                if fix_range:
                    fix_range_file(sim_file_name)

                if rename:
                    old_access_file_name = "{}_access.csv".format(sim_file_name)
                    old_range_file_name = "{}_range.txt".format(sim_file_name)
                    new_access_file_name = "{}_access_edit.csv".format(sim_file_name)
                    new_range_file_name = "{}_range_edit.txt".format(sim_file_name)

                    assert os.path.exists(new_access_file_name)
                    assert os.path.exists(new_range_file_name)

                    os.rename(new_access_file_name, access_file_name)
                    os.rename(new_range_file_name, range_file_name)
                    
                processed_sims.append(sim_file_name)

    print("Number of Sims: {}".format(num_sims))
    print("Number of Missing Sims: {}".format(missing_sims))


if __name__ == '__main__':
    main(fix_access=True, fix_range=True, rename=True)

