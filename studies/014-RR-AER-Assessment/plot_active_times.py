""""
Author: Rohan Ramasamy
Date: 12/10/2018

This scripts plots the active times for different satellite constellations
"""

from numpy import unravel_index

from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignFunctions import *
from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *


def plot_constellation_active_times():
    kwargs = {
        "resolutions": np.array((50.0, 100.0, 100.0, 250.0)),
        "thresholds": np.array((1000.0, 1500.0, 2500.0)),
        "min_perigee": 800.0,
        "max_perigee": 5000.0,
        "max_apogee": 5000.0
    }


    fig, ax = plt.subplots(2, 3, sharex='col', sharey='row', figsize=(15, 8))
    for i, study_name in enumerate(['NorthPole_IncrementedRes_Generic', 'Equatorial_IncrementedRes_Generic']):
        study = study_initialization(study_name, **kwargs)
        for num_sps in range(1, 4):
            # Get pathway to main Lunar_SPS directory
            current_folder = os.getcwd()
            parent_folder = os.path.dirname(current_folder)
            main_directory = os.path.dirname(parent_folder)
            stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

            # --- READ IN DATA FILES ---
            data_set = read_in_processed_data_reports(stk_data_path, study_name, num_sps)

            data_set['total_active_time'] = [100.0 * j / study['duration'] for j in data_set['total_active_time']]

            # Reorganize the data lists into 2D arrays
            sorted_data_set, perigee_altitudes, unique_perigees, apogee_altitudes, unique_apogees = sort_data_lists(data_set, study['orbits'], study_name, **kwargs)

            im = ax[i, num_sps - 1].contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['total_active_time'], 500)
            fig.colorbar(im, ax=ax[i, num_sps - 1])
            ax[i, num_sps - 1].set_xlabel("Apolune Radius")
        name = "Polar" if "NorthPole" in study_name else "Equatorial"
        ax[i, 0].set_ylabel("{}\n Perilune Radius".format(name))
    plt.savefig('active_times')
    plt.show()
    


if __name__ == '__main__':
    plot_constellation_active_times()

