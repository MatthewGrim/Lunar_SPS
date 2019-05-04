""""
Author: Rohan Ramasamy
Date: 12/10/2018

This scripts plots the active times for different satellite constellations
"""

from numpy import unravel_index

from Lunar_SPS.pysrc.sps_design_tool.SPS_Constellation_DesignFunctions import *
from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def plot_constellation_active_times():
    kwargs = {
        "resolutions": np.array((50.0, 100.0, 100.0, 250.0)),
        "thresholds": np.array((1000.0, 1500.0, 2500.0)),
        "min_perigee": 800.0,
        "max_perigee": 5000.0,
        "max_apogee": 5000.0
    }

    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    study_name = 'Equatorial_IncrementedRes_Generic'
    total_duration = (end - start).total_seconds()
    current_folder = os.getcwd()
    parent_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(parent_folder)
    stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

    target_lighting_raw = os.path.join(stk_data_path, 'DVP_{}_Target_Lighting.csv'.format(study_name))
    target_lighting = parse_csv_to_array(target_lighting_raw, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)
    eclipse_sum = np.sum(target_eclipse[2])

    fig, ax = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(15, 8))
    for i, study_name in enumerate(['Equatorial_IncrementedRes_Generic']):
        study = study_initialization(study_name, **kwargs)
        for num_sps in range(1, 4):
            # Get pathway to main Lunar_SPS directory
            current_folder = os.getcwd()
            parent_folder = os.path.dirname(current_folder)
            main_directory = os.path.dirname(parent_folder)
            stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

            # --- READ IN DATA FILES ---
            data_set = read_in_processed_data_reports(stk_data_path, study_name, num_sps, use_storage=True)

            data_set['total_active_time'] = [100.0 * j / eclipse_sum for j in data_set['total_active_time']]

            # Reorganize the data lists into 2D arrays
            sorted_data_set, perigee_altitudes, unique_perigees, apogee_altitudes, unique_apogees = sort_data_lists(data_set, study['orbits'], study_name, **kwargs)

            im = ax[num_sps - 1].contourf(apogee_altitudes, perigee_altitudes, sorted_data_set['total_active_time'], 500)
            im.set_clim(0.0, 100.0)
            fig.colorbar(im, ax=ax[num_sps - 1])
            ax[num_sps - 1].set_title('Number of Satellites: {}'.format(num_sps ))
            ax[num_sps - 1].set_xlabel("Apolune altitude [$km$]", fontsize=font_size)
        name = "Polar" if "NorthPole" in study_name else "Equatorial"
        ax[0].set_ylabel("Perilune altitude [$km$]".format(name), fontsize=font_size)
    fig.tight_layout()
    plt.savefig('active_times')
    plt.show()
    


if __name__ == '__main__':
    font_size = 14
    plot_constellation_active_times()

