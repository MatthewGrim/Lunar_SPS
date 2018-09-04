"""
Author: Rohan Ramasamy
Date: 04/09/2018

This script processes the data obtained programmatically from STK for evaluating
the optimal SPS orbit for the polar lunar SPS configuration based on feasible pointing and power
constraints.

"""

import numpy as np

from Lunar_SPS.pysrc.post_process_functions.process_stk_sims import process_stk_data
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import convert_string_to_datetime

if __name__ == '__main__':
    # Set main variables of study
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    study_name = 'Equatorial_IncrementedRes_Generic'
    constellation_variable = 'arg_peri'

    # Set kwargs for orbital variables
    kwargs = {
        "resolutions": np.array((50.0, 100.0, 100.0, 250.0)),
        "thresholds": np.array((1000.0, 1500.0, 2500.0)),
        "max_perigee": 5000.0,
        "max_apogee": 5000.0,
        "min_apogee": 800.0
    }

    process_stk_data(3, study_name, constellation_variable, start, end, **kwargs)

