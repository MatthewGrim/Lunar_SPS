"""
Author: Rohan Ramasamy
Date: 03/12/2018

This file contains physical constants to be used in this project
"""

class UnitConversions(object):
    km_to_m = 1e3
    m_to_km = 1e-3

    W_to_kW = 1e-3
    kW_to_W = 1e3
    J_to_Whr = 1.0 / 3600.0
    Whr_to_J = 3600.0

    radian_to_microradian = 1e6
    microradian_to_radian = 1e-6

    def __init__(self):
        raise RuntimeError('This is a static class!')

