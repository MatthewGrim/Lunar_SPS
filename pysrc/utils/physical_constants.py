"""
Author: Rohan Ramasamy
Date: 03/12/2018

This file contains physical constants to be used in this project
"""

class PhysicalConstants(object):
    r_moon = 1737.0
    mu_moon = 6.674e-11 * 7.347673e22
    solar_irradiance = 1367.0

    def __init__(self):
        raise RuntimeError('This is a static class!')

