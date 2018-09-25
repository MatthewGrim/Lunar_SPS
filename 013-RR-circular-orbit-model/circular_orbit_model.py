"""
Author: Rohan Ramasamy
Date: 25/09/2018

This script contains functions to approximate the necessary access time, orbit altitude and power requirements for a
particular target, assuming the delivering satellite has a circular orbit.
"""

import numpy as np
from matplotlib import pyplot as plt


def get_orbit_period(altitude):
    radius = altitude
    omega = np.sqrt(G * m_moon / radius ** 3)

    return 2 * np.pi / omega


def get_access_time(r_orbit, omega):
    theta_access = 2 * np.arccos(r_moon / r_orbit)
    theta_eclipse = 2 * np.arcsin(r_moon / r_orbit)
    t_access = (theta_access - theta_eclipse) / omega

    return t_access


if __name__ == '__main__':
    # Define constants
    r_moon = 1737e3
    m_moon = 7.34767309e22
    G = 6.67408e-11
    seconds_to_hours = 1 / 3600.0
    m_to_km = 1e-3

    # Set parameter space
    r_orbit = np.linspace(800e3, 5000e3, 1000) + r_moon

    # Get access time
    T = get_orbit_period(r_orbit)
    t_access = get_access_time(r_orbit, 2 * np.pi / T)

    # Get required access time
    p_operation = 100.0
    p_hibernation = 7.2
    p_ratio = p_hibernation / (p_operation - p_hibernation)

    # Plot orbital periods
    T *= seconds_to_hours
    t_access *= seconds_to_hours
    r_orbit -= r_moon
    r_orbit *= m_to_km
    plt.figure()
    plt.plot(r_orbit, t_access / (T - t_access))
    plt.axhline(p_ratio, linestyle="--")
    plt.show()

