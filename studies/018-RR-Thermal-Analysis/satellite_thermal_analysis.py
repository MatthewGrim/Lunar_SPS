"""
Author: Rohan Ramasamy
Date: 10/04/2020

This file contains a thermal assessment of a solar power satellite. The calculations are based on the original work from 
study 017.
"""

import numpy as np
import matplotlib.pyplot as plt


class ThermalAnalysisConstants(object):
    """
    Structure to store the constants used in thermal analysis
    """
    stefan_boltzmann = 5.670374e-8
    solar_irradiance = 1370.0
    r_moon = 1731.0

    def __init__(self, panel_area=10.0, 
                 front_emissivity=0.78, front_absorptivity=0.05, 
                 back_emissivity=0.87, back_absorptivity=0.2, 
                 inclination=np.pi/3.0,
                 orbit_altitude=1300.0,
                 p_las=19.725e3, laser_efficiency=0.4, power_fraction=0.1,
                 targ_theta=3*np.pi/2):
        """
        Initialise design based constants
        """
        # Design based
        self.panel_area = panel_area                    # Area of individual radiator panels
        self.front_emissivity = front_emissivity        # Infrared emissivity of Moon facing side
        self.front_absorptivity = front_absorptivity    # Solar absorptivity of Moon facing side
        self.back_emissivity = back_emissivity          # Infrared emissivity of panel back
        self.back_absorptivity = back_absorptivity      # Solar absorptivity of panel back
        self.inclination = inclination                  # Tilt of panel from the lunar equatorial plane

        self.p_las = p_las
        self.laser_efficiency = laser_efficiency
        self.power_fraction = power_fraction

        self.targ_theta = targ_theta

        self.orbit_altitude = orbit_altitude            # Orbit altitude in km


def get_solar_heat(plot_result=False):
    # calculate surface projection of incident sunlight onto panel
    phi_tilt = ta_struct.inclination
    cos_elevation_angle = np.abs(np.sin(theta) * np.cos(phi_tilt))

    # calculate angles where Sun is blocked by Moon
    r_orbit = ta_struct.r_moon + ta_struct.orbit_altitude
    eclipse_theta = np.arcsin(ta_struct.r_moon / r_orbit)

    # Calculate thermal emissivity
    solar_absorptivity = np.zeros(theta.shape)
    solar_absorptivity[theta < np.pi] = ta_struct.back_absorptivity
    solar_absorptivity[theta >= np.pi] = ta_struct.front_absorptivity
    solar_absorptivity[np.logical_and(theta > 3 * np.pi / 2 - eclipse_theta, theta < 3 * np.pi / 2 + eclipse_theta)] = 0.0

    # calculate power
    p_sol = ta_struct.solar_irradiance * ta_struct.panel_area * solar_absorptivity * cos_elevation_angle

    if plot_result:
        plt.figure()
        plt.plot(theta, p_sol)
        plt.show()
    
    return p_sol

def get_lunar_heat(plot_result=False):
    phi_tilt = ta_struct.inclination
    surface_projection = np.cos(phi_tilt)

    # Approximate lunar temperature based on thermal equilibrium and known night temperature
    lunar_temperature = np.zeros(theta.shape)
    lunar_temperature[theta < np.pi] = (1370.0 * np.sin(theta[theta < np.pi]) / ta_struct.stefan_boltzmann) ** 0.25
    lunar_temperature[lunar_temperature < 100.0] = 100.0
    
    lunar_irradiance = ta_struct.stefan_boltzmann * lunar_temperature ** 4 * (ta_struct.r_moon / (ta_struct.r_moon + ta_struct.orbit_altitude)) ** 2
    p_lunar = lunar_irradiance * ta_struct.panel_area * ta_struct.front_emissivity * surface_projection

    if plot_result:
        plt.figure()
        plt.plot(theta, p_lunar)
        plt.show()

    return p_lunar

def get_laser_power(plot_result=False):
    # Get theta angles when laser is beaming power
    targ_theta = ta_struct.targ_theta
    power_fraction = ta_struct.power_fraction
    min_power_theta = targ_theta - np.pi * power_fraction
    max_power_theta = targ_theta + np.pi * power_fraction
    assert min_power_theta > 0.0
    assert max_power_theta < 2 * np.pi

    p_las = np.zeros(theta.shape)
    p_las[np.logical_and(theta > min_power_theta, theta < max_power_theta)] = ta_struct.p_las * ta_struct.laser_efficiency

    return p_las


if __name__ == '__main__':
    theta = np.linspace(0.0, 2 * np.pi, 1001)
    ta_struct = ThermalAnalysisConstants()

    p_sol = get_solar_heat(plot_result=False)
    p_lun = get_lunar_heat(plot_result=False)
    p_las = get_laser_power()

    fig, ax = plt.subplots(2)
    ax[0].plot(theta, p_sol)
    ax[0].plot(theta, p_lun)
    ax[1].plot(theta, p_las)
    plt.show()
