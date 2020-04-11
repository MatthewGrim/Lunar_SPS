"""
Author: Rohan Ramasamy
Date: 10/04/2020

This file contains a thermal assessment of a solar power satellite. The calculations are based on the original work from 
study 017.

Deployable Panel Radiator - 47th International Conference on Environmental Systems
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import integrate


class ThermalAnalysisConstants(object):
    """
    Structure to store the constants used in thermal analysis
    """
    stefan_boltzmann = 5.670374e-8
    solar_irradiance = 1370.0
    r_moon = 1731.0
    aluminium_heat_capacity = 920.0

    def __init__(self, panel_area=10.0, 
                 front_emissivity=0.78, front_absorptivity=0.05, 
                 back_emissivity=0.87, back_absorptivity=0.2, 
                 panel_inclination=np.pi/3.0,
                 orbit_altitude=1300.0,
                 p_las=19.725e3, laser_efficiency=0.4, power_fraction=0.1,
                 targ_theta=3*np.pi/2,
                 laser_mass=25.0, radiator_mass=100.0, number_of_radiators=2):
        """
        Initialise design based constants
        """
        # Design based
        self.panel_area = panel_area                    # Area of individual radiator panels
        self.front_emissivity = front_emissivity        # Infrared emissivity of Moon facing side
        self.front_absorptivity = front_absorptivity    # Solar absorptivity of Moon facing side
        self.back_emissivity = back_emissivity          # Infrared emissivity of panel back
        self.back_absorptivity = back_absorptivity      # Solar absorptivity of panel back
        self.panel_inclination = panel_inclination      # Tilt of panel from the lunar equatorial plane

        self.p_las = p_las
        self.laser_efficiency = laser_efficiency
        self.power_fraction = power_fraction

        self.targ_theta = targ_theta

        self.orbit_altitude = orbit_altitude            # Orbit altitude in km

        self.laser_heat_capacity = laser_mass * ThermalAnalysisConstants.aluminium_heat_capacity
        self.radiator_heat_capacity = radiator_mass * ThermalAnalysisConstants.aluminium_heat_capacity
        self.number_of_radiators = number_of_radiators


def get_solar_heat(plot_result=False):
    # calculate surface projection of incident sunlight onto panel
    phi_tilt = ta_struct.panel_inclination
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
    phi_tilt = ta_struct.panel_inclination
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
    p_las[np.logical_and(theta > min_power_theta, theta < max_power_theta)] = ta_struct.p_las * (1.0 - ta_struct.laser_efficiency)

    return p_las


if __name__ == '__main__':
    # Generate periodic domain and variable struct
    num_pts = 1001
    theta = np.linspace(0.0, 2 * np.pi, num_pts)
    ta_struct = ThermalAnalysisConstants()
    p_sol = get_solar_heat(plot_result=False)
    p_lun = get_lunar_heat(plot_result=False)
    p_las = get_laser_power()
    p_sol_interp = interp1d(theta, p_sol)
    p_lun_interp = interp1d(theta, p_lun)
    p_las_interp = interp1d(theta, p_las)

    # Get laser heating period for heat switch in integrator
    targ_theta = ta_struct.targ_theta
    power_fraction = ta_struct.power_fraction
    min_power_theta = targ_theta - np.pi * power_fraction
    max_power_theta = targ_theta + np.pi * power_fraction
    assert min_power_theta > 0.0
    assert max_power_theta < 2 * np.pi

    # Set up full simulation time and additional structures
    n_periods = 3
    period = 250.0 * 60.0
    t_end = n_periods * period
    t_sim = 0.0
    t = np.linspace(0.0, t_end, num_pts)
    angle_to_time = period / (2 * np.pi)
    phi_tilt = ta_struct.panel_inclination
    T_las_final = np.zeros(t.shape)
    T_rad_final = np.zeros(t.shape)
    heater_on = False
    laser_coupled=False
    T_las= 273.0
    T_rad = 238.0
    start_idx = 0
    num_iter = 0
    start_idx_glob_prev = 0
    break_points = list()
    while t_sim < t_end:
        start_idx_glob_prev += start_idx

        def integration_model(y, time, laser_coupled=False, heater_on=False):
            T_las = y[0]
            T_rad = y[1]
            local_theta = (time / angle_to_time) % (2 * np.pi)

            # Calculate laser change in temperature
            if (local_theta > min_power_theta and local_theta < max_power_theta):
                p_las = ta_struct.p_las * (1.0 - ta_struct.laser_efficiency)
            else:
                p_las = 0.0

            if (0.5 * (T_las + T_rad) > 273.0):
                p_las2rad = (T_las - T_rad) * 60.0
            else:
                p_las2rad = 0.0
            
            p_tot = p_las - ta_struct.number_of_radiators * p_las2rad
            dTlas_dt = p_tot / ta_struct.laser_heat_capacity
            
            p_sol = p_sol_interp(local_theta)
            p_lun = p_lun_interp(local_theta)
            if (heater_on and T_rad <= 233.0):
                p_heater = 2500.0
            elif (heater_on and T_rad > 233.0):
                p_heater = 0.0
            elif (not heater_on and T_rad > 223.0):
                p_heater = 0.0
            elif (not heater_on and T_rad <= 223.0):
                p_heater = 2500.0
            else:
                print(T_rad, heater_on)

            p_rad = (ta_struct.front_emissivity + ta_struct.back_emissivity) * ta_struct.stefan_boltzmann * T_rad ** 4 * ta_struct.panel_area
            p_tot_rad = p_sol + p_lun + p_heater + p_las2rad - p_rad
            
            dTRad_dt = p_tot_rad / ta_struct.radiator_heat_capacity

            return [dTlas_dt, dTRad_dt]
        if start_idx_glob_prev == 0:
            T = integrate.odeint(integration_model, [T_las, T_rad], t[start_idx_glob_prev + 1:], args=(laser_coupled, heater_on), hmin=1.0e-6, hmax=1000.0)
        else:
            T = integrate.odeint(integration_model, [T_las, T_rad], t[start_idx_glob_prev + 1:], args=(laser_coupled, heater_on), hmin=1.0e-6, hmax=1000.0)
        T_las = T[:, 0]
        T_rad = T[:, 1]

        # Test if simulation has finished
        dT_las = np.abs(np.diff(T_las))
        dT_rad = np.abs(np.diff(T_rad))
        las_zero_indices = np.where(dT_las > 25.0)[0]
        rad_zero_indices = np.where(dT_rad > 25.0)[0]
        if len(rad_zero_indices) == 0:
            T_las_final[start_idx_glob_prev+1:] = T_las
            T_rad_final[start_idx_glob_prev+1:] = T_rad
            break

        # Plot intermediate integrations, removing large number from results outside of integration
        start_idx = min(las_zero_indices[0], rad_zero_indices[0])
        T_las[start_idx+1:] = 0.0
        T_rad[start_idx+1:] = 0.0
        
        # Find new start index and store results from integration
        t_sim = t[start_idx_glob_prev + start_idx + 1]
        T_las_final[start_idx_glob_prev:start_idx_glob_prev + start_idx+1] = T_las[:start_idx+1]
        T_rad_final[start_idx_glob_prev:start_idx_glob_prev + start_idx+1] = T_rad[:start_idx+1]

        # Determine switch that has been triggered
        T_las = T[start_idx, 0]
        T_rad = T[start_idx, 1]
        assert 400.0 > T_las >= 222.99, T_las
        assert 400.0 > T_rad >= 222.99, T_rad
        T_averaged = 0.5 * (T_rad + T_las)
        if np.isclose(T_rad, 223.0, atol=0.25):
            print("Heater on", T_las, T_rad, T_averaged, t_sim)
            heater_on = True
        elif np.isclose(T_rad, 233.0, atol=0.25):
            print("Heater off", T_las, T_rad, T_averaged, t_sim)
            heater_on = False
        elif np.isclose(T_averaged, 273.0, atol=0.5):
            if T_averaged > 273.0:
                print("Laser de-coupled", T_las, T_rad, T_averaged, t_sim)
                laser_coupled=False
            else:
                print("Laser Coupled", T_las, T_rad, T_averaged, t_sim)
                laser_coupled=True
        else:
            # This is triggered by discontinuities in the solar and lunar irradiance
            pass

        # store break points
        break_points.append(t_sim)

    print(len(break_points))

    t /= 60.0
    fig, ax = plt.subplots()
    ax.plot(t, T_las_final - 273.0, c='b')
    ax.plot(t, T_rad_final - 273.0, c='r')
    ax.plot(t, 0.5 * (T_rad_final + T_las_final) - 273.0, linestyle='--', c='grey')
    for t_bp in break_points:
        ax.axvline(t_bp / 60.0, linestyle=':', c='grey', linewidth=0.5)
    ax.grid()
    ax.set_xlim([0.0, t[-1]])
    ax2 = ax.twinx()
    ax2.plot(t, p_las_interp((t * 60.0 / angle_to_time) % (2 * np.pi)))
    ax2.plot(t, p_lun_interp((t * 60.0 / angle_to_time) % (2 * np.pi)))
    ax2.plot(t, p_sol_interp((t * 60.0 / angle_to_time) % (2 * np.pi)))
    plt.show()

