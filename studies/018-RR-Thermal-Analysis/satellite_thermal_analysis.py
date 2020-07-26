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


def rKN(x, fx, n, hs, **kwargs):
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    xk = []
    for i in range(n):
        k1.append(fx[i](x, **kwargs)*hs)
    for i in range(n):
        xk.append(x[i] + k1[i]*0.5)
    for i in range(n):
        k2.append(fx[i](xk, **kwargs)*hs)
    for i in range(n):
        xk[i] = x[i] + k2[i]*0.5
    for i in range(n):
        k3.append(fx[i](xk, **kwargs)*hs)
    for i in range(n):
        xk[i] = x[i] + k3[i]
    for i in range(n):
        k4.append(fx[i](xk, **kwargs)*hs)
    for i in range(n):
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6
    return x


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
    num_pts = 10001
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

    def laser_model(y, **kwargs):
        time = kwargs.get("time")
        
        T_las = y[0]
        T_rad = y[1]
        local_theta = (time / angle_to_time) % (2 * np.pi)

        # Calculate laser change in temperature
        if (local_theta > min_power_theta and local_theta < max_power_theta):
            p_las = ta_struct.p_las * (1.0 - ta_struct.laser_efficiency)
        else:
            p_las = 0.0

        if (0.5 * (T_las + T_rad) > 273.3):
            p_las2rad = (T_las - T_rad) * 60.0
        else:
            p_las2rad = 0.0
        
        p_tot = p_las - ta_struct.number_of_radiators * p_las2rad
        dTlas_dt = p_tot / ta_struct.laser_heat_capacity

        return dTlas_dt

    def radiator_model(y, **kwargs):
        time = kwargs.get("time")
        heater_on = kwargs.get("heater_on")
        
        T_las = y[0]
        T_rad = y[1]
        local_theta = (time / angle_to_time) % (2 * np.pi)

        if (0.5 * (T_las + T_rad) > 273.3):
            p_las2rad = (T_las - T_rad) * 60.0
        else:
            p_las2rad = 0.0
        
        p_sol = p_sol_interp(local_theta)
        p_lun = p_lun_interp(local_theta)
        if (heater_on):
            p_heater = 2500.0
        else:
            p_heater = 0.0

        p_rad = (ta_struct.front_emissivity + ta_struct.back_emissivity) * ta_struct.stefan_boltzmann * T_rad ** 4 * ta_struct.panel_area
        p_tot_rad = p_sol + p_lun + p_heater + p_las2rad - p_rad
        
        dTRad_dt = p_tot_rad / ta_struct.radiator_heat_capacity

        return dTRad_dt

    # Set up full simulation time and additional structures
    n_periods = 3
    period = 250.0 * 60.0
    t_end = n_periods * period
    t_sim = 0.0
    t = np.linspace(0.0, t_end, n_periods * num_pts)
    dt = t[1] - t[0]
    angle_to_time = period / (2 * np.pi)
    
    phi_tilt = ta_struct.panel_inclination
    T_l= 273.3
    T_r = 238.0
    T_las = np.zeros(t.shape)
    T_rad = np.zeros(t.shape)
    p_heater = np.zeros(t.shape)
    T_las[0] = T_l
    T_rad[0] = T_r
    kwargs = {
        "time" : 0.0,
        "heater_on" : False
    }

    for t_idx, t_curr in enumerate(t):
        if t_idx == 0:
            continue
        
        kwargs["time"] = t_curr
        T = rKN([T_l, T_r], [laser_model, radiator_model], 2, dt, **kwargs)
        
        T_l = T[0]
        T_r = T[1]
        T_las[t_idx] = T_l
        T_rad[t_idx] = T_r
        
        # Determine switch that has been triggered
        # assert 400.0 > T_l >= 222.99, T_l
        # assert 400.0 > T_r >= 222.99, T_r
        T_averaged = 0.5 * (T_r + T_l)
        if T_r < 223.0:
            print("Heater on", T_l, T_r, T_averaged, t_sim)
            kwargs["heater_on"] = True
        if T_r > 233.0:
            print("Heater off", T_l, T_r, T_averaged, t_sim)
            kwargs["heater_on"] = False

        if kwargs["heater_on"]:
            p_heater[t_idx] = 2500.0


    t /= 3600.0
    fig, ax = plt.subplots(figsize=(14, 7))
    lns1 = ax.plot(t, T_las - 273.3, c='r', label='$T_{laser}$')
    lns2 = ax.plot(t, T_rad - 273.3, c='b', label='$T_{radiator}$')
    lns3 = ax.plot(t, 0.5 * (T_rad + T_las) - 273.3, linestyle='--', c='k', label='$T_{average}$')
    ax.grid(linestyle='--')
    ax.legend()
    ax.set_xlim([6.5, 10.8])
    ax.set_ylabel('$T\ [^\circ C]$', fontsize=14)
    ax.set_xlabel('$Time\ [hrs]$', fontsize=14)

    ax2 = ax.twinx()
    lns4 = ax2.plot(t, p_las_interp((t * 3600.0 / angle_to_time) % (2 * np.pi)), label='$P_{laser}$', c='g')
    lns5 = ax2.plot(t, p_lun_interp((t * 3600.0 / angle_to_time) % (2 * np.pi)), label='$P_{lunar}$', c='grey')
    lns6 = ax2.plot(t, p_sol_interp((t * 3600.0 / angle_to_time) % (2 * np.pi)), label='$P_{solar}$', c='orange')
    lns7 = ax2.plot(t, p_heater, label='$P_{heater}$', c='c')
    
    ax2.set_ylabel('$P\ [W]$', fontsize=14)

    lns = lns1+lns2+lns3+lns4+lns5+lns6+lns7
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, prop={'size': 14})
    
    fig.tight_layout()
    plt.show()

