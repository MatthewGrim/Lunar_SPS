"""
Author: Rohan Ramasamy
Date: 10/04/2020

This file contains a thermal assessment of a rover powered by laser light. The calculations are based on 
the original work from study 017.
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

    def __init__(self, panel_area=1.83, 
                 front_emissivity=0.84, front_absorptivity=0.5, 
                 back_emissivity=0.12,
                 panel_inclination=0.0,
                 p_las=248.0,
                 targ_theta = 2.0 * np.pi / 2.0,
                 T_service=0.1,
                 panel_density=3.5,
                 panel_heat_capacity=500.0,
                 p_heater=26.7):
        """
        Initialise design based constants
        """
        # Design based
        self.panel_area = panel_area                    # Area of individual radiator panels
        self.front_emissivity = front_emissivity        # Emissivity of space facing side
        self.front_absorptivity = front_absorptivity    # Emissivity of space facing side
        self.back_emissivity = back_emissivity          # Emissivity of non_space facing side
        self.panel_inclination = panel_inclination      # Tilt of panel from the lunar equatorial plane

        self.p_las = p_las
        self.T_service = T_service
        
        self.targ_theta = targ_theta
        self.panel_heat_capacity = panel_heat_capacity
        self.panel_density = panel_density
        self.p_heater = p_heater

def get_laser_power(ta_struct, theta, plot_result=False):
    # Get theta angles when laser is beaming power
    targ_theta = ta_struct.targ_theta
    T_service = ta_struct.T_service
    min_power_theta = targ_theta - np.pi * T_service
    max_power_theta = targ_theta + np.pi * T_service
    assert min_power_theta > 0.0
    assert max_power_theta < 2 * np.pi

    p_las = np.zeros(theta.shape)
    p_las[np.logical_and(theta > min_power_theta, theta < max_power_theta)] = ta_struct.p_las * ta_struct.front_absorptivity

    return p_las

def run_temperature_analysis(p_las=250, p_heater=26.7, n_periods=3, num_pts=10001, plot_result=False):
    # Generate periodic domain and variable struct
    theta = np.linspace(0.0, 2 * np.pi, n_periods * num_pts)
    ta_struct = ThermalAnalysisConstants(p_las=p_las, p_heater=p_heater)
    p_las = get_laser_power(ta_struct, theta)
    p_las_interp = interp1d(theta, p_las)

    # Get laser heating period for heat switch in integrator
    targ_theta = ta_struct.targ_theta
    T_service = ta_struct.T_service
    min_power_theta = targ_theta - np.pi * T_service
    max_power_theta = targ_theta + np.pi * T_service
    assert min_power_theta > 0.0
    assert max_power_theta < 2 * np.pi

    # Set up full simulation time and additional structures
    period = 250.0 * 60.0
    t_end = n_periods * period
    t = np.linspace(0.0, t_end, n_periods * num_pts)
    dt = t[1] - t[0]
    angle_to_time = period / (2 * np.pi)
    
    phi_tilt = ta_struct.panel_inclination
    T_panel = np.zeros(t.shape)
    P_heater = np.zeros(t.shape)
    kwargs = {
        "heater_on" : False,
        "time" :  0.0
    }
    T_p= 273.3 - 146.0
    T_panel[0] = T_p

    def integration_model(T_p, **kwargs):
        heater_on = kwargs.get("heater_on")
        time = kwargs.get("time")

        local_theta = (time / angle_to_time) % (2 * np.pi)

        # Calculate laser change in temperature
        if (local_theta > min_power_theta and local_theta < max_power_theta):
            p_las = ta_struct.p_las * ta_struct.front_absorptivity * ta_struct.panel_area
        else:
            p_las = 0.0

        dTpanel_dt = p_las / (ta_struct.panel_heat_capacity * ta_struct.panel_density * ta_struct.panel_area)
        
        if (heater_on):
            p_heater = ta_struct.p_heater
        else:
            p_heater = 0.0

        p_rad = (ta_struct.front_emissivity + ta_struct.back_emissivity) * ta_struct.stefan_boltzmann * T_p[0] ** 4 * ta_struct.panel_area
        p_tot_rad = p_rad - p_heater
        
        dTpanel_dt -= p_tot_rad / (ta_struct.panel_heat_capacity * ta_struct.panel_density * ta_struct.panel_area)

        return dTpanel_dt

    for t_idx, t_curr in enumerate(t):
        if t_idx == 0:
            continue
        
        kwargs["time"] = t_curr
        T_p = rKN([T_p], [integration_model], 1, dt, **kwargs)[0]
        T_panel[t_idx] = T_p
        
        # Determine switch that has been triggered
        assert 400.0 > T_p >= 119.3, T_p
        if (T_p <= 127.3):
            print("Heater on", T_p, t_curr)
            kwargs["heater_on"] = True
        if (T_p > 133.3):
            print("Heater off", T_p, t_curr)
            kwargs["heater_on"] = False

        if kwargs["heater_on"]:
            P_heater[t_idx] = ta_struct.p_heater

    P_las = p_las_interp((t / angle_to_time) % (2 * np.pi))
    if plot_result:
        t /= 60.0
        fig, ax = plt.subplots()
        ax.plot(t, T_panel - 273.0, c='b')
        ax.grid()
        ax.set_xlim([0.0, t[-1]])
        ax2 = ax.twinx()
        ax2.plot(t, P_las, c='g')
        ax2.plot(t, P_heater, c='r')
        plt.show()
    
    return t, T_panel, P_las, P_heater

if __name__ == '__main__':
    fig, ax = plt.subplots(2, figsize=(14, 7), sharex=True)

    # Average temperature case
    t, T_panel, P_las, P_heater = run_temperature_analysis()

    t /= 3600
    lns1 = ax[0].plot(t, T_panel - 273.0, c='b', label='$T_{panel}$')
    ax[0].grid(linestyle='--')
    ax[0].set_xlim([0.0, t[-1]])
    ax[0].set_ylim([-200, 50])
    ax[0].set_ylabel('$T\ [^\circ C]$', fontsize=14)
    ax2 = ax[0].twinx()
    lns2 = ax2.plot(t, P_las, c='g', label='$P_{laser}$')
    lns3 = ax2.plot(t, P_heater, c='r', label='$P_{heater}$')
    ax2.set_ylim([0, 450])
    
    ax2.set_ylabel('$P\ [W]$', fontsize=14)
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    ax[0].legend(lns, labs, prop={'size': 14})

    # Peak temperature case
    t, T_panel, P_las, P_heater = run_temperature_analysis(p_las=800.0)

    t /= 3600
    lns1 = ax[1].plot(t, T_panel - 273.0, c='b', label='$T_{panel}$')
    ax[1].grid(linestyle='--')
    ax[1].set_xlim([0.0, t[-1]])
    ax[1].set_ylim([-200, 50])
    ax[1].set_ylabel('$T\ [^\circ C]$', fontsize=14)
    ax2 = ax[1].twinx()
    lns2 = ax2.plot(t, P_las, c='g', label='$P_{laser}$')
    lns3 = ax2.plot(t, P_heater, c='r', label='$P_{heater}$')
    ax2.set_ylim([0, 450])
    ax2.set_ylabel('$P\ [W]$', fontsize=14)
    ax[1].set_xlabel('$Time\ [hrs]$', fontsize=14)

    fig.tight_layout()
    plt.show()