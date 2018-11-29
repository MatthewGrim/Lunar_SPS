"""
Author: Rohan Ramasamy
Date: 01/08/2018

This script generates the steady state temperature as a function of total sps efficiency
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


def get_steady_state_temperature():
    eta = np.linspace(0.0, 1.0, 1000)
    epsilon = 1.0
    sigma = 5.670367e-8
    phi = 1360.0

    T = (phi * (1 - eta) / (2 * epsilon * sigma)) ** 0.25
    T_interp = interpolate.interp1d(eta, T)
    eta_estimate = 0.35 * 0.30
    print(eta_estimate)
    print(T_interp(eta_estimate) - 273)

    plt.figure(figsize=(8, 8))

    plt.plot(eta * 100.0, T)
    plt.axvline(eta_estimate * 100.0, linestyle="--")

    plt.ylabel("Steady State Temperature [K]")
    plt.xlabel("$\eta_{pv} \eta_{laser}$ [%]")
    plt.xlim([0.0, 100.0])
    plt.ylim([np.min(T), np.max(T)])
    plt.tight_layout()
    plt.savefig("steady_state_temperature")
    plt.show()

if __name__ == '__main__':
    get_steady_state_temperature()

