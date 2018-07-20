
from DVP_Programmatic_Functions import calculate_orbital_perturbations, read_data_from_file, write_data_to_file
from DVP_general_SPS_functions import convert_string_to_datetime
import numpy as np
import matplotlib.pyplot as plt
import os


def main():

    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    duration = (end - start).total_seconds()
    time_step = 3600.0

    r_moon = 1737.0
    perigee_altitude = 10.0 + r_moon
    apogee_altitude = [10.0, 20.0, 50.0, 100.0, 300.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0]
    apogee_altitude = [i + r_moon for i in apogee_altitude]
    inclination = 90.0
    arg_perigee = 90.0

    eccentricity = [(j - perigee_altitude) / (j + perigee_altitude) for j in apogee_altitude]
    semi_maj_axis = [perigee_altitude / (1 - i) for i in eccentricity]

    # Calculate drift rate with analytical model
    # Argument of perigee [0], eccentricity [1], inclination [2]
    change_in_orbital_elements = calculate_orbital_perturbations(semi_maj_axis, eccentricity, inclination, arg_perigee)

    # Read in drift rate from STK data reports
    current_folder = os.getcwd()
    stk_data_path = '{}\STK Data\SouthPole_ArgPeriDriftRate'.format(current_folder)
    norm_diff = []

    # Read in processed STK data for drift rate and total drift of argument of perigee
    stk_data_drift_rate = read_data_from_file(stk_data_path, 'ArgPeri_Stability_Comparison', 'STK_DriftRate_ArgPeri')
    stk_data_total_drift = read_data_from_file(stk_data_path, 'ArgPeri_Stability_Comparison', 'STK_TotalDrift_ArgPeri')

    # Calculate normalized difference between STK data and analytical model
    for i in range(len(stk_data_drift_rate)):
        norm_diff.append(abs(abs(stk_data_drift_rate[i]) - abs(change_in_orbital_elements[0][i] * 180.0 / np.pi)) / abs(stk_data_drift_rate[i]))

    apogee_altitude = [i - r_moon for i in apogee_altitude]

    # Show comparison
    plt.plot(apogee_altitude, change_in_orbital_elements[0] * (180.0 / np.pi), label='Analytic Model')
    for i in range(len(apogee_altitude)):
        plt.scatter(apogee_altitude[i], stk_data_drift_rate[i], marker='x', c='k', label='STK Data')
    plt.legend()
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('$\dot{\omega}$ [deg/s]')
    plt.title('Drift Rate in $\omega$ for 10 km Perigee Altitude')
    plt.show()

    # Normalized difference
    plt.plot(apogee_altitude, norm_diff)
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('$\dot{\omega}$ [deg/s]')
    plt.title('Normalized Difference Between Models')
    plt.show()

    # Total drift
    plt.plot(apogee_altitude, stk_data_total_drift, label='STK Data')
    plt.plot(apogee_altitude, duration * change_in_orbital_elements[0] * (180.0 / np.pi), label='Analytical Model')
    plt.legend()
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('$\Delta \omega$ [deg]')
    plt.title('Total $\Delta \omega$ for 10 km Perigee Altitude')
    plt.show()


def parse_orbital_elements_data(file_name):

    the_file = open(file_name, "r")
    arg_perigee_data = list()

    for i, line in enumerate(the_file):
        if i == 0:
            continue
        # Split line into table components, assuming comma delimited
        components = line.split(",")
        # Break once the final line is reached (which is a blank 'return')
        if components[0] == '\n':
            break
        else:
            arg_perigee_data.append(float(components[5]))

    return arg_perigee_data


main()
