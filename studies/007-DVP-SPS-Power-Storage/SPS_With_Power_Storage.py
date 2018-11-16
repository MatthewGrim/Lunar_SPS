""""
29/05/2018
Author: Darian van Paridon

This script evaluates the times for which the SPS would be required to use stored power are calculated (defined as
events when the SPS in in range, while the satellite and target are simultaneously eclipsed). The necessary battery
size is also calculated, as a function of transmitter power and maximum 'stored power' event duration. Data is available
for equatorial high altitude orbit, and polar lower altitude orbit.

"""

from Lunar_SPS.pysrc.STK_functions.DVP_Programmatic_Functions import *
from Lunar_SPS.pysrc.post_process_functions.DVP_general_SPS_functions import *


def main():

    # Select which data set to access:
    study_name = 'Equatorial_IncrementedRes'
    # study_name = 'Equatorial_IncrementedRes'

    # Set bounds on parametric scan
    max_perigee = 5000.0
    max_apogee = 5000.0

    # Get orbit data set
    semi_maj_axis, eccentricity, orbit_data = vary_orbital_elements_incrementing_resolution(max_perigee, max_apogee)

    # Initialize simulation times, and refer to altitude corresponding to orbit selection
    start = convert_string_to_datetime(['2018', '05', '17', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '17', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    main_directory = os.path.dirname(os.path.dirname(current_folder))

    # Set file path for data
    stk_data_path = os.path.join(main_directory, 'STK Data', study_name)

    max_stored_power_time = sort_incremented_resolution_data(orbit_data, read_data_from_file(stk_data_path, study_name, 'MaxStoredPowerEvent'))
    total_stored_power_time = sort_incremented_resolution_data(orbit_data, read_data_from_file(stk_data_path, study_name, 'TotalStoredPowerEvent'))

    # Estimate battery size
    trans_power = 10e3
    trans_eff = 0.35
    li_ion_energy_density = 140.0  # Watt hours per kilogram
    fuel_cell_energy_density = 500.0
    battery_mass = [trans_power * i / (trans_eff * 3600.0 * li_ion_energy_density) for i in max_stored_power_time]
    fuel_cell_mass = [trans_power * i / (trans_eff * 3600.0 * fuel_cell_energy_density) for i in max_stored_power_time]

    # Extract unique perigees and apogees tested for plotting
    unique_perigees = [orbit_data[1][0]]
    unique_apogees = [orbit_data[1][1]]
    r_moon = 1737.0
    for i in range(1, len(orbit_data)):
        if orbit_data[i][0] > max(unique_perigees):
            unique_perigees.append(orbit_data[i][0])
        if orbit_data[i][1] > max(unique_apogees):
            unique_apogees.append(orbit_data[i][1])

    PERI, APO = np.meshgrid(np.asarray(unique_perigees), np.asarray(unique_apogees), indexing='ij')
    SMA = (PERI + APO) * 1e3 / 2
    mu_moon = 6.674e-11 * 7.347673e22
    mask = PERI > APO
    # Charge time of the battery per orbit period
    orbit_period = 2 * np.pi * np.sqrt(SMA ** 3 / mu_moon) / 3600.0
    charge_time = orbit_period / 2
    number_of_cycles = 24 / orbit_period * 365 * 10
    number_of_cycles[mask] = np.nan

    # Reduce perigee and apogee to altitudes instead of radii
    perigee_altitudes = [i - r_moon for i in unique_perigees]
    apogee_altitudes = [i - r_moon for i in unique_apogees]

    fig, ax = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(7, 7))
    im = ax[0, 0].contourf(apogee_altitudes, perigee_altitudes, total_stored_power_time / (3600.0 * 24 * 365 * 2), 500)
    ax[0, 0].set_ylabel('Perigee Altitude [km]')
    fig.colorbar(im, ax=ax[0, 0])
    ax[0, 0].set_title('Total Beaming Time [%]')
    
    im = ax[1, 0].contourf(apogee_altitudes, perigee_altitudes, np.asarray(battery_mass) * li_ion_energy_density * 1e-3, 500)
    fig.colorbar(im, ax=ax[1, 0])
    ax[1, 0].set_title('Energy Capacity [kWhrs]')
    ax[1, 0].set_ylabel('Perigee Altitude [km]')
    ax[1, 0].set_xlabel('Apogee Altitude [km]')
    
    im = ax[0, 1].contourf(apogee_altitudes, perigee_altitudes, number_of_cycles, 500)
    fig.colorbar(im, ax=ax[0, 1])
    ax[0, 1].set_title('$N_{cycles}$ [over 10 years]')

    im = ax[1, 1].contourf(apogee_altitudes, perigee_altitudes, battery_mass, 500)
    fig.colorbar(im, ax=ax[1, 1])
    ax[1, 1].set_title('LiPo Battery Mass [kg]')
    ax[1, 1].set_xlabel('Apogee Altitude [km]')
    fig.suptitle("Battery mass requirements for a {} kW laser".format(trans_power * 1e-3))
    plt.savefig('battery_mass_requirements')
    plt.show()
    
    plt.figure()
    plt.title('Charge time per orbit period [hrs]')
    plt.contourf(apogee_altitudes, perigee_altitudes, charge_time, 500)
    plt.xlabel('Apogee Altitude [km]')
    plt.ylabel('Perigee Altitude [km]')
    plt.colorbar()
    plt.savefig('available_charge_time')
    plt.show()
    

main()
