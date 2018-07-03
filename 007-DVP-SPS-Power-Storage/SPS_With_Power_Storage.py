""""
29/05/2018
Author: Darian van Paridon

This script evaluates the times for which the SPS would be required to use stored power are calculated (defined as
events when the SPS in in range, while the satellite and target are simultaneously eclipsed). The necessary battery
size is also calculated, as a function of transmitter power and maximum 'stored power' event duration. Data is available
for equatorial high altitude orbit, and polar lower altitude orbit.

"""

from DVP_general_SPS_functions import *


def main():

    # Select orbit type which is studied. Options are Equatorial or Polar
    orbit = 'Polar'

    # Initialize simulation times, and refer to altitude corresponding to orbit selection
    if orbit == 'Equatorial':
        altitude = 30
        start = convert_string_to_datetime(['2008', '07', '01', '10', '0', '0.0'])
        end = convert_string_to_datetime(['2010', '06', '30', '10', '0', '0.0'])
        total_duration = (end - start).total_seconds()

    elif orbit == 'Polar':
        altitude = 3.5
        start = convert_string_to_datetime(['2018', '05', '18', '10', '0', '0.0'])
        end = convert_string_to_datetime(['2020', '05', '18', '10', '0', '0.0'])
        total_duration = (end - start).total_seconds()

    # Import SPS illumination events
    sps_lighting_raw = '{}/SPS({}k)_Lighting.csv'.format(orbit, altitude)
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    sps_eclipse = invert_events_list(sps_lighting, total_duration)

    # Import access data for sps
    raw_access = "{}/SPS({}k)_Access.csv".format(orbit, altitude)
    sps_access = parse_csv_to_array(raw_access, start)

    # Import target lighting data, and invert to find eclipse events
    raw_target_lighting = "{}/Target_Lighting.csv".format(orbit)
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    print('Total time spent in black-out w/o SPS: {}% per year'.format(
        round(100.0 * np.sum(target_eclipse[2]) / total_duration, 2)))

    stored_power_events = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)

    print('Total time available for using stored power for WPT: {}% per year'.format(round(100.0 * np.sum(stored_power_events[2]) / total_duration, 2)))

    # Estimate battery size
    trans_power = 100e3
    li_ion_energy_density = 270.0  # Watt hours per kilogram
    battery_capacity = trans_power * max(stored_power_events[2]) / 3600.0  # in Watt hours
    battery_mass = battery_capacity / li_ion_energy_density

    print("To power a {} kW transmitter for the longest event where store power is required ({} hours), a {} kg Li-ion battery is required".format(round(trans_power / 1000.0, 2), round(max(stored_power_events[2]) / 3600.0, 2), round(battery_mass, 2)))


main()
