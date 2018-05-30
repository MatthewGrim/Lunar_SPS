""""
29/05/2018
Author: Darian van Paridon

This script evaluates the total active and blackout times for an SPS configuration. In addition, the times for
which the SPS would be required to use stored power are calculated. The resulting reduction in total and max
blackout durations are printed to evaluate the benefit of power storage onboard the SPS. The necessary battery size is
also calculated, as a function of transmitter power and duty cycle while in eclipse.
"""

from DVP_general_SPS_functions import *


def determine_battery_size(trans_power, duty_cycle, longest_event):

    li_ion_energy_density = 270  # Watt hours per kilogram
    duration = duty_cycle * longest_event
    battery_capacity = trans_power * duration  # in Watt hours
    battery_size = battery_capacity / li_ion_energy_density
    print("To power the SPS for {}% of this duration, a {} kg Li-ion battery is required".format(active_percentage*100.0, round(battery_size, 2)))
    return battery_size


def main():
    start = convert_string_to_datetime(['2018', '05', '18', '10', '0', '0.0'])
    end = convert_string_to_datetime(['2020', '05', '18', '10', '0', '0.0'])
    total_duration = (end - start).total_seconds()

    # Import SPS illumination events
    sps_lighting_raw = 'SPS(3.5k)_Lighting.csv'
    sps_lighting = parse_csv_to_array(sps_lighting_raw, start)
    sps_eclipse = invert_events_list(sps_lighting, total_duration)

    # Import access data for sps
    raw_access = "SPS(3.5k)_Access.csv"
    sps_access = parse_csv_to_array(raw_access, start)

    # Import target lighting data, and invert to find eclipse events
    raw_target_lighting = "Target_Lighting.csv"
    target_lighting = parse_csv_to_array(raw_target_lighting, start)
    target_eclipse = invert_events_list(target_lighting, total_duration)

    print('Maximum black-out duration w/o SPS: {} hrs'.format(round(max(target_eclipse[2]) / 3600.0, 2)))
    print('Total time spent in black-out w/o SPS: {}% per year'.format(
        round(100.0 * np.sum(target_eclipse[2]) / (2.0 * 3600.0 * 365.0 * 24.0), 2)))

    stored_power_events = determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)
    duty_cycle = 0.5
    trans_power = 100e3
    determine_battery_size(trans_power, duty_cycle, max(stored_power_events[2]) / 3600.0)
    determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)



main()
