

from DVP_general_SPS_functions import *


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

    determine_SPS_storedpower_time(sps_eclipse, target_eclipse, sps_access)
    determine_SPS_active_time(sps_lighting, target_eclipse, sps_access)


main()
