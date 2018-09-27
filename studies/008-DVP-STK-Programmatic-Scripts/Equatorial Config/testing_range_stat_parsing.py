
import os


def import_range_data_statistics(file_name, stk_data_path):
    import re
    import os

    the_file = open("{}/{}.txt".format(stk_data_path, file_name), "r")
    new_file = open("{}/{}_Semi_Parsed.txt".format(stk_data_path, file_name), "w+")
    mean_range = []
    max_range = []
    min_range = []

    # Remove lines which do not contain data, and sort columns
    # by removing white space
    for i, line in enumerate(the_file):
        if line[0] == "M":
            removed_white_space = re.split('\s{2,}', line.strip())
            new_file.write("{}\n".format(removed_white_space))
    the_file.close()
    new_file.close()

    last_file = open("{}/{}_Semi_Parsed.txt".format(stk_data_path, file_name), "r")
    # Sort data into minimum maximum and mean arrays
    for j, new_line in enumerate(last_file):
        new_line_split = new_line.split(",")
        if new_line_split[0][2:-1] == "Min Range":
            min_range.append(float(new_line_split[2][2:-3]))
        elif new_line_split[0][2:-1] == "Max Range":
            max_range.append(float(new_line_split[2][2:-3]))
        elif new_line_split[0][2:-1] == "Mean Range":
            mean_range.append(float(new_line_split[1][2:-3]))
    new_file.close()

    range_stats = [min_range, max_range, mean_range]

    return range_stats


def main():

    resolution = 1000.0

    # Get pathway to main Lunar_SPS directory
    current_folder = os.getcwd()
    issue_folder = os.path.dirname(current_folder)
    main_directory = os.path.dirname(issue_folder)

    # Name of study
    study_name = 'Brandhorst_{}kmRes'.format(resolution)

    # File path
    stk_data_path = r'{}\STK Data\{}'.format(main_directory, study_name)

    file_name = "SPS_Range_Stats"

    range_stats = import_range_data_statistics(file_name, stk_data_path)

    print(range_stats)

main()
