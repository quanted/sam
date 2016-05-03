import os
import csv
import re
from collections import defaultdict

def map_dir(d):
    out_map = defaultdict(dict)
    for f in os.listdir(d):
        match = re.match(".+?_(\d+?)_(\d{4})", f)
        if match:
            reach_id, year = match.groups()
            out_map[reach_id][year] = os.path.join(d, f)
    return out_map

def read_file(path):
    contents = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            date = row['Date']
            contents[date] = row
    return contents

def compare(dict1, dict2):
    for date in sorted(dict1.keys()):
        if date in dict2:
            row1, row2 = dict1[date], dict2[date]
            for field, val1 in row1.items():
                val2 = row2.get(field)
                if val2:
                    if val1 != val2:
                        print("{}, {}: With Erosion: {}, Without: {}".format(date, field, val1, val2))
                else:
                    print("Field {} not in both".format(field))


def compare_folders(with_dir, without_dir):
    with_map = map_dir(with_dir)
    without_map = map_dir(without_dir)
    for reach_id in without_map:
        for year, without_path in sorted(without_map[reach_id].items()):
            with_path = with_map[reach_id][year]
            if with_path and without_path:
                print(reach_id, year)
                with_data = read_file(with_path)
                without_data = read_file(without_path)
                compare(with_data, without_data)
            else:
                print("{}/{} not found for both datasets".format(reach_id, year))




def main():
    with_dir = r"T:\pySAM\bin\Outputs\Python\ErosionTest\With"
    without_dir = r"T:\pySAM\bin\Outputs\Python\ErosionTest\Without"

    compare_folders(with_dir, without_dir)

if __name__ == "__main__":
    main()