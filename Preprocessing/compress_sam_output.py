import pickle
import os
import numpy as np
import datetime
import re

def map_sam_output(sam_output_dir, sam_output_format, year_filter=2010):

    # Map output files to dictionary
    output_files = {}
    for f in os.listdir(sam_output_dir):
        match = re.match(sam_output_format, f)
        if match:
            reach_id, year = map(int, match.groups())
            if year == year_filter:
                output_files[reach_id] = os.path.join(sam_output_dir, f)
    return output_files



def read_sam_output(output_files):

    # Written for header [Date, Conc(ug/L), RMass(kg), Runoff(m), RConc(ug/L), TotalFlow(m3), baseflow(m3)]
    # Looking for RMass (col 2), Runoff (col 3), baseflow (col 6)
    # sam_output is a 3d array with dimensions [attributes, reaches, dates]
    sam_lookup = {}
    initial = True
    for i, reach in enumerate(output_files.items()):
        reach_id, reach_file = reach
        data = np.genfromtxt(reach_file, delimiter=' ', skip_header=1, usecols=(2,3,6))  # (5479, 3) [ndates, attributes]
        dates = np.genfromtxt(reach_file, dtype=np.str, delimiter=' ', skip_header=1, usecols=(0))  # (5479, 3) [ndates, attributes]
        if initial:
            size = (len(output_files), data.shape[0], data.shape[1])  # (920, 5479, 3)
            sam_output = np.ndarray(size)
            initial = False
        sam_output[i] = data[:]
        sam_lookup[reach_id] = i
    sam_output = np.rollaxis(sam_output, 2)
    dates = [datetime.date(*map(int, date.split("-"))) for date in dates]
    return sam_output, sam_lookup, dates


def crush_cube(materials, collection_id, cube_dir, cube_format):
    if not os.path.isdir(cube_dir):
        os.mkdir(cube_dir)
    outfile = os.path.join(cube_dir, cube_format.format(collection_id))
    with open(outfile, 'wb') as f:
        pickle.dump(materials, f)


def compress_sam_output_into_a_cube(sam_output_dir, sam_output_format, cube_dir, cube_format, collection_id, year_filter):

    output_files = map_sam_output(sam_output_dir, sam_output_format, year_filter)

    materials = read_sam_output(output_files)

    crush_cube(materials, collection_id, cube_dir, cube_format)

    print("Your data has been crushed into a cube\nYou have one hour to remove your cube")

def main():

    sam_output_dir = r"T:\SAM\Outputs\Python"

    sam_output_format = "Eco_(\d+?)_(\d{4})_daily.out"

    cube_dir = r"T:\SAM\Preprocessed\OutputCubes"

    cube_format = r"cube_{}.p"

    collection_id = "mark_twain"

    year_filter = 2010

    compress_sam_output_into_a_cube(sam_output_dir, sam_output_format, cube_dir, cube_format, collection_id, year_filter)


if __name__ == "__main__":
    main()
