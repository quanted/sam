import numpy as np
import pickle
import os
import datetime

def crush_cube(materials, outfile):
    if not os.path.isdir(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))
    with open(outfile, 'wb') as f:
        pickle.dump(materials, f)

def read_template(template_file):
    with open(template_file, 'rb') as f:
        _, map_dict, dates = pickle.load(f)
    return map_dict, dates

def create_dummy_output(out_file, template, n_dates, events, baseflow):

    map_dict, dates = read_template(template)

    dates = [min(dates) + datetime.timedelta(days=i) for i in range(n_dates)]

    n_reaches = len(map_dict.keys())

    output_array = np.zeros((3, n_reaches, n_dates))

    for day, mass, runoff in events:

        output_array[0, :, day] = mass
        output_array[1, :, day] = runoff

    for start_day, flow in baseflow:
        output_array[2, :, start_day:] = flow

    crush_cube((output_array, map_dict, dates), out_file)


def main():

    out_file = r"T:\SAM\Preprocessed\OutputCubes\cube_dummy_mtb.p"

    # Template is used to get reach names and map
    template = r"T:\SAM\Preprocessed\OutputCubes\cube_mark_twain.p"

    n_dates = 100

    # Each event is a tuple with the form (day, mass, runoff)
    events = [( 5,  1,  2),
              (15,  2,  1),
              (30, 10,  5),
              (50,  5, 10),
              (80,  1,  0),
              (90,  0,  1)]

    # Each baseflow entry is a tuple with the form (start day, flow), and carries to the next start day
    baseflow = [( 0, 1),
                (10, 2),
                (20, 3),
                (30, 4),
                (40, 5),
                (60, 4),
                (70, 3),
                (80, 2),
                (90, 1)]
    create_dummy_output(out_file, template, n_dates, events, baseflow)


if __name__ == "__main__":
    main()