import os
import numpy as np
import travel_time_functions as functions
import datetime

import read
import output

# @@@ - to do:
# class for indexed array?  class for input file with dir and format?  OR structured arrays?
# how to set up?  will SAM be run for whole regions?
# At the moment, reaches with stream_calc = 0 are not processed or included as upstream. Probably should be...

def time_of_travel(lake_file, upstream_file, sam_output_file, output_dir, output_format):

    # Read in compressed SAM output
    sam_output, sam_lookup, start_date = read.unpickle(sam_output_file)

    # Read lake file
    volume_dict, outflow_dict, outlet_dict, residence_times, waterbody_dict = read.unpickle(lake_file)

    # Get upstream path data
    upstream, path_map = read.unpickle(upstream_file)

    # Separate flowing reaches from reaches that are in a reservoir
    lotic_reaches, upstream_lentics = functions.classify_reaches(sam_lookup.keys(), outlet_dict, waterbody_dict)

    # Identify the dates that correspond to the time series
    time_period = sam_output.shape[2]
    dates = [datetime.date(*map(int, start_date.split("-"))) + datetime.timedelta(days=i) for i in range(time_period)]

    # Loop through all reaches with SAM output and run convolution
    for i, reach in enumerate(lotic_reaches):

        reach_address = path_map.get(reach)

        if reach_address:

            # Pull upstream path data for this reach
            upstream_paths, upstream_times, upstream_lakes, upstream_reaches = \
                functions.trim_to_reach(upstream, *reach_address)

            all_upstream = set(np.unique(upstream_paths))

            # Check to make sure that the reach was found in the path matrix.
            if all_upstream:

                # Get all applicable SAM output time series for this reach
                upstream_output, upstream_baseflow, upstream_lookup = \
                    functions.trim_to_upstream(sam_output, sam_lookup, all_upstream)

                # Convolve paths that pass through reservoirs
                upstream_output = functions.convolve_reservoirs(upstream_output, upstream_lookup, upstream_paths,
                                                                upstream_lakes, residence_times, time_period)

                # Snap the travel times to an interval (portions of a day).  Default is 1-day
                upstream_output = functions.convolve_flowing(upstream_reaches, upstream_output, upstream_lookup,
                                                             time_period)

                # Unpack time series and compute concentration
                runoff_mass, total_runoff, baseflow, total_flow, total_conc, runoff_conc = \
                    functions.compute_concentration(upstream_output, upstream_baseflow)

                if runoff_mass.sum():
                    print(reach)

                # If the reach is an outlet for a lake, write the output to the reach and all reaches in the lake
                for reach in {reach} | upstream_lentics[reach]:
                    output.daily(output_dir, output_format, reach, total_conc, runoff_conc, runoff_mass, dates,
                                 total_flow, baseflow, total_runoff)



            else:
                print("{} not found in path map".format(reach))

def main():
    region = '07'
    sam_output_id = "mark_twain"

    lakefile_dir = r"T:\SAM\Preprocessed\LakeFiles"
    upstream_repository = r"T:\SAM\Preprocessed\UpstreamPaths"
    sam_output_dir = r"T:\SAM\Preprocessed\OutputCubes"
    convolution_output_dir = r"T:\SAM\Outputs\Convolved"

    lakefile_format = "region_{}.p".format(region)
    upstream_format = "upstream_{}.p".format(region)
    sam_output_format = "cube_{}.p".format(sam_output_id)
    convolution_output_format = "conv_{}_{{}}.csv".format(sam_output_id)

    lake_file = os.path.join(lakefile_dir, lakefile_format)
    upstream_file = os.path.join(upstream_repository, upstream_format)
    sam_output_file = os.path.join(sam_output_dir, sam_output_format)

    time_of_travel(lake_file, upstream_file, sam_output_file, convolution_output_dir, convolution_output_format)

if __name__ == "__main__":
    main()
    #import cProfile
    #cProfile.run("main()")