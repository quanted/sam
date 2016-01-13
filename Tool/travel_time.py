import os

import travel_time_functions as functions
import read
import write

# @@@ - to do:
# class for indexed array?  class for input file with dir and format?  OR structured arrays?
# how to set up?  will SAM be run for whole regions?
# still not clear on flattening hydrograph for lake outlets

def time_of_travel(lake_file, upstream_file, sam_output_file, output_dir, output_format, diagnostics=False):

    # Read in compressed SAM output
    sam_output, sam_lookup, dates = read.unpickle(sam_output_file)

    # Read lake file
    volume_dict, outflow_dict, outlet_dict, residence_times, waterbody_dict = read.unpickle(lake_file)

    # Get upstream path data
    upstream, path_map = read.unpickle(upstream_file)

    # Separate flowing reaches from reaches that are in a reservoir
    lentic_reaches, lotic_reaches, upstream_lentics = \
        functions.classify_reaches(sam_lookup.keys(), outlet_dict, waterbody_dict)

    # Loop through all reaches with SAM output and run convolution
    for i, reach in enumerate(lotic_reaches):

        reach_address = path_map.get(reach)

        if reach_address:

            # Pull upstream path data for this reach
            upstream_paths, upstream_lakes, upstream_times = functions.trim_to_reach(upstream, *reach_address)

            # Get a set of all unique reaches upstream of the reach,
            if upstream_paths.any():  # Check that reach was found in the path matrix.

                # Get baseflow for the reach.  This is not convolved or accumulated upstream
                baseflow = sam_output[2, sam_lookup[reach]]

                # Get all applicable SAM output time series for this reach
                local_output, local_lookup = \
                    functions.trim_to_upstream(sam_output, sam_lookup, upstream_paths)

                # Create an un-convolved output dataset for comparison
                if diagnostics:
                    functions.preconvolution_report(reach, dates, output_dir, output_format, local_output, baseflow)

                # Convolve paths that pass through reservoirs
                local_output = \
                    functions.convolve_reservoirs(local_output, local_lookup, upstream_paths, upstream_lakes,
                                                  residence_times)

                # Snap the travel times to an interval (portions of a day).  Default is 1-day
                runoff_mass, total_runoff = \
                    functions.convolve_flowing(upstream_paths, upstream_times, local_output, local_lookup)

                # Unpack time series and compute concentration
                runoff_mass, total_runoff, baseflow, total_flow, total_conc, runoff_conc = \
                    functions.compute_concentration(runoff_mass, total_runoff, baseflow)

                # If the reach is an outlet for a lake, write the output to the reach and all reaches in the lake
                for reach in {reach} | upstream_lentics[reach]:
                    write.daily(output_dir, output_format, reach, total_conc, runoff_conc, runoff_mass, dates,
                                 total_flow, baseflow, total_runoff)

            else:
                print("{} not found in path map".format(reach))


def main():

    diagnostics = True

    region = '07'
    #sam_output_id = "mark_twain"
    sam_output_id = "mtb_test01"

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

    time_of_travel(lake_file, upstream_file, sam_output_file, convolution_output_dir, convolution_output_format,
                   diagnostics)

if __name__ == "__main__":
    main()