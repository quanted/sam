import os
from collections import defaultdict
import numpy as np
import pandas as pd

from Preprocessing.utilities import get_nhd, read_dbf


def read_tables(nhd_folder):
    # Set paths
    erom_dir = os.path.join(nhd_folder, "EROMExtension")
    attribute_file = os.path.join(nhd_folder, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")

    # Initialize output table variable
    out_table = None

    # Read flows
    for month in range(1, 13):
        print(month)

        # Read monthly flow file
        erom_file = os.path.join(erom_dir, "EROM_{}0001.dbf".format(str(month).zfill(2)))
        monthly = \
            pd.DataFrame(data=read_dbf(erom_file, ["COMID", "Q0001E", "V0001E"]),
                         columns=["COMID", "Q{}".format(month), "V{}".format(month)])
        # Adjust units
        monthly["Q{}".format(month)] *= 2446.58  # cfs -> cmd
        monthly["V{}".format(month)] *= 26334.7  # f/s -> md

        # Merge monthly file into master
        if out_table is None:
            out_table = monthly
        else:
            out_table = out_table.merge(monthly)

    # Read lengths
    lengths = pd.DataFrame(data=read_dbf(attribute_file, ["ComID", "LengthKM"]), columns=["COMID", "L"])
    lengths.L *= 1000.  # km -> m
    out_table = out_table.merge(lengths).sort("COMID")

    return out_table


def get_surface_area(master_table):
    # Stream channel geometry
    stream_channel_a = 4.28
    stream_channel_b = 0.55

    for month in range(1, 13):
        cross_section = master_table["Q{}".format(month)] / master_table["V{}".format(month)]
        width = stream_channel_a * np.power(cross_section, stream_channel_b)
        master_table["SA{}".format(month)] = width * master_table.L
    return master_table


def write_table(master_table, outfile, keyfile):
    print("Writing...")
    import time

    # Save COMID index to keyfile
    np.save(keyfile, master_table.COMID.as_matrix())

    # Save necessary fields to matrix
    array = np.memmap(outfile, dtype='float32', mode='w+', shape=(master_table.shape[0], 2, 12))
    index = ["Q{}".format(month) for month in range(1, 13)] + ["SA{}".format(month) for month in range(1, 13)]
    for i, (_, row) in enumerate(master_table[index].iterrows()):
        array[i] = row.reshape((2, 12))

def generate_flow_files(nhd_path, out_folder):
    # Iterate through regions
    for region, region_dir in get_nhd(nhd_path).items():
        print(region)

        # Set output paths
        outfile = os.path.join(out_folder, "region_{}.dat".format(region))
        keyfile = os.path.join(out_folder, "region_{}_key.npy".format(region))

        # Create output directory if it doesn't exist
        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))

        # Read flow data
        master_table = read_tables(region_dir)

        # Compute surface areas from flow and length
        master_table = get_surface_area(master_table)

        # Save master table
        write_table(master_table, outfile, keyfile)

def main():
    nhd_path = r'T:\NationalData\NHDPlusV2'
    out_folder = r'..\bin\Preprocessed\FlowFiles'

    generate_flow_files(nhd_path, out_folder)


if __name__ == "__main__":
    main()
