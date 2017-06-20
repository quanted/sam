import os
from collections import defaultdict
import numpy as np

from Preprocessing.utilities import get_nhd, read_dbf


def get_flows(nhd_folder):
    flows = defaultdict(lambda: np.zeros(24))
    erom_dir = os.path.join(nhd_folder, "EROMExtension")
    for month in range(1, 13):
        erom_file = os.path.join(erom_dir, "EROM_{}0001.dbf".format(str(month).zfill(2)))
        if os.path.isfile(erom_file):
            for comid, q, v in read_dbf(erom_file, ["COMID", "Q0001E", "V0001E"]):
                flows[comid][month - 1] = q * 2446.58  # cfs -> cmd
                flows[comid][month + 11] = v * 26334.7  # f/s -> md
    return flows


def get_lengths(nhd_folder):
    attribute_file = os.path.join(nhd_folder, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
    return dict(read_dbf(attribute_file, ["ComID", "LengthKM"]))


def compile_and_compute(flows, lengths):
    all_found = set(flows.keys()) & set(lengths.keys())
    not_found = (flows.keys() | lengths.keys()) - all_found
    if not_found:
        print("{} of {} reaches could not be matched".format(len(not_found), len(all_found) + len(not_found)))
    out_data = []
    for comid in all_found:
        l = lengths[comid]
        row = np.concatenate([np.array([comid, l]), flows[comid]])
        out_data.append(row)
    out_data = np.array(out_data)
    return out_data


def generate_flow_files(nhd_path, out_folder):
    for region, region_dir in get_nhd(nhd_path).items():
        print(region)
        outfile = os.path.join(out_folder, "region_{}.dat".format(region))
        keyfile = os.path.join(out_folder, "region_{}_key.npy".format(region))
        if not os.path.exists(os.path.dirname(outfile)):
            os.makedirs(os.path.dirname(outfile))
        flows = get_flows(region_dir)
        lengths = get_lengths(region_dir)
        all_data = compile_and_compute(flows, lengths)
        array = np.memmap(outfile, dtype='float32', mode='w+', shape=all_data.shape)
        array[:] = all_data
        shape = np.array(all_data.shape)
        key = np.concatenate((shape, np.int32(all_data[:, 0])))
        np.save(keyfile, key)


def main():
    nhd_path = r'T:\NationalData\NHDPlusV2'
    out_folder = r'..\bin\Preprocessed\FlowFiles'
    generate_flow_files(nhd_path, out_folder)


if __name__ == "__main__":
    main()
