import numpy as np
import os
import pickle
from collections import defaultdict

from parameters import paths
from read_tools import get_nhd, read_dbf

def parse(path):
    print("\"COMID\" = " + " OR \"COMID\" = ".join(map(str, path)))

def get_upstream_lakes(upstream_file, flowlines, outfile, overwrite=True):

    def open_upstream(path):
        with open(path, 'rb') as f:
            data = np.load(f)
            return data['paths'], data['path_map'], data['conversion_array']

    if not os.path.exists(outfile) or overwrite:
        upstream_lake_dict = defaultdict(set)

        # Unpack upstream file
        paths, path_map, convert = open_upstream(upstream_file)
        convert_back = {val: i for i, val in enumerate(convert)}

        # Create dictionary of reaches containing waterbodies
        waterbodies = np.zeros(len(convert))
        for reach, lake in read_dbf(flowlines, ["ComID", "WBAREACOMI"]):
            if lake:
                alias = convert_back.get(reach)
                if alias:
                    waterbodies[alias] = lake

        # Iterate through paths and assign upstream lakes to reaches
        for alias, address in enumerate(path_map):
            if alias and not alias % 10000:
                print(alias)
            lake_id = waterbodies[alias]
            if any(address) and lake_id:  # If the reach is in a lake and the address is valid...
                try:
                    start_row, end_row, col = address
                    all_upstream = np.int64(paths[start_row:end_row, col:])
                    upstream_lake_dict[lake_id] |= set(np.unique(all_upstream[all_upstream > 0]))
                except Exception as e:
                    print(e)

        with open(outfile, 'wb') as f:
            pickle.dump(dict(upstream_lake_dict), f)



def main():
    region_filter = None
    nhd = get_nhd()
    for region in sorted(nhd.keys() & set(region_filter)) if region_filter else nhd.keys():
        print(region)
        flowlines = os.path.join(nhd[region], "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
        upstream_file = paths.upstream_path.format(region)
        outfile = paths.lentics_path.format(region)
        if os.path.isfile(upstream_file):
            get_upstream_lakes(upstream_file, flowlines, outfile)
        else:
            print("No upstream file found for region {}".format(region))

main()