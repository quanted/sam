import numpy as np
import time
from collections import defaultdict
from read_tools import *
from Tool.travel_time_functions import report_progress
import pickle
import os

# At the moment, reaches with stream_calc = 0 are not processed or included as upstream. Probably should be...
#   HERE'S MY IDEA ON THAT: Include all reaches, regardless of stream_calc.  BUT, only use stream_calc == 1 for
#   intermediate reaches/times.  So, if paths 1234568 and 1234578 exist, but stream_calc[7] == 0,
#   then don't use 7 for 8's time

def snapshot(nodes, outlets, max_length=2500, max_paths=500000):
    diagnoses = []
    diagnostic_i = 0
    all_paths = np.zeros((max_paths, max_length), dtype=np.int32)
    path_cursor = 0
    for i, start_node in enumerate(outlets):
        report_progress(i, len(outlets), 100)
        queue = np.zeros((nodes.shape[0], 2))
        active_node = start_node
        active_reach = np.zeros(max_length, dtype=np.int32)
        queue_cursor = 0
        active_reach_cursor = 0
        while True:
            upstream = nodes[nodes[:,0]==active_node]
            active_reach[active_reach_cursor] = active_node
            active_reach_cursor += 1
            if upstream.size:
                active_node = upstream[0][1]
                for i in range(1, upstream.shape[0]):
                    queue[queue_cursor] = upstream[i]
                    queue_cursor += 1
            else:
                #all_paths[path_cursor] = active_reach
                if 5025559 in active_reach:
                    diagnosis = list(filter(lambda x: x != 0, active_reach))
                    if not diagnosis in diagnoses:
                        diagnoses.append(diagnosis)
                        diagnostic_i += 1
                        print(diagnostic_i)
                # Accumulate up reach
                queue_cursor -= 1
                path_cursor += 1
                last_node, active_node = queue[queue_cursor]
                if not any((last_node, active_node)):
                    break
                active_reach_cursor = np.where(active_reach == last_node)[0][0] + 1

                active_reach[active_reach_cursor:] = 0.


    return all_paths[~np.all(all_paths == 0, axis=1)]

def get_outlets(flow_table, flow_lines, vaa_table):

    # Identify reaches that empty into another NHD region
    in_region = set(read_dbf(flow_lines, ["COMID"]))
    nodes = read_dbf(flow_table, ["FROMCOMID", "TOCOMID"])
    _, ds_reaches = zip(*nodes)
    basin_outlets = set(ds_reaches) - in_region
    basin_outlets = {x[0] for x in nodes if x[1] in basin_outlets}

    # Identify all the reaches that are designated as a terminal path (has to be converted from HydroSeq)
    c, h, t = zip(*read_dbf(vaa_table, ["ComID", "HydroSeq", "TerminalPa"]))
    terminal_paths = set(map(dict(zip(h, c)).get, set(t) & set(h)))

    return (terminal_paths | basin_outlets) - {0}

def get_nodes(flow_table):
    nodes = read_dbf(flow_table, ["TOCOMID", "FROMCOMID"])  # Get all nodes
    nodes = filter(all, nodes)  # Filter out nodes where from or to value is zero
    return np.array(sorted(list(nodes)), dtype=np.int32)  # Return a sorted array of nodes

def get_times(paths, q_table, vaa_table):

    # Read in velocity and length data into lookup tables
    velocity_dict = dict(read_dbf(q_table, ["COMID", "V0001E"]))
    length_dict = dict(read_dbf(vaa_table, ["ComID", "LengthKM"]))

    # Create lookup dictionary with times for each reach
    time_dict = {}
    for reach in np.unique(paths):  # loop through all reaches
        length = length_dict.get(reach, 0) * 1000.  # km -> m
        velocity = velocity_dict.get(reach, 0) * 0.3048  # ft/sec -> m/sec
        if length and velocity:
            time_dict[reach] = (length / velocity)  / 86400.  # seconds -> days
    time_dict[0] = 0.

    # Assign times to each reach in paths
    times = np.zeros_like(paths, dtype=np.float)
    for i, path in enumerate(paths):
        path_end = np.argmax(path == 0)
        times[i][:path_end] = np.cumsum([time_dict.get(r, 0.0) for r in path[:path_end]])
    return times


def get_lakes(paths, flow_lines):
    wb_dict = dict(read_dbf(flow_lines, ["COMID", "WBAREACOMI"]))
    map_lakes = np.vectorize(lambda x: wb_dict.get(x, 0))
    lake_map = map_lakes(paths)
    return lake_map


def map_paths(paths):
    path_dict = defaultdict(lambda: defaultdict(list))
    for i, path in enumerate(paths):
        for j, val in enumerate(path):
            if val:
                path_dict[val][j].append(i)
    return path_dict

def write_to_outfiles(output_dir, region, paths, times, lakes, path_map):
    upstream_cube = np.dstack((paths, times, lakes))
    outfile = os.path.join(output_dir, "upstream_{}.p".format(region))
    with open(outfile, 'wb') as f:
        pickle.dump((upstream_cube, path_map), f)


def main(nhd_dir, output_directory, region_filter='all'):

    for region, region_dir in get_nhd(nhd_dir).items():
        if region == region_filter or region_filter == 'all':
            print("Processing region {}...".format(region))
            start = time.time()

            # Designate paths
            flow_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlow.dbf")
            flow_lines = os.path.join(region_dir, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
            vaa_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
            q_table = os.path.join(region_dir, "EROMExtension", "EROM_MA0001.dbf")
            assert all(map(os.path.isfile, (flow_table, flow_lines, vaa_table, q_table)))

            # Do the work
            nodes = get_nodes(flow_table)
            outlets = get_outlets(flow_table, flow_lines, vaa_table)
            paths = snapshot(nodes, outlets)
            times = get_times(paths, q_table, flow_lines)
            lakes = get_lakes(paths, flow_lines)
            path_map = map_paths(paths, times)

            # Write to files
            write_to_outfiles(output_directory, region, paths, times, lakes, path_map)

            total_time = int(time.time() - start)
            print("\tComplete in {}:{}".format(int(total_time / 60), str(total_time % 60).zfill(2)))

if __name__ == "__main__":
    output_directory = r"T:\SAM\Preprocessed\UpstreamPaths"
    nhd_dir = r"T:\NationalData\NHDPlusV2"
    region_filter = '07'
    main(nhd_dir, output_directory, region_filter)