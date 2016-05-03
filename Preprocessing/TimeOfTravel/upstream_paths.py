import numpy as np
import os
import sys
import pickle
from collections import defaultdict

from read_tools import get_nhd, read_dbf

def snapshot(nodes, outlets, time_dict, max_length=2500, max_paths=500000):
    all_paths = np.zeros((max_paths, max_length), dtype=np.int32)
    all_times = np.zeros((max_paths, max_length), dtype=np.float)
    path_cursor = 0
    for index, start_node in enumerate(outlets):
        queue = np.zeros((nodes.shape[0], 2))
        active_node = start_node
        active_reach = np.zeros(max_length, dtype=np.int32)
        active_times = np.zeros(max_length, dtype=np.float)
        start_cursor = 0
        queue_cursor = 0
        active_reach_cursor = 0
        while True:
            upstream = nodes[nodes[:,0]==active_node]
            active_reach[active_reach_cursor] = active_node
            active_times[active_reach_cursor] = time_dict.get(active_node, 0.0)
            active_reach_cursor += 1
            if upstream.size:
                active_node = upstream[0][1]
                for j in range(1, upstream.shape[0]):
                    queue[queue_cursor] = upstream[j]
                    queue_cursor += 1
            else:
                # Accumulate up reach
                all_paths[path_cursor, start_cursor:] = active_reach[start_cursor:]
                all_times[path_cursor] = np.cumsum(active_times) * (all_paths[path_cursor] > 0)
                queue_cursor -= 1
                path_cursor += 1
                last_node, active_node = queue[queue_cursor]
                if not any((last_node, active_node)):
                    break
                active_reach_cursor = np.where(active_reach == last_node)[0][0] + 1
                start_cursor = active_reach_cursor
                active_reach[active_reach_cursor:] = 0.
                active_times[active_reach_cursor:] = 0.

    return all_paths[:path_cursor], all_times[:path_cursor]


def get_outlets(flow_table, flow_lines, vaa_table, conversion_dict):

    # Identify reaches that empty into another NHD region
    in_region = set(read_dbf(flow_lines, ["COMID"]))
    nodes = read_dbf(flow_table, ["FROMCOMID", "TOCOMID"])
    _, ds_reaches = zip(*nodes)
    basin_outlets = set(ds_reaches) - in_region
    basin_outlets = {x[0] for x in nodes if x[1] in basin_outlets}

    # Identify all the reaches that are designated as a terminal path (has to be converted from HydroSeq)
    c, h, t = zip(*read_dbf(vaa_table, ["ComID", "HydroSeq", "TerminalPa"]))
    terminal_paths = set(map(dict(zip(h, c)).get, set(t) & set(h)))
    outlets = (terminal_paths | basin_outlets) - {0}

    return set(filter(None, map(conversion_dict.get, outlets)))


def get_nodes(flow_table, vaa_table):
    # Get a set of nodes
    flows = read_dbf(flow_table, ["TOCOMID", "FROMCOMID"])  # Get all nodes
    sc = dict(read_dbf(vaa_table, ["ComID", "StreamCalc"]))
    nodes = list(filter(all, flows))  # Filter out nodes where from or to value is zero

    # Convert to indices
    unique_comids = np.unique(nodes)
    conversion_dict = {comid: i + 1 for i, comid in enumerate(unique_comids)}
    nodes = np.vectorize(lambda x: conversion_dict.get(x))(nodes)
    sc = dict(zip(map(conversion_dict.get, sc.keys()), sc.values()))

    sc_nodes = np.array([bool(sc.get(n)) for n in nodes.flat]).reshape(nodes.shape)
    node_sums = np.sum(sc_nodes, axis=1)
    active_nodes = nodes[node_sums == 2]
    connecting_nodes = nodes[(sc_nodes[:,0] == True) & (sc_nodes[:, 1] == False)]
    passive_nodes = nodes[node_sums == 0]


    return active_nodes, connecting_nodes, passive_nodes, conversion_dict


def map_paths(paths):

    # Get starting row and column for each value
    column_numbers = np.tile(np.arange(paths.shape[1]) + 1, (paths.shape[0], 1)) * (paths > 0)
    path_begins = np.argmax(column_numbers > 0, axis=1)

    path_map = defaultdict(list)
    n_paths = paths.shape[0]
    for i, path in enumerate(paths):
        for j, val in enumerate(path):
            if val:
                if i == n_paths:
                    end_row = 0
                else:
                    next_row = (path_begins[i + 1:] <= j)
                    if next_row.any():
                        end_row = np.argmax(next_row)
                    else:
                        end_row = n_paths - i - 1
                path_map[val] = (i, i + end_row + 1, j)

    return path_map


def get_times(all_nodes, q_table, vaa_table, conversion_dict):

    # Read in velocity and length data into lookup tables
    velocities = tuple(zip(*read_dbf(q_table, ["COMID", "V0001E"])))
    lengths = tuple(zip(*read_dbf(vaa_table, ["ComID", "LengthKM"])))

    velocity_dict = dict(zip(map(conversion_dict.get, velocities[0]), velocities[1]))
    length_dict = dict(zip(map(conversion_dict.get, lengths[0]), lengths[1]))

    velocity_dict.pop(None, None)
    length_dict.pop(None, None)

    # Create lookup dictionary with times for each reach
    time_dict = {}
    all_reaches = {reach for node in all_nodes for reach in node if reach}
    for reach in all_reaches:  # loop through all reaches
        length = length_dict.get(reach, 0) * 1000.  # km -> m
        velocity = velocity_dict.get(reach, 0) * 0.3048  # ft/sec -> m/sec
        if length and velocity:
            time_dict[reach] = (length / velocity) / 86400.  # seconds -> days

    return time_dict


def write_to_outfiles(output_dir, region, paths, times, path_map, conversion_dict):
    outfile = os.path.join(output_dir, "upstream_{}.p".format(region))
    with open(outfile, 'wb') as f:
        pickle.dump((paths, times, path_map, conversion_dict), f)


def trim_paths(paths):
    longest_path = np.max(paths.nonzero()[1])
    return paths[:, :longest_path]


def create_upstream_paths(nhd_dir, output_directory, region_filter='all'):

    for region, region_dir in get_nhd(nhd_dir).items():
        if region == region_filter or region_filter == 'all':
            print("Processing region {}...".format(region))

            # Designate paths
            flow_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlow.dbf")
            flow_lines = os.path.join(region_dir, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf")
            vaa_table = os.path.join(region_dir, "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
            q_table = os.path.join(region_dir, "EROMExtension", "EROM_MA0001.dbf")

            # Do the work
            active_nodes, connectors, passive_nodes, conversion_dict = get_nodes(flow_table, vaa_table)
            time_dict = get_times(active_nodes, q_table, vaa_table, conversion_dict)
            outlets = get_outlets(flow_table, flow_lines, vaa_table, conversion_dict)
            paths, times = snapshot(active_nodes, outlets, time_dict)
            path_map = map_paths(paths)
            write_to_outfiles(output_directory, region, paths, times, path_map, conversion_dict)


def main():
    output_directory = r"T:\pySAM\bin\Preprocessed\UpstreamPaths"
    nhd_dir = r"T:\NationalData\NHDPlusV2"
    region_filter = '07'
    create_upstream_paths(nhd_dir, output_directory, region_filter)

if __name__ == "__main__":
    time_it = False
    if time_it:
        import cProfile
        cProfile.run('main()')
    else:
        main()
