import pandas as pd
import numpy as np
import os

from Tool.functions import Navigator
from Preprocessing.utilities import read_dbf


class NHDTable(object):
    def __init__(self, region, path):
        from fields import gridcode_fields, flowline_fields, vaa_fields
        self.region = region
        self.path = path.format(region)
        self.tables = \
            [(os.path.join(self.path, "NHDPlusAttributes", "PlusFlow.dbf"), None),
             (os.path.join(self.path, "NHDSnapshot", "Hydrography", "NHDFlowline.dbf"), flowline_fields),
             (os.path.join(self.path, "NHDPlusAttributes", "PlusFlowlineVAA.dbf"), vaa_fields),
             (os.path.join(self.path, "NHDPlusCatchment", "featureidgridcode.dbf"), gridcode_fields)]
        self.erom_dir = os.path.join(self.path, "EROMExtension", "EROM_{}0001.dbf")

        self.matrix = self.extract_tables()
        self.add_erom()
        self.calculate_attributes()

    def extract_tables(self):
        all_tables = None
        for table, fields in self.tables:
            table = pd.read_csv(table)
            if fields is not None:
                table = table.rename(columns=fields.convert)[fields.new]
            if all_tables is None:
                all_tables = table
            else:
                all_tables = all_tables.merge(table, on='comid')

        return all_tables[all_tables.comid != 0]

    def add_erom(self):
        master_table = None
        months = list(map(str, range(1, 13))) + ['ma']
        for month in months:
            table_path = self.erom_dir.format(month.zfill(2))
            fields = ["comid", "q0001e", "v0001e"] if month == 'ma' else ["comid", "q0001e"]
            erom_table = read_dbf(table_path)[fields]
            erom_table.columns = [c.replace("0001e", month) for c in erom_table.columns]
            if master_table is None:
                master_table = erom_table
            else:
                master_table = pd.merge(master_table, erom_table, on='comid')
        return self.matrix.merge(master_table, on='comid')

    def calculate_attributes(self):
        months = list(map(str, range(1, 13)))

        # Convert units and calculate travel times
        self.matrix['length'] = self.matrix.pop('lengthkm') * 1000.  # km -> m
        for month in months:
            self.matrix["q{}".format(month)] *= 2446.58  # cfs -> cmd
        self.matrix["vma"] *= 26334.7  # f/s -> md
        self.matrix["travel_time"] = self.matrix.length / self.matrix.vma

        # Calculate surface area
        stream_channel_a = 4.28
        stream_channel_b = 0.55
        cross_section = self.matrix.qma / self.matrix.pop('vma')
        self.matrix['surface_area'] = stream_channel_a * np.power(cross_section, stream_channel_b)

        # Indicate whether reaches are coastal
        self.matrix['coastal'] = np.int16(self.matrix.pop('fcode') == 56600)


class NavigatorBuilder(object):
    def __init__(self, nhd_table, output_path):
        nodes, times, conversion_dict = self.get_nodes(nhd_table)
        attributes = np.array([times])

        outlets = self.identify_outlets(nhd_table, conversion_dict)

        paths, attribute_matrix = self.upstream_trace(nodes, outlets, attributes)

        path_map = self.map_paths(paths)

        paths, attribute_matrix, start_cols = self.collapse_array(paths, attribute_matrix)

        conversion_array = self.dict_to_array(conversion_dict)

        self.write_outfile(output_path, paths, path_map, conversion_array, attribute_matrix[0])

    def identify_outlets(self, comid_table):
        """ Identify the outlet reach corresponding to each reservoir """
        # Filter the reach table down to only outlet reaches by getting the highest hydroseq for each wb_comid
        outlets = comid_table.loc[comid_table.groupby(["wb_comid"])["hydroseq"].idxmin()]
        outlets = outlets[[f for f in outlets.columns if f != 'hydroseq']].rename(columns={'comid': 'outlet_comid'})
        return outlets

    def get_nodes(self, nhd_table):
        nodes = np.int32(nhd_table[["tocomid", "comid"]].as_matrix())
        comids, index = np.unique(nodes.T[1], return_index=True)
        times = nhd_table["travel_time"].as_matrix()[index]
        conversion_dict = dict(zip(comids, np.arange(comids.size)))
        nodes = np.vectorize(lambda x: conversion_dict.get(x, -1))(nodes)

        return nodes, times, conversion_dict

    def upstream_trace(self, nodes, outlets, attributes=np.array([]), max_length=3000, max_paths=500000):
        n_attributes = attributes.shape[0]

        # Output arrays
        all_paths = np.zeros((max_paths, max_length), dtype=np.int32)
        all_attributes = np.zeros((n_attributes, max_paths, max_length), dtype=np.float)

        # Bounds
        path_cursor = 0
        longest_path = 0

        for index, start_node in enumerate(outlets):

            # Reset everything except the master. Trace is done separately for each outlet
            # Holders
            queue = np.zeros((nodes.shape[0], 2), dtype=np.int32)
            active_reach = np.zeros(max_length, dtype=np.int32)
            active_attributes = np.zeros((n_attributes, max_length), dtype=np.float)

            # Cursors
            start_cursor = 0
            queue_cursor = 0
            active_reach_cursor = 0
            active_node = start_node

            while True:
                upstream = nodes[nodes[:, 0] == active_node]
                active_reach[active_reach_cursor] = active_node
                for i, attribute_array in enumerate(attributes):
                    active_attributes[i, active_reach_cursor] = attribute_array[active_node]
                active_reach_cursor += 1
                if active_reach_cursor > longest_path:
                    longest_path = active_reach_cursor
                if upstream.size:
                    active_node = upstream[0][1]
                    for j in range(1, upstream.shape[0]):
                        queue[queue_cursor] = upstream[j]
                        queue_cursor += 1
                else:
                    all_paths[path_cursor, start_cursor:] = active_reach[start_cursor:]
                    for i in range(n_attributes):
                        all_attributes[i, path_cursor] = np.cumsum(active_attributes[i]) * (all_paths[path_cursor] > 0)
                    queue_cursor -= 1
                    path_cursor += 1
                    last_node, active_node = queue[queue_cursor]
                    if not any((last_node, active_node)):
                        break
                    active_reach_cursor = np.where(active_reach == last_node)[0][0] + 1
                    start_cursor = active_reach_cursor
                    active_reach[active_reach_cursor:] = 0.
                    for i in range(n_attributes):
                        active_attributes[i, active_reach_cursor:] = 0.

        return all_paths[:path_cursor, :longest_path + 1], all_attributes[:, :path_cursor, :longest_path + 1]

    def map_paths(self, paths):
        # Get starting row and column for each value
        column_numbers = np.tile(np.arange(paths.shape[1]) + 1, (paths.shape[0], 1)) * (paths > 0)
        path_begins = np.argmax(column_numbers > 0, axis=1)
        max_reach = np.max(paths)
        path_map = np.zeros((max_reach + 1, 3))
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
                    values = np.array([i, i + end_row + 1, j])
                    path_map[val] = values

        return path_map

    def collapse_array(self, paths, attribute_matrix):
        out_paths = []
        out_attributes = [[] for _ in range(attribute_matrix.shape[0])]
        path_starts = []
        for i, row in enumerate(paths):
            active_path = (row > 0)
            path_starts.append(np.argmax(active_path))
            out_paths.append(row[active_path])
            for j in range(attribute_matrix.shape[0]):
                out_attributes[j].append(attribute_matrix[j][i][active_path])
        out_attributes = np.array(out_attributes)
        return np.array(out_paths), np.array(out_attributes), np.array(path_starts)

    def write_outfile(self, outfile, paths, path_map, conversion_array, times):
        np.savez_compressed(outfile, paths=paths, path_map=path_map, alias_index=conversion_array, time=times)

    def dict_to_array(self, cd):
        comids, aliases = map(np.array, zip(*cd.items()))
        out_array = np.zeros(max(aliases) + 1)
        out_array[aliases[aliases > 0]] = comids[aliases > 0]
        return out_array


class LakeFileBuilder(object):
    def __init__(self, nhd_table, volume_table_path, nav, outfile_path):
        # Get a table of all lentic reaches, with the COMID of the reach and waterbody
        reservoir_table = nhd_table[["comid", "wb_comid", "hydroseq", "qma"]].rename(columns={'qma': 'flow'})

        # Get the outlets for each reservoir
        reservoir_table = self.identify_outlets(reservoir_table)

        # Get residence times
        reservoir_table = self.get_residence_times(reservoir_table, volume_table_path)

        # Count number of reservoirs upstream of each reservoir
        reservoir_table = self.make_lake_bins(nav, reservoir_table)

        # Save table
        self.save_table(reservoir_table, outfile_path)

    def identify_outlets(self, nhd_table, conversion_dict):
        # Identify reaches that drain to a coastline or a reach not in the region
        non_coastal = nhd_table[nhd_table.coastal == 0]
        outlets = non_coastal[~non_coastal['tocomid'].isin(non_coastal.comid)].comid.as_matrix()
        outlets = np.vectorize(lambda x: conversion_dict.get(x))(outlets)

        return outlets

    def get_residence_times(self, reservoir_table, volume_path):
        # Read and reformat volume table
        volume_table = read_dbf(volume_path)[["comid", "volumecorr"]]
        volume_table = volume_table.rename(columns={"comid": "wb_comid", "volumecorr": "volume"})

        # Join reservoir table with volumes
        joined_table = pd.merge(reservoir_table, volume_table, on="wb_comid")
        joined_table['residence_time'] = joined_table.volume / (joined_table.flow * 0.0283168) / 86400.

        del joined_table['volume']

        return joined_table

    def make_lake_bins(self, nav, reservoir_table):
        lentic_reaches = set(reservoir_table.outlet_comid)
        reservoir_table['n_upstream'] = 0
        for index, row in reservoir_table.iterrows():
            upstream_reaches, warning = nav.upstream_watershed(row.outlet_comid, return_times=False)
            upstream_lentics = set(map(float, upstream_reaches)) & lentic_reaches
            if upstream_lentics:
                reservoir_table.loc[index, 'n_upstream'] = \
                    reservoir_table['outlet_comid'].isin(upstream_lentics).sum()
        return reservoir_table

    def save_table(self, reservoir_table, outfile_path):
        reservoir_table[["outlet_comid", "wb_comid"]] = \
            np.int32(reservoir_table[["outlet_comid", "wb_comid"]].as_matrix())
        np.savez_compressed(outfile_path, table=reservoir_table.as_matrix(), key=reservoir_table.columns)


class FlowFileBuilder(object):
    def __init__(self, nhd_table, out_table):
        from fields import flow_file_fields
        self.fields = flow_file_fields
        self.table = nhd_table[self.fields]
        self.out_table = out_table

        self.save()

    def save(self):
        np.savez_compressed(self.out_table, table=self.table.as_matrix(), key=self.table.columns.tolist())


def main():
    from Preprocessing.utilities import nhd_states

    # Input paths
    nhd_path = os.path.join("..", "bin", "Preprocessed", "CondensedNHD")
    volume_path = os.path.join("..", "bin", "Tables", "LakeMorphometry", "region_{}.dbf")

    # Output paths
    lake_output_path = os.path.join("..", "bin", "Preprocessed", "LakeFiles", "region_{}.npz")
    flow_output_path = os.path.join("..", "bin", "Preprocessed", "FlowFiles", "region_{}.npz")
    navigator_output_path = os.path.join("..", "bin", "Preprocessed", "Navigators", "region_{}.npz")

    # Set run parameters
    build_navigator = True
    build_lake_file = True
    build_flow_file = True

    # Loop through regions
    for region in nhd_states.keys():

        nhd_table = NHDTable(region, nhd_path).matrix

        if build_flow_file:
            FlowFileBuilder(nhd_table, flow_output_path.format(region))

        if build_navigator:
            NavigatorBuilder(nhd_table, navigator_output_path.format(region))

        if build_lake_file:
            nav = Navigator(region, navigator_output_path)
            LakeFileBuilder(nhd_table, volume_path.format(region), nav, lake_output_path.format(region))


main()
