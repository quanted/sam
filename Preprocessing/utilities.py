import os
import numpy as np
from collections import OrderedDict, defaultdict


class Navigator(object):
    def __init__(self, region, paths_dir=r"T:\TripTools\WatershedTools\Paths\Collapsed", form="upstream_{}.npz"):
        self.region = region
        data = np.load(os.path.join(paths_dir, form.format(self.region)))
        self.paths, self.path_map, self.start_cols, self.alias_to_comid = \
            data['paths'], data['path_map'], data['start_cols'], data['conversion_array']
        self.n_reaches = self.alias_to_comid.size
        self.comid_to_alias = dict(zip(self.alias_to_comid, np.arange(self.n_reaches)))
        self.reaches = np.sort(self.alias_to_comid)
        self.aliases = np.arange(self.reaches.size)
        self.last_outlet = None

    def all_upstream(self, reach, mode='comid'):
        reach = self._format_input(reach, mode)
        if reach:
            start_row, end_row, col = map(int, self.path_map[reach])
            start_col = list(self.paths[start_row]).index(reach)
            upstream_reaches = list(self.paths[start_row:end_row])
            upstream_reaches.append(self.paths[start_row][start_col:])
            output = np.concatenate(upstream_reaches)
            return self._format_output(output, mode)
        else:
            return np.array([])

    def all_downstream(self, reach, mode='comid'):
        if not self.last_outlet:
            a = (self.start_cols == 0)
            self.last_outlet = np.where(a)[0][a.cumsum() - 1]
        reach = self._format_input(reach, mode)
        start_row, _, _ = map(int, self.path_map[reach])
        last_outlet = self.last_outlet[start_row]
        start_cols = self.start_cols[last_outlet:start_row+1]
        active_paths = np.where(start_cols == np.minimum.accumulate(start_cols[::-1])[::-1])[0] + last_outlet
        output = np.zeros((start_cols[-1] + len(self.paths[start_row])) * 1.5, dtype=np.int32)
        for i, start, end in zip(active_paths, self.start_cols[active_paths], self.start_cols[active_paths][1:]):
            output[start:end] = self.paths[i][:end - start]
        return self._format_output(output[output > 0], mode)

    def upstream_paths(self, reach, mode='comid'):
        reach = self._format_input(reach, mode)
        start_row, end_row, _ = map(int, self.path_map[reach])
        path = np.zeros(max(map(lambda x: self.start_cols[x] + self.paths[x].size, range(start_row, end_row))), dtype=np.int32)
        baseline = self.start_cols[start_row]
        for i in range(start_row, end_row - 1):
            stub = self.paths[i]
            start_col = self.start_cols[i] - baseline
            path[start_col:] = 0
            path[start_col:start_col + stub.size] = stub
            output = path[path != 0]
            yield self._format_output(output, mode)

    def _format_input(self, reach_id, mode):
        return reach_id if mode == 'alias' else self.comid_to_alias.get(reach_id)

    def _format_output(self, output, mode):
        return output if mode == 'alias' else self.alias_to_comid[output]

    def __repr__(self):
        return "NHD Region {} Navigator".format(self.region)

def gis_it(comids, field="COMID", lookup_dict=None):

    if lookup_dict:
        comids = filter(None, map(lookup_dict.get, comids))

    print("\"{}\" = ".format(field) + " OR \"{}\" = ".format(field).join(map(str, comids)))

# Reads the contents of a dbf table
def read_dbf(dbf_file, out_fields=None):
    import ogr

    # Initialize file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(dbf_file)
    layer = data_source.GetLayer(0)

    # Set fields
    ld = layer.GetLayerDefn()
    fields = {ld.GetFieldDefn(i).GetName() for i in range(ld.GetFieldCount())}
    if out_fields:
        missing_fields = set(out_fields) - fields
        remedies = [lambda x: x.upper(), lambda x: x.capitalize()]
        while missing_fields and remedies:
            new_fields = map(remedies.pop(0), out_fields)
            out_fields = [nf if nf in fields else out_fields[i] for i, nf in enumerate(new_fields)]
            missing_fields = set(out_fields) - fields
        if missing_fields:
            print("Fields {} not found in {}".format(out_fields, dbf_file))
        fields = [field for field in out_fields if not field in missing_fields]

    # Read data
    if len(fields) > 1:
        table = [[row.GetField(f) for f in fields] for row in layer]
    else:
        table = [row.GetField(list(fields)[0]) for row in layer]
    # noinspection PyUnusedLocal
    data_source = None
    return table


def gis_it(comids, field="COMID", lookup_dict=None):

    if lookup_dict:
        comids = filter(None, map(lookup_dict.get, comids))

    print("\"{}\" = ".format(field) + " OR \"{}\" = ".format(field).join(map(str, comids)))

# Assembles a dictionary of NHD Plus directory structure indexed by region
def get_nhd(nhd_dir=r"T:\NationalData\NHDPlusV2", region_filter='all'):
    all_paths = defaultdict()
    regions = list(nhd_states.keys())
    region_dirs = {"NHDPlus{}".format(region) for region in regions}
    for root_dir, sub_dirs, _ in os.walk(nhd_dir):
        if set(sub_dirs) & region_dirs:
            for sub_dir in sub_dirs:
                region = sub_dir.lstrip("NHDPlus")
                if region in regions:
                    all_paths[sub_dir.lstrip("NHDPlus")] = os.path.join(root_dir, sub_dir)
    return OrderedDict(sorted(all_paths.items()))


nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                          ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                          ('03N', {"VA", "NC", "SC", "GA"}),
                          ('03S', {"FL", "GA"}),
                          ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                          ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                          ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                          ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                          ('07', {"MN", "WI", "SD", "IA", "IL", "MO"}),
                          ('08', {"MO", "KY", "TN", "AR", "MS", "LA"}),
                          ('09', {"ND", "MN", "SD"}),
                          ('10U', {"MT", "ND", "WY", "SD", "MN", "NE", "IA"}),
                          ('10L', {"CO", "WY", "MN", "NE", "IA", "KS", "MO"}),
                          ('11', {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"}),
                          ('12', {"NM", "TX", "LA"}),
                          ('13', {"CO", "NM", "TX"}),
                          ('14', {"WY", "UT", "CO", "AZ", "NM"}),
                          ('15', {"NV", "UT", "AZ", "NM", "CA"}),
                          ('16', {"CA", "OR", "ID", "WY", "NV", "UT"}),
                          ('17', {"WA", "ID", "MT", "OR", "WY", "UT", "NV"}),
                          ('18', {"OR", "NV", "CA"})))

