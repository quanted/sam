import numpy as np
import os
from collections import namedtuple


class Hydroregion:
    """
    Contains all datasets and functions related to the NHD Plus region, including all hydrological features and links
    between them, as well as the configuration of all reach catchments (recipes) and the scenarios that they contain.
    Contains many file input functions.
    """

    def __init__(self, upstream_dir, region):

        from parameters import time_of_travel
        self.id = region

        self.minimum_residence_time = time_of_travel.minimum_residence_time
        self.round_func = np.int32 if time_of_travel.round_down else np.vectorize(lambda x: np.int32(np.round(x)))

        # Read hydrological input files
        self.upstream = self.read_upstream_file(upstream_dir)

    def upstream_watershed(self, reach_id, mode='reach', return_times=True):

        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        # Look up reach ID and fetch address from upstream object
        reach = reach_id if mode == 'alias' else self.upstream.reach_to_alias.get(reach_id)
        start_row, end_row, col = map(int, self.upstream.map[reach])
        start_col = list(self.upstream.paths[start_row]).index(reach)

        # Fetch upstream reaches and times
        aliases = unpack(self.upstream.paths)
        reaches = aliases if mode == 'alias' else np.int32(self.upstream.alias_to_reach[aliases])
        if not return_times:
            return reaches
        else:
            times = unpack(self.upstream.times)
            adjusted_times = self.round_func(times - self.upstream.times[start_row][start_col])
            return reaches, adjusted_times

    def cascade(self):
        run_reaches = set()
        confined_lake_table = self.lake_table[self.lake_table.index.isin(self.active_lakes)]
        for _, lake in confined_lake_table.iterrows():
            upstream_reaches = set(self.upstream_watershed(lake.OutletID, mode='reach', return_times=False))
            reaches = upstream_reaches & self.active_reaches
            yield reaches, lake
            run_reaches |= reaches
        remaining_reaches = self.active_reaches - run_reaches
        yield remaining_reaches, None

    def read_upstream_file(self, upstream_path):
        upstream_file = upstream_path.format(self.id)
        assert os.path.isfile(upstream_file), "Upstream file {} not found".format(upstream_file)
        Upstream = namedtuple("Upstream", ["paths", "times", "map", "alias_to_reach", "reach_to_alias"])
        data = np.load(upstream_file, mmap_mode='r')
        conversion_array = data['alias_index']
        reverse_conversion = dict(zip(conversion_array, np.arange(conversion_array.size)))
        return Upstream(data['paths'], data['time'], data['path_map'], conversion_array, reverse_conversion)


    def upstream_watershed(self, reach_id, mode='reach', return_times=True):
        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        # Look up reach ID and fetch address from upstream object
        reach = reach_id if mode == 'alias' else self.upstream.reach_to_alias.get(reach_id)
        start_row, end_row, col = map(int, self.upstream.map[reach])
        start_col = list(self.upstream.paths[start_row]).index(reach)

        # Fetch upstream reaches and times
        aliases = unpack(self.upstream.paths)
        reaches = aliases if mode == 'alias' else np.int32(self.upstream.alias_to_reach[aliases])
        if not return_times:
            return reaches
        else:
            times = unpack(self.upstream.times)
            adjusted_times = self.round_func(times - self.upstream.times[start_row][start_col])
            return reaches, adjusted_times


hydro_dir = r"S:\bin\Preprocessed\Upstream\upstream_{}.npz"
region = '07'
test_reach = 4867727

region = Hydroregion(hydro_dir, region)

paths, times = region.upstream_watershed(test_reach)

print(np.unique(times))