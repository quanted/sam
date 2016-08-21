import math
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import copy

from Tool import read, write
from Tool.parameters import time_of_travel as tot


class ImpulseResponse:
    def __init__(self):
        self.common_calls = {}

    def make(self, alpha, beta, length, keep=True):
        key = (alpha, beta, length)
        out_series = self.common_calls.get(key, np.array([]))
        if not out_series.any():
            out_series = self.impulse_response_function(alpha, beta, length)
            if keep:
                self.common_calls[key] = out_series
        return out_series

    @staticmethod
    def impulse_response_function(alpha, beta, length):
        def gamma_distribution(t, a, b):
            a, b = map(float, (a, b))
            tau = a * b
            return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)

        return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


class Waterbody:
    def __init__(self, region, comid, name):
        self.region = region
        self.id = comid
        self.name = name


class Reach(Waterbody):
    def __init__(self, region, reach_id=None, path_address=None, actual_id=None):

        super(Reach, self).__init__(region, reach_id, actual_id)
        self.mode = region.mode
        self.path_address = path_address
        self.n_dates = region.n_dates

        self.n_paths = path_address[1] - path_address[0]

        self.wb = 0
        self.outlet = None
        self.run = False

        self._type = None

        self.daysheds = 0
        self.times_convolved = 0

    def aliased(self, alias):
        new_copy = copy.copy(self)
        new_copy.name = alias
        return new_copy

    @staticmethod
    def compute_concentration(totals, baseflow):
        mass_time_series, runoff_time_series = totals
        total_flow = runoff_time_series + baseflow
        concentration = np.nan_to_num(mass_time_series / total_flow)
        runoff_concentration = np.nan_to_num(mass_time_series / runoff_time_series)

        return total_flow, concentration, runoff_concentration

    def local_sam(self, upstream_reaches):
        mass_and_runoff = self.region.sam_output[:2, upstream_reaches]
        baseflow = self.region.sam_output[2][self.id]
        return mass_and_runoff, baseflow

    def local_upstream(self):

        start_row, end_row, column = self.path_address

        all_reaches = self.region.paths[start_row:end_row, column:]

        all_times = self.region.times[start_row:end_row, column:]

        start_time = all_times[0, 0]

        addresses = (all_reaches != 0)
        reaches = all_reaches[addresses]
        # self.region.count += len(reaches)  # Commented out by Jon F. (Trip said it was old diagnostic code)
        reach_times = all_times[addresses].copy()
        reach_times -= float(start_time)
        reach_times = self.region.round_func(reach_times / tot.interval) * tot.interval

        return reaches, reach_times

    def process(self):

        upstream_reaches, upstream_times = self.local_upstream()

        mass_and_runoff, baseflow = self.local_sam(upstream_reaches)

        if self.mode in ("convolved", "unconvolved"):
            totals = self.upstream_flowing(self.region.irf, upstream_times, upstream_reaches, mass_and_runoff)
        elif self.mode == "aggregated":
            totals = np.sum(mass_and_runoff, axis=1)
        else:
            sys.exit("Invalid processing mode {}".format(self.mode))

        total_flow, concentration, runoff_conc = self.compute_concentration(totals, baseflow)

        if self.region.write_to_file:
            write.daily_tot(self.region.output_path, self.region.id, self.name, self.mode, self.region.dates,
                            total_flow, totals, baseflow, concentration, runoff_conc)

        self.run = True

    def upstream_flowing(self, impulse_response, upstream_times, upstream_reaches, mass_and_runoff):
        """ Get all unique upstream 'daysheds' and match with reaches """

        output_time_series = np.zeros((2, self.n_dates))  # (mass/runoff, dates)

        self.daysheds = np.max(upstream_times) + 1  # Log the number of daysheds for this reaach
        # Get all tanks upstream
        # for tank in range(np.max(upstream_times) + 1):
        for tank in range(self.daysheds):

            reaches_in_tank = upstream_reaches[upstream_times == tank]

            tank_time_series = self.region.sam_output[:2, reaches_in_tank].sum(axis=1)

            if tank > 0:
                self.times_convolved += 1  # Number of times convolution is run counter
                if self.mode == "convolved":  # Only perform convolution if timestep is not 0
                    irf = impulse_response.make(tank + 1, 1, self.n_dates)  # Get the convolution function
                    tank_time_series[0] = np.convolve(tank_time_series[0], irf)[:self.n_dates]  # Convolve mass
                    tank_time_series[1] = np.convolve(tank_time_series[1], irf)[:self.n_dates]  # Convolve runoff

                elif self.mode == "unconvolved":  # If convolution is off, just offset the time series for the tank.
                    offset = np.zeros((2., tank))
                    tank_time_series = np.hstack((offset, tank_time_series))[:, :self.n_dates]
                else:
                    sys.exit("Invalid convolution mode: {}".format(self.mode))

            output_time_series += tank_time_series  # Add the convolved tank time series to the total for the reach

        return output_time_series


class Reservoir(Waterbody):
    def __init__(self, region, comid, name, max_flow, outlet, volume, residence_time):
        super(Reservoir, self).__init__(region, comid, name)
        self.max_flow = max_flow
        self.outlet = outlet
        self.volume = volume
        self.residence_time = residence_time

        self.upstream_reaches = np.array([], dtype=np.int32)

    def process(self, convolve_runoff=False):
        """Convolve each reach upstream of a reservoir"""

        # Impulse response function for reservoir
        irf = self.region.irf.make(1, self.residence_time, self.region.n_dates)
        for reach_id in self.upstream_reaches:

            self.region.sam_output[0, reach_id] = \
                np.convolve(self.region.sam_output[0, reach_id], irf)[:self.region.n_dates]

            if convolve_runoff:
                self.region.sam_output[1, reach_id] = \
                    np.convolve(self.region.sam_output[1, reach_id], irf)[:self.region.n_dates]
            else:
                mean_runoff = np.mean(self.region.sam_output[0, reach_id])
                self.region.sam_output[1, reach_id] = np.repeat(mean_runoff, self.region.n_dates)


class Region:
    def __init__(self, inputs, paths, irf, write_to_file=False):

        self.id = inputs.region

        self.irf = irf

        self.mode = inputs.mode

        self.write_to_file = write_to_file

        self.output_path = paths.tot_output_path

        # Read in upstream paths
        self.paths, self.times, self.path_map, self.conversion_dict = read.upstream(paths.upstream_path.format(self.id))

        # Backward conversion dict
        self.convert_back = dict(zip(self.conversion_dict.values(), self.conversion_dict.keys()))

        # Read in compressed SAM output, lake file, and upstream path data
        self.sam_output, self.dates = self.read_sam_output(paths.sam_output_path.format(inputs.sam_output_id))

        # Populate reach index
        self.reaches = self.map_reaches()

        # Read in data from lake file
        self.waterbodies, self.wb_lookup, self.wb_conversion = self.read_lake_file(paths.lakefile_path)

        # Upstream lentics
        self.read_upstream_lentics(paths.lentics_path)

        # Set whether travel times will be rounded up or down to the nearest day
        # self.round_func = np.int32 if tot.round_down else np.vectorize(lambda x: np.int32(np.round(x)))

    @staticmethod
    def round_func(x):
        if tot.round_down:
            return np.int32(x)
        else:
            return np.int32(np.round(x))

    def cascade(self):

        # Sort lakes into bins based on how many lakes lie upstream
        lake_bins = defaultdict(set)
        for lake in self.waterbodies.values():
            upstream_lakes = np.unique(self.wb_lookup[lake.upstream_reaches])
            upstream_lakes = upstream_lakes[(upstream_lakes > 0) & (upstream_lakes != lake.id)]
            lake_bins[upstream_lakes.size].add(lake)

        # TODO: Commented out by Jon F. and replaced with code below next TODO statement
        # for lake_bin, lakes in sorted(lake_bins.items()):
        #     for lake in lakes:
        #         reaches = filter(None, map(self.reaches.get, lake.upstream_reaches))
        #         yield reaches, lake

        # TODO: Altered by Jon F. to split up processing over each Lake Bin
        for lake_bin, lakes in sorted(lake_bins.items()):
            yield lake_bin, lakes, len(lake_bins[lake_bin])

        yield None, None, 0  # No more lake_bins or lakes, return None for both to process the remaining reaches

    def map_reaches(self):

        reaches = {}

        # Populate reach index
        for reach_id, path_address in enumerate(self.path_map):
            if reach_id and path_address.any():
                actual_id = self.convert_back.get(reach_id)
                reaches[reach_id] = Reach(self, reach_id, path_address, actual_id)
        reaches.pop(0, None)

        return reaches

    def read_lake_file(self, lakefile_path):

        lake_file = lakefile_path.format(self.id)

        # Read table from file
        waterbody_table = pd.read_csv(lake_file, header=0).as_matrix()

        # Trim table to lakes that exceed the minimum residence time
        residence_time_col = 5
        active_lakes = (waterbody_table[:, residence_time_col] >= tot.minimum_residence_time)
        waterbody_table = waterbody_table[active_lakes, :]

        # 1-d array where waterbody_lookup[reach_id] = waterbody_id
        waterbody_lookup = np.zeros(max(self.reaches.keys()) + 1, dtype=np.int32)

        # Dictionary where waterbody_conversion[comid] = row num in waterbody_table
        # JCH - maybe get rid of this in future versions
        waterbody_conversion = dict(zip(waterbody_table[:, 0], range(waterbody_table.shape[0])))

        # Dictionary where waterbody_dict[row_num] = Reservoir object
        waterbody_dict = {}
        for i, waterbody in enumerate(waterbody_table):
            comid, reaches, maxflow, outlet_id, volume, residence_time = waterbody
            outlet_id = self.conversion_dict.get(outlet_id)
            new_waterbody = Reservoir(self, i, comid, maxflow, outlet_id, volume, residence_time)
            waterbody_dict[i] = new_waterbody
            for reach_id in map(int, reaches.split(",")):
                reach_alias = self.conversion_dict.get(reach_id)
                reach = self.reaches.get(reach_alias)
                if reach:
                    waterbody_lookup[reach_alias] = i
                    self.reaches[reach_alias].wb = new_waterbody

        return waterbody_dict, waterbody_lookup, waterbody_conversion

    def read_sam_output(self, sam_output_file):

        output_matrix, lookup_dict, dates = read.unpickle(sam_output_file)

        # Rearrange output matrix such that the row number corresponds to the reach alias
        # JCH - can we do this in preprocessing?
        max_reach_id = max(self.conversion_dict.values())
        new_matrix = np.zeros((3, max_reach_id + 1, output_matrix.shape[2]))  # n_dates
        for comid, row_num in lookup_dict.items():
            alias = self.conversion_dict.get(comid)
            if alias:
                new_matrix[:, alias, :] = output_matrix[:, row_num, :]

        return new_matrix, dates

    def read_upstream_lentics(self, lentics_path):
        lentics_file = lentics_path.format(self.id)
        data = read.unpickle(lentics_file)
        for waterbody, upstream_reaches in data.items():
            wb_index = self.wb_conversion.get(int(waterbody))
            if wb_index:
                self.waterbodies[wb_index].upstream_reaches = np.array(list(upstream_reaches), dtype=np.int32)

    @property
    def n_dates(self):
        return self.sam_output.shape[2]


def gis_it(comids, field="COMID"):
    print("\"{}\" = ".format(field) + " OR \"{}\" = ".format(field).join(map(str, comids)))
