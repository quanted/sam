import os
import json
import math

import numpy as np
import pandas as pd

from collections import Iterable, OrderedDict

from tempfile import mkstemp
from json import encoder
from numba import guvectorize, njit
import logging


class MemoryMatrix(object):
    """ A wrapper for NumPy 'memmap' functionality which allows the storage and recall of arrays from disk """

    def __init__(self, dimensions, dtype=np.float32, path=None, existing=False):
        self.dtype = dtype
        self.path = path
        self.existing = existing

        # Initialize dimensions of array
        self.dimensions = tuple(np.array(d) if isinstance(d, Iterable) else d for d in dimensions)
        self.labels = tuple(d if isinstance(d, Iterable) else None for d in self.dimensions)
        self.shape = tuple(int(d.size) if isinstance(d, Iterable) else d for d in self.dimensions)
        self.lookup = None if self.labels[0] is None else dict(zip(self.labels[0], np.arange(self.shape[0])))

        self.initialize_array()

    def fetch(self, index, aliased=False, copy=False, verbose=True):
        array = self.copy if copy else self.reader
        try:
            output = array[index if aliased else self.lookup.get(index)]
        except IndexError:
            if verbose:
                print("{} not found".format(index))
            output = None
        del array
        return output

    def initialize_array(self):

        # Raise an error if there should be an existing matrix but isn't
        if self.existing:  #
            if not os.path.exists(self.path):
                raise Exception('Specified MemoryMatrix {} not found'.format(self.path))

        # If a path is given but not expected to exist, create a matrix at that location.
        # Otherwise, use a temporary file
        else:
            if self.path is None:
                self.path = mkstemp(suffix=".dat", dir=os.path.join("..", "bin", "temp"))[1]
            if not os.path.exists(self.path):
                np.memmap(self.path, dtype=self.dtype, mode='w+', shape=self.shape)  # Allocate memory

    def fetch_multiple(self, indices, copy=False, verbose=False, aliased=True, return_index=False, columns=None):

        # If selecting by aliases, get indices for aliases
        index = None
        if aliased:
            addresses = np.int32([self.lookup.get(x, -1) for x in indices])
            found = np.where(addresses >= 0)[0]
            not_found = indices.size - found.size
            if not_found:
                index = found
                if verbose:
                    print("Missing {} of {} indices in {} matrix".format(not_found, len(addresses), self.name))
            indices = addresses[found]

        # Fetch data from memory map
        array = np.memmap(self.path, dtype='float32', mode='c' if copy else 'r+', shape=self.shape)
        out_array = array[indices] if columns is None else array[np.ix_(indices, columns)]
        del array

        if return_index:
            return out_array, index
        else:
            return out_array

    def update(self, key, value, aliased=True):
        array = self.writer
        array_index = self.lookup.get(key) if aliased else key
        if array_index is not None:
            array[array_index] = value
        else:
            print("Index {} not found in {} array".format(key, self.name))
        del array

    @property
    def reader(self):
        return np.memmap(self.path, dtype=self.dtype, mode='r+', shape=self.shape)

    @property
    def copy(self):
        return np.memmap(self.path, dtype=self.dtype, mode='c', shape=self.shape)

    @property
    def writer(self):
        mode = 'r+' if os.path.isfile(self.path) else 'w+'
        return np.memmap(self.path, dtype=self.dtype, mode=mode, shape=self.shape)


class Geometry(object):
    def __init__(self, region, geometry_dir):
        self.region = region
        self.flowline_path = os.path.join(geometry_dir, "region_{}_flowlines".format(self.region))
        self.intake_path = os.path.join(geometry_dir, "region_{}_intakes".format(self.region))
        self.key_path = self.flowline_path + "_key.npz"

        self._intakes = None
        self._map = None
        self._shape = None

    def load_key(self):
        data = np.load(self.key_path)
        data_map = {comid: (start_row, end_row) for comid, start_row, end_row in data['map']}
        shape = tuple(data['shape'])
        return data_map, shape

    def fetch(self, comid, feature_type, verbose=False):
        if feature_type == "flowline":
            address = self.map.get(comid)
            if address is None and verbose:
                print("Reach {} not found in geometry file for Region {}".format(comid, self.region))
            elif address is not None:
                array = np.memmap(self.flowline_path + ".dat", dtype=np.float32, shape=self.shape, mode='r+')
                start_row, end_row = address
                coordinates = array[start_row:end_row].tolist()
                del array
                return coordinates
        elif feature_type == "intake":
            return self.intakes.get(comid)
        else:
            raise ValueError("Invalid feature_type {} provided".format(feature_type))

    def read_intakes(self):
        intake_table = pd.read_csv(self.intake_path + ".csv").set_index("COMID")
        self._intakes = {comid: row.to_dict() for comid, row in intake_table.iterrows()}

    def index(self, feature_type):
        if feature_type == "flowline":
            return set(self.map.keys())
        elif feature_type == "intake":
            return set(self.intakes.keys())
        else:
            raise ValueError("Invalid feature_type {} provided".format(feature_type))


    @property
    def intakes(self):
        if self._intakes is None:
            self.read_intakes()
        return self._intakes

    @property
    def map(self):
        if self._map is None:
            self._map, self._shape = self.load_key()
        return self._map

    @property
    def shape(self):
        if self._shape is None:
            self._map, self._shape = self.load_key()
        return self._shape


class HydroTable(pd.DataFrame):
    def __init__(self, region, path):
        super().__init__()
        self.region = region
        self.path = path
        self.table_path = os.path.join(path, "region_{}.npz".format(region))

        data, header = self.read_table()
        super(HydroTable, self).__init__(data=data, columns=header)

        index_col = 'wb_comid' if 'wb_comid' in self.columns else 'comid'
        self.set_index(index_col, inplace=True)
        self.index = np.int32(self.index)

    def read_table(self):
        assert os.path.isfile(self.table_path), "Table not found for region {} in {}".format(self.region, self.path)
        data = np.load(self.table_path)
        return data['table'], data['key']

    def flows(self, reach, month_index):
        return self.loc[reach].as_matrix()[month_index]


class Hydroregion(object):
    """
    Contains all datasets and functions related to the NHD Plus region, including all hydrological features and links
    """

    def __init__(self, region, sim_type, map_dir, flowfile_dir, upstream_dir, lakefile_dir, geometry_dir):
        self.id = region
        self.sim_type = sim_type
        self.feature_type = "flowline" if self.sim_type == 'eco' else "intake"

        # Read hydrological input files
        self.geometry = Geometry(self.id, geometry_dir)
        if self.id == 'mtb':
            self.id = '07'
        self.flow_file = HydroTable(self.id, flowfile_dir)
        self.lake_table = HydroTable(self.id, lakefile_dir)
        self.nav = Navigator(self.id, upstream_dir)
        self.recipe_map = RecipeMap(self.id, map_dir)

        # Confine to available reaches and assess what's missing
        self.geometry_index = self.geometry.index(self.feature_type)
        self.active_reaches = self.confine()
        self.run_reaches = set()

    def cascade(self):

        from .parameters import time_of_travel as tot

        # Identify all outlets (and therefore, reservoirs) that exist in the current run scope
        confined_outlets = self.lake_table.outlet_comid.isin(self.active_reaches)

        # Confine lake table to these outlets and sort by number of upstream reservoirs
        confined_lake_table = self.lake_table[confined_outlets]
        confined_lake_table = confined_lake_table[confined_lake_table.residence_time >= tot.minimum_residence_time]
        confined_lake_table.sort_values('n_upstream', inplace=True)

        # Loop through reservoirs in a downstream direction and process all upstream reaches
        for i, lake in confined_lake_table.iterrows():
            upstream_reaches, warning = self.nav.upstream_watershed(lake.outlet_comid, mode='reach', return_times=False)
            upstream_reaches = set(upstream_reaches)
            reaches = upstream_reaches & self.active_reaches - {lake.outlet_comid}
            yield reaches - self.run_reaches, lake
            self.run_reaches |= reaches
        remaining_reaches = self.active_reaches - self.run_reaches
        yield remaining_reaches, None

    def confine(self):

        # Recipe/reaches that are (1) in the upstream file and (2) have a recipe file in at least 1 yr
        map_reaches = {reach_id for reach_id, _ in self.recipe_map.map.keys()}
        region_reaches = self.nav.reach_ids
        active_reaches = map_reaches & region_reaches

        # Get reaches from geometry file and all associated upstream
        upstream_reaches = \
            {us for r in self.geometry.index(self.feature_type)
             for us in self.nav.upstream_watershed(r, return_times=False)[0]}
        active_reaches &= upstream_reaches
        active_reaches.discard(0)

        return active_reaches


class ImpulseResponseMatrix(MemoryMatrix):
    """ A matrix designed to hold the results of an impulse response function for 50 day offsets """

    def __init__(self, n_dates, size=50):
        self.n_dates = n_dates
        self.size = size
        super(ImpulseResponseMatrix, self).__init__([size, n_dates])
        for i in range(size):
            irf = impulse_response_function(i, 1, self.n_dates)
            self.update(i, irf)

    def fetch(self, index):
        if index <= self.size:
            irf = super(ImpulseResponseMatrix, self).fetch(index, verbose=False)
        else:
            irf = impulse_response_function(index, 1, self.n_dates)
        return irf


class InputParams(object):
    """
    User-specified parameters and parameters derived from them.
    This class is used to hold parameters and small datasets that are global in nature and apply to the entire model
    run including Endpoints, Crops, Dates, Intake reaches, and Impulse Response Functions
    """

    def __init__(self, input_dict):
        from .parameters import time_of_travel

        # Read input dictionary
        self.__dict__.update(input_dict)

        # JCH - enables running just the Mark Twain
        if self.region == "Mark Twain Demo":
            self.region = 'mtb'
        elif self.region in ('7.0', 7.0):
            self.region = '07'

        # JCH - this will be used to run multiple regions eventually (change frontend)
        self.active_regions = [self.region]

        # Endpoints
        self.endpoints = pd.DataFrame(self.endpoints)

        # Crops
        self.crops = set(self.applications[:, 0])

        # Dates
        self.dates = pd.date_range(self.sim_date_start, self.sim_date_end)
        self.n_dates = len(self.dates)
        self.year_index = self.dates.year - self.dates.year[0]
        self.month_index = self.dates.month - 1
        self.unique_years, self.year_length = np.unique(self.dates.year, return_counts=True)
        self.new_year = np.int32([(np.datetime64("{}-01-01".format(year)) - self.sim_date_start).astype(int)
                                  for year in self.unique_years])

        # Initialize an impulse response matrix if convolving timesheds
        self.irf = None if not time_of_travel.gamma_convolve else ImpulseResponseMatrix(self.dates.size)

        # Read token
        try:
            self.token = self.csrfmiddlewaretoken
        except AttributeError:
            print("No token provided: output will be written to directory \'dummy\'")
            self.token = 'dummy'

        # To be added to front end
        self.write_contributions = True
        self.write_exceedances = True
        self.write_time_series = False
        self.read_overlay = False


class Navigator(object):
    def __init__(self, region_id, upstream_path):
        self.file = os.path.join(upstream_path, "region_{}.npz".format(region_id))
        self.paths, self.times, self.map, self.alias_to_reach, self.reach_to_alias = self.load()
        self.reach_ids = set(self.reach_to_alias.keys())

    def load(self):
        assert os.path.isfile(self.file), "Upstream file {} not found".format(self.file)
        data = np.load(self.file, mmap_mode='r')
        conversion_array = data['alias_index']
        reverse_conversion = dict(zip(conversion_array, np.arange(conversion_array.size)))
        return data['paths'], data['time'], data['path_map'], conversion_array, reverse_conversion

    def upstream_watershed(self, reach_id, mode='reach', return_times=True):

        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        # Look up reach ID and fetch address from upstream object
        reach = reach_id if mode == 'alias' else self.reach_to_alias.get(reach_id)
        reaches, adjusted_times, warning = np.array([]), np.array([]), None
        try:
            start_row, end_row, col = map(int, self.map[reach])
            start_col = list(self.paths[start_row]).index(reach)
        except TypeError:
            warning = "Reach {} not found in region".format(reach)
        except ValueError:
            warning = "{} not in upstream lookup".format(reach)
        else:
            # Fetch upstream reaches and times
            aliases = unpack(self.paths)
            reaches = aliases if mode == 'alias' else np.int32(self.alias_to_reach[aliases])
        if not return_times:
            return reaches, warning
        else:
            if warning is None:
                times = unpack(self.times)
                adjusted_times = np.int32(times - self.times[start_row][start_col])
            return reaches, adjusted_times, warning


class RecipeMap(MemoryMatrix):
    def __init__(self, region_id, recipe_path):
        self.region = region_id
        self.path = os.path.join(recipe_path, "region_{}".format(self.region))
        self.key_path = self.path + "_key.npz"

        self.map, self.shape, self.scenarios = self.load_key()

    def load_key(self):
        data = np.load(self.key_path)
        scenario_index, comid_table, shape = data['scenarios'], data['map'], data['shape']
        comid_map = {(comid, year): (start_row, end_row) for year, comid, start_row, end_row in comid_table}

        return comid_map, tuple(shape), scenario_index

    def fetch(self, comid, year, aliased=False, verbose=False):
        address = self.map.get((comid, year))
        if address is None and verbose:
            print("Reach {} not found in recipe map for Region {}".format(comid, self.region))
        elif address is not None:
            array = np.memmap(self.path + ".dat", np.int32, shape=self.shape)
            start_row, end_row = address
            areas, aliases = array[start_row:end_row].T
            scenario_ids = self.scenarios[aliases]
            del array
            return scenario_ids, areas


class Recipes(object):
    def __init__(self, i, o, year, region, scenarios, output_path, active_reaches):
        self.i = i
        self.o = o
        self.year = year
        self.region = region
        self.output_dir = os.path.join(output_path, i.token)
        self.scenario_matrix = scenarios.processed_matrix
        self.recipe_ids = sorted(region.active_reaches)
        self.outlets = set(self.recipe_ids) & set(self.region.lake_table.outlet_comid)
        self.active_reaches = active_reaches
        self.processed = set()

        # Initialize local matrix: matrix of local runoff and mass, for rapid internal recall
        self.local = MemoryMatrix([self.recipe_ids, 2, self.i.n_dates])

    def burn_reservoir(self, lake, upstream_reaches):

        from .parameters import time_of_travel as time_of_travel_params

        if lake is not None and upstream_reaches:

            # Get the convolution function
            irf = impulse_response_function(1, lake.residence_time, self.i.n_dates)

            # Pull mass and runoff time series for all upstream reaches and add together
            if upstream_reaches:

                upstream_reaches = np.array(list(upstream_reaches))

                old_mass, old_runoff = self.local.fetch_multiple(upstream_reaches).sum(axis=0)

                # Modify combined time series to reflect reservoir
                new_mass = np.convolve(old_mass, irf)[:self.i.n_dates]
                if time_of_travel_params.convolve_runoff:  # Convolve runoff
                    new_runoff = np.convolve(old_runoff, irf)[:self.i.n_dates]
                else:  # Flatten runoff
                    new_runoff = np.repeat(np.mean(old_runoff), self.i.n_dates)

                # Add all lake mass and runoff to outlet
                self.local.update(lake.outlet_comid, np.array([new_mass, new_runoff]))

    def fetch_scenarios(self, recipe_id):
        """  Fetch all scenarios and multiply by area.  For erosion, area is adjusted. """
        try:
            scenarios, areas = self.region.recipe_map.fetch(recipe_id, self.year, aliased=True).T
        except AttributeError:
            return None, None
        else:
            # Pull data from memmap. Config is (scenarios, vars, dates)
            data = self.scenario_matrix.fetch_multiple(scenarios, copy=True, aliased=False)

            # Adjust axes so that configuration is (var, date, scenario)
            data = np.rollaxis(data, 0, 3)

            # Modify values based on area of scenario
            data[:2] *= areas  # runoff, runoff_mass.
            data[2:] *= np.power(areas / 10000., .12)  # erosion, erosion_mass

            return scenarios, data

    def local_loading(self, recipe_id, cumulative, process_benthic=False):
        """ Pull all scenarios in recipe from scenario matrix and adjust for area """

        # Sum time series from all scenarios
        runoff, erosion, runoff_mass, erosion_mass = cumulative

        # Run benthic/water column partitioning
        benthic_conc = self.partition_benthic(recipe_id, erosion, erosion_mass) if process_benthic else None

        # If it's a lake outlet, do this
        if recipe_id in self.outlets:
            runoff_mass, runoff = np.array([runoff_mass, runoff]) + self.local.fetch(recipe_id)

        return runoff_mass, runoff, benthic_conc

    def partition_benthic(self, recipe_id, erosion, erosion_mass):
        """ Compute concentration in the benthic layer based on mass of eroded sediment """

        from .parameters import benthic

        surface_area = self.region.flow_file.loc[recipe_id]["surface_area"]
        soil_volume = benthic.depth * surface_area
        pore_water_volume = soil_volume * benthic.porosity
        benthic_mass = benthic_loop(erosion, erosion_mass, soil_volume)
        return benthic_mass / pore_water_volume

    def process_recipes(self, recipe_ids, progress_interval=1000):

        for recipe_id in recipe_ids:

            # Determine whether to do additional analysis on recipe
            active_recipe = recipe_id in self.active_reaches

            self.processed.add(recipe_id)
            if not len(self.processed) % progress_interval:
                print("Processed {} of {} recipes".format(len(self.processed), len(self.recipe_ids)))

            # Fetch time series data from all scenarios in recipe
            scenarios, time_series = self.fetch_scenarios(recipe_id)  # (var, date, scenario)

            if scenarios is not None:

                # Assess the contributions to the recipe from each source (runoff/erosion) and crop
                self.o.update_contributions(recipe_id, scenarios, time_series[[1, 3]].sum(axis=1))

                # Process local contributions
                local_mass, local_runoff, benthic_conc = \
                    self.local_loading(recipe_id, time_series.sum(axis=2), active_recipe)

                # Update local array with mass and runoff
                self.local.update(recipe_id, np.array([local_mass, local_runoff]))

                # Upstream processing and output generation only done if the recipe is in the write list
                if active_recipe:

                    # Process upstream contributions
                    total_flow, total_runoff, total_mass, total_conc = self.upstream_loading(recipe_id)

                    if total_conc is not None:
                        # Calculate exceedances
                        self.o.update_exceedances(recipe_id, total_conc)

                        # Store results in output array
                        self.o.update_time_series(recipe_id, total_flow, total_runoff, total_mass, total_conc,
                                                  benthic_conc)

    def upstream_loading(self, reach):
        """ Identify all upstream reaches, pull data and offset in time """
        from .parameters import time_of_travel

        # Fetch all upstream reaches and corresponding travel times
        reaches, times, warning = self.region.nav.upstream_watershed(reach)
        indices = np.int16([i for i, r in enumerate(reaches) if r not in self.region.run_reaches])
        reaches, times = reaches[indices], times[indices]

        if len(reaches) > 1:  # Don't need to do this if it's a headwater
            # JCH - Check here to confirm that all upstream reaches have been processed?
            mass_and_runoff, index = self.local.fetch_multiple(reaches, return_index=True)  # (reaches, vars, dates)
            if index is not None:
                reaches, times = reaches[index], times[index]
            totals = np.zeros((2, self.i.n_dates))  # (mass/runoff, dates)
            for tank in range(np.max(times) + 1):
                in_tank = mass_and_runoff[times == tank].sum(axis=0)
                if tank > 0:
                    if time_of_travel.gamma_convolve:
                        irf = self.region.irf.fetch(tank)  # Get the convolution function
                        in_tank[0] = np.convolve(in_tank[0], irf)[:self.i.n_dates]  # mass
                        in_tank[1] = np.convolve(in_tank[1], irf)[:self.i.n_dates]  # runoff
                    else:
                        in_tank = np.pad(in_tank[:, :-tank], ((0, 0), (tank, 0)), mode='constant')
                totals += in_tank  # Add the convolved tank time series to the total for the reach
            mass, runoff = totals
        else:
            mass, runoff = self.local.fetch(reach)

        flow = self.region.flow_file.flows(reach, self.i.month_index)
        if flow is not None:
            total_flow, (concentration, runoff_conc) = \
                compute_concentration(mass, runoff, self.i.n_dates, flow)
            return total_flow, runoff, mass, concentration
        else:
            return None, None, None, None


class Scenarios(object):
    def __init__(self, i, region, input_memmap_path, active_reaches='all', recipe_map=None, retain=None):
        self.i = i

        # JCH - temporary, for demo
        if region == 'mtb':
            region = '07'

        self.path = os.path.join(input_memmap_path, "region_" + region)
        self.keyfile_path = self.path + "_key.txt"
        self.active_reaches = active_reaches
        self.recipe_map = recipe_map

        # Load key file for interpreting input matrices
        self.arrays, self.variables, self.names, self.array_shape, self.variable_shape, self.start_date, \
        self.n_dates = self.load_key()

        # Calculate date offsets
        self.end_date = self.start_date + np.timedelta64(self.n_dates, 'D')
        self.start_offset, self.end_offset = self.date_offsets()

        # Confine to active reaches
        if active_reaches != 'all' and self.recipe_map is not None:
            self.names = self.confine()

        # Initialize input matrices
        self.array_matrix = \
            MemoryMatrix([self.names, self.arrays, self.n_dates], path=self.path + "_arrays.dat", existing=True)
        self.variable_matrix = \
            MemoryMatrix([self.names, self.variables], path=self.path + "_vars.dat", existing=True)

        # Initialize empty matrix for procesed scenarios
        process = retain is None or not os.path.exists(retain)
        self.processed_matrix = MemoryMatrix([self.names, 4, i.n_dates], path=retain)
        if process:
            self.process_scenarios()

    def confine(self):
        years = sorted((k[1] for k in self.recipe_map.map.keys()))
        names = {s for y in years for r in self.active_reaches for sc, _ in self.recipe_map.fetch(r, y) for s in sc}
        return sorted(names)

    def date_offsets(self):
        if self.start_date > self.i.sim_date_start:
            self.i.sim_date_start = self.start_date
            print("Simulation dates start earlier than the dates in the scenario matrix. Adjusting run dates")
        if self.end_date < self.i.sim_date_end:
            self.i.sim_date_end = self.end_date
            print("Simulation dates end after the dates in scenario matrix. Adjusting run dates")
        start_offset = (self.i.sim_date_start - self.start_date).astype(int)
        end_offset = (self.i.sim_date_end - self.end_date).astype(int) + 1
        if not end_offset:
            end_offset = self.n_dates + 1
        return start_offset, end_offset

    def load_key(self):
        with open(self.keyfile_path) as f:
            arrays, variables, scenarios = (next(f).strip().split(",") for _ in range(3))
            start_date = np.datetime64(next(f).strip())
            shape = np.array([int(val) for val in next(f).strip().split(",")])
        return arrays, variables, np.array(scenarios), shape[:3], shape[3:], start_date, int(shape[2])

    def process_scenarios(self, chunk=5000, progress_interval=2500):

        from .parameters import soil, plant

        # Initialize readers and writers
        # Reminder: array_matrix.shape, processed_matrix.shape = (scenario, variable, date)
        array_reader, variable_reader = self.array_matrix.reader, self.variable_matrix.reader
        processed_writer = self.processed_matrix.writer

        # Iterate scenarios
        for n, scenario_id in enumerate(self.names):

            # Report progress and reset readers/writers at intervals
            if not (n + 1) % progress_interval:
                print("{}/{}".format(n + 1, len(self.names)))

            # Open and close read/write cursors at intervals. This seems to help
            if not n % chunk:
                if n != 0:
                    del array_reader, variable_reader, processed_writer
                    array_reader, variable_reader = self.array_matrix.reader, self.variable_matrix.reader
                    processed_writer = self.processed_matrix.writer

            # Get crop ID of scenario and find all associated crops in group
            crop = int(scenario_id.split("cdl")[1])

            if crop in self.i.crops:

                # Read non-sequential variables from scenario and process
                if self.i.read_overlay:
                    covmax, org_carbon, bulk_density, overlay, *plant_dates = variable_reader[n]
                else:
                    overlay = None

                    covmax, org_carbon, bulk_density, *plant_dates = variable_reader[n]

                # Extract arrays
                leaching, runoff, erosion, soil_water, plant_factor, rain = \
                    array_reader[n][:, self.start_offset:self.end_offset]

                # Set kd depending on settings
                kd = self.i.koc * org_carbon if self.i.kd_flag else self.koc

                # Assert that all data is the proper shape for use in the functions
                assert len(plant_dates) == 5, "Looking for 5 planting dates, found {}".format(len(plant_dates))
                assert self.i.applications.shape[1] == 11, "Invalid application matrix, should have 11 columns"

                # Calculate the daily input of pesticide to the soil and plant canopy
                application_mass = \
                    pesticide_to_field(self.i.applications, self.i.new_year, crop, plant_dates, rain, n == 47646)

                # Calculate the daily mass of pesticide in the soil
                pesticide_mass_soil = pesticide_to_soil(application_mass, rain, plant_factor, soil.cm_2,
                                                        plant.deg_foliar, plant.washoff_coeff, covmax)

                # Determine the loading of pesticide into runoff and eroded sediment
                runoff_mass, erosion_mass = \
                    pesticide_to_water(pesticide_mass_soil, runoff, erosion, leaching, bulk_density, soil_water, kd,
                                       self.i.deg_aqueous, soil.runoff_effic, soil.delta_x, soil.erosion_effic,
                                       soil.soil_depth)

                # If the scenario is an overlay, runoff and erosion are not to be added to totals
                if overlay == 1:
                    runoff[:] = 0.
                    erosion[:] = 0.

                # Write runoff, erosion, runoff mass and erosion mass to processed matrix
                processed_writer[n] = runoff, runoff_mass, erosion, erosion_mass

            else:
                # Write runoff and erosion
                processed_writer[n, [0, 2]] = array_reader[n][1:3, self.start_offset:self.end_offset]

        del array_reader, variable_reader


class Outputs(object):
    def __init__(self, i, scenario_ids, output_path, geometry, feature_type, demo_mode=False):

        logging.info("SAM TASK Generating Outputs... ")
        logging.info("SAM Outputs part 1")
        self.i = i
        self.geometry = geometry
        self.scenario_ids = scenario_ids
        self.output_dir = os.path.join(output_path, i.token)
        self.feature_type = feature_type
        self.recipe_ids = sorted(self.geometry.index(self.feature_type))
        self.demo_mode = demo_mode
        logging.info("SAM Outputs part 2")
        # Initialize output matrices
        self.output_fields = ['total_flow', 'total_runoff', 'total_mass', 'total_conc', 'benthic_conc']
        self.time_series = MemoryMatrix([self.recipe_ids, self.output_fields, self.i.n_dates])
        logging.info("SAM Outputs part 3")

        # Initialize contributions matrix: loading data broken down by crop and runoff v. erosion source
        self.exceedances = MemoryMatrix([self.recipe_ids, self.i.endpoints.shape[0]])
        logging.info("SAM Outputs part 4")

        # Initialize contributions matrix: loading data broken down by crop and runoff v. erosion source
        self.contributions = MemoryMatrix([self.recipe_ids, 2, self.i.crops])
        self.contributions.columns = np.int32(sorted(self.i.crops))
        self.contributions.header = ["cls" + str(c) for c in self.contributions.columns]
        logging.info("SAM Outputs Completed")


    def update_contributions(self, recipe_id, scenario_index, loads):
        """ Sum the total contribution by land cover class and add to running total """

        scenario_names = self.scenario_ids[scenario_index]
        classes = [int(name.split("cdl")[1]) for name in scenario_names]
        contributions = np.zeros((2, 255))
        for i in range(2):  # Runoff Mass, Erosion Mass
            contributions[i] += np.bincount(classes, weights=loads[i], minlength=255)
        self.contributions.update(recipe_id, contributions[:, self.contributions.columns])

    def update_exceedances(self, recipe_id, concentration):
        durations, endpoints = self.i.endpoints[["duration", "endpoint"]].as_matrix().T
        exceed = exceedance_probability(concentration, *map(np.int16, (durations, endpoints, self.i.year_index)))
        self.exceedances.update(recipe_id, exceed)

    def update_time_series(self, recipe_id, total_flow=None, total_runoff=None, total_mass=None, total_conc=None,
                           benthic_conc=None):

        writer = self.time_series.writer
        index = self.time_series.index.index(recipe_id)

        # This must match self.fields as designated in __init__
        rows = [total_flow, total_runoff, total_mass, total_conc, benthic_conc]
        for i, row in enumerate(rows):
            if row is not None:
                writer[index, i] = row
        del writer

    def write_json(self, write_exceedances=False, write_contributions=False):

        # Initialize JSON output
        encoder.FLOAT_REPR = lambda o: format(o, '.4f')
        out_file = os.path.join(self.output_dir, "out_json.csv")
        out_json = OrderedDict((("type", "FeatureCollection"), ("features", [])))

        # Iterate through recipes
        for recipe_id in self.recipe_ids:

            # Initialize new feature
            new_feature = self.geometry.fetch(recipe_id, self.feature_type)
            if self.feature_type == "flowline":
                new_feature = {"type": "Feature",
                               "geometry": {"type": "LineString", "coordinates": new_feature},
                               "properties": {"COMID": int(recipe_id)}}
            else:
                coordinates = [float(new_feature.pop(f)) for f in ("y_coord", "x_coord")]
                new_feature = {"type": "Feature",
                               "geometry": {"type": "Point", "coordinates": coordinates},
                               "properties": new_feature}

            # Add exceedance probabilities
            if write_exceedances:
                exceedance_dict = dict(zip(self.i.endpoints.short_name, map(float, self.exceedances.fetch(recipe_id))))
                new_feature['properties'].update(exceedance_dict)

            # Add percent contributions
            if write_contributions:
                contributions = self.contributions.fetch(recipe_id)
                for i, category in enumerate(("runoff", "erosion")):
                    labels = ["{}_load_{}".format(category, label) for label in self.contributions.header]
                    contribution_dict = dict(zip(labels, map(float, contributions[i])))
                    new_feature['properties'].update(contribution_dict)

            # Append new feature to output
            out_json["features"].append(new_feature)

        # Convert output dict to JSON object
        out_json = json.dumps(out_json, sort_keys=False, separators=(',', ':'))

        # Write to file
        with open(out_file, 'w') as f:
            f.write(out_json)

    def write_demo(self, fields):
        # Initialize JSON output
        encoder.FLOAT_REPR = lambda o: format(o, '.4f')
        out_file = os.path.join(self.output_dir, "out_json.csv")
        out_json = OrderedDict((("type", "FeatureCollection"), ("features", [])))

        # Iterate through recipes
        for recipe_id in self.recipe_ids:

            # Initialize new feature
            new_feature = self.geometry.fetch(recipe_id, self.feature_type)
            if self.feature_type == "flowline":
                new_feature = {"type": "Feature",
                               "geometry": {"type": "LineString", "coordinates": new_feature},
                               "properties": {"COMID": int(recipe_id)}}
            else:
                coordinates = [float(new_feature.pop(f)) for f in ("y_coord", "x_coord")]
                new_feature = {"type": "Feature",
                               "geometry": {"type": "Point", "coordinates": coordinates},
                               "properties": new_feature}

            # Add exceedance probabilities
            for field_name, sigma in fields:
                new_feature["properties"][field_name] = min((abs(np.random.normal(0, sigma)), 1))

            # Append new feature to output
            out_json["features"].append(new_feature)

        # Convert output dict to JSON object
        out_json = json.dumps(out_json, sort_keys=False, separators=(',', ':'))

        # Write to file
        with open(out_file, 'w') as f:
            f.write(out_json)

    def write_output(self):

        # Create output directory
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)

        # Write JSON output
        if self.demo_mode:

            demo_fields = [['acute_human', 0.1],
                           ['acute_fw_fish', 0.1],
                           ['chronic_fw_fish', 0.2],
                           ['acute_fw_inv', 0.1],
                           ['chronic_fw_inv', 0.3],
                           ['acute_em_fish', 0.05],
                           ['chronic_em_fish', 0.2],
                           ['acute_em_inv', 0.2],
                           ['chronic_em_inv', 0.3],
                           ['acute_nonvasc_plant', 0.01],
                           ['acute_vasc_plant', 0.1],
                           ['pct_corn', 0.4],
                           ['pct_soy', 0.2],
                           ['pct_row', 0.1]]

            self.write_demo(demo_fields)
        else:
            self.write_json(self.i.write_exceedances, self.i.write_contributions)

        # Write time series
        if self.i.write_time_series:
            self.write_time_series()

    def write_time_series(self, fields='all'):
        if fields == 'all':
            fields = self.output_fields
            field_indices = np.arange(len(self.output_fields))
        else:
            field_indices = [self.output_fields.index(field) for field in fields]

        heading_lookup = {'runoff_conc': 'RunoffConc(ug/L)',
                          'local_runoff': 'LocalRunoff(m3)',
                          'total_runoff': 'TotalRunoff(m3)',
                          'local_mass': 'LocalMass(m3)',
                          'total_flow': 'TotalFlow(m3)',
                          'baseflow': 'Baseflow(m3)',
                          'total_conc': 'TotalConc(ug/L)',
                          'total_mass': 'TotalMass(kg)',
                          'wc_conc': 'WC_Conc(ug/L)',
                          'erosion': 'Erosion(kg)',
                          'erosion_mass': 'ErodedMass(kg)',
                          'runoff_mass': 'RunoffMass(kg)',
                          'benthic_conc': 'BenthicConc(ug/L)'}

        headings = [heading_lookup.get(field, "N/A") for field in fields]

        for recipe_id in self.recipe_ids:
            out_file = os.path.join(self.output_dir, "time_series_{}.csv".format(recipe_id))
            out_data = self.time_series.fetch(recipe_id)[field_indices].T
            df = pd.DataFrame(data=out_data, index=self.i.dates, columns=headings)
            df.to_csv(out_file)


def initialize():
    # Make sure needed subdirectories exist
    preexisting_subdirs = ["Results", "temp"]
    for subdir in preexisting_subdirs:
        d = os.path.join("..", "bin", subdir)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            if subdir == 'temp':
                for f in os.listdir(d):
                    try:
                        os.remove(os.path.join(d, f))
                    except PermissionError:
                        print("Unable to get permission to delete temp file")


@njit
def benthic_loop(eroded_soil, erosion_mass, soil_volume):
    benthic_mass = np.zeros(erosion_mass.size, dtype=np.float32)
    benthic_mass[0] = erosion_mass[0]
    for i in range(1, erosion_mass.size):
        influx_ratio = eroded_soil[i] / (eroded_soil[i] + soil_volume)
        benthic_mass[i] = (benthic_mass[i - 1] * (1. - influx_ratio)) + (erosion_mass[i] * (1. - influx_ratio))
    return benthic_mass


def compute_concentration(transported_mass, runoff, n_dates, q):
    """ Concentration function for time of travel """
    mean_runoff = runoff.mean()  # m3/d
    baseflow = np.subtract(q, mean_runoff, out=np.zeros(n_dates), where=(q > mean_runoff))
    total_flow = runoff + baseflow
    concentration = np.divide(transported_mass, total_flow, out=np.zeros(n_dates),
                              where=(total_flow != 0))
    runoff_concentration = np.divide(transported_mass, runoff, out=np.zeros(n_dates),
                                     where=(runoff != 0))

    return total_flow, map(lambda x: x * 1000000., (concentration, runoff_concentration))  # kg/m3 -> ug/L


@guvectorize(['void(float64[:], int16[:], int16[:], int16[:], float64[:])'], '(p),(o),(o),(p)->(o)')
def exceedance_probability(time_series, window_sizes, endpoints, years_since_start, res):
    # Count the number of times the concentration exceeds the test threshold in each year
    n_years = years_since_start.max()
    for test_number in range(window_sizes.size):
        window_size = window_sizes[test_number]
        threshold = endpoints[test_number]
        window_sum = np.sum(time_series[:window_size])
        exceedances = np.zeros(n_years)
        for day in range(window_size, len(time_series)):
            year = years_since_start[day]
            window_sum += time_series[day] - time_series[day - window_size]
            avg = window_sum / window_size
            if avg > threshold:
                exceedances[year] = 1
        res[test_number] = exceedances.sum() / n_years


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        try:
            return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)
        except Exception as e:
            print(a, b, tau, t, length)
            print(e)
            raise

    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


@njit
def pesticide_to_field(applications, new_years, active_crop, event_dates, rain, diagnostic=False):
    """ Simulate timing of pesticide appplication to field """

    application_mass = np.zeros((2, rain.size))
    for i in range(applications.shape[0]):

        crop, event, offset, canopy, step, window1, pct1, window2, pct2, effic, rate = applications[i]
        diagnostic = False
        if diagnostic:
            print("sam.functions.pesticide to field, diagnostic is set to true, crop = active crop")
            print("print applications.shape")
            print(i, applications.shape[0])
            print("print active crop")
            print(crop, crop == active_crop)
        if crop == active_crop:
            event_date = int(event_dates[int(event)])
            daily_mass_1 = rate * (pct1 / 100.) / window1
            if diagnostic:
                print("sam.functions.pesticide to field, diagnostic is set to true, calculating mass")
            if step:
                daily_mass_2 = rate * (pct2 / 100.) / window2
            # if diagnostic:
            #    print("b")
            for year in range(new_years.size):
                new_year = new_years[year]
                # print(year, new_year, int(window1), int(window2))
                for k in range(int(window1)):
                    date = int(new_year + event_date + offset + k)
                    # print(k, int(canopy), new_year, event_date, offset, date, daily_mass_1)
                    # print(application_mass[int(canopy), date])
                    application_mass[int(canopy), date] = daily_mass_1
                if step:
                    for l in range(int(window2)):
                        date = int(new_year + event_date + window1 + offset + l)
                        application_mass[int(canopy), date] = daily_mass_2
            if diagnostic:
                print("diagnostic = true, finished with this application")
    return application_mass


@njit
def pesticide_to_soil(application_mass, rain, plant_factor, soil_2cm, foliar_degradation, washoff_coeff, covmax):
    """ Calcluate pesticide in soil and simulate movement of pesticide from canopy to soil """

    # Initialize output
    pesticide_mass_soil = np.zeros(rain.size)
    canopy_mass, last_application = 0, 0  # Running variables

    # Determine if any pesticide has been applied to canopy
    canopy_applications = application_mass[1].sum() > 0

    # Loop through each day
    for day in range(plant_factor.size):
        # Start with pesticide applied directly to soil
        pesticide_mass_soil[day] = application_mass[0, day] * soil_2cm
        # If pesticide has been applied to the canopy, simulate movement from canopy to soil
        if canopy_applications:
            if application_mass[1, day] > 0:  # Pesticide applied to canopy on this day
                canopy_pesticide_additions = application_mass[1, day] * plant_factor[day] * covmax
                pesticide_mass_soil[day] += (application_mass[1, day] - canopy_pesticide_additions) * soil_2cm
                canopy_mass = canopy_pesticide_additions + \
                              canopy_mass * np.exp((day - last_application) * foliar_degradation)
                last_application = day
            if rain[day] > 0:  # Simulate washoff
                canopy_mass *= np.exp((day - last_application) * foliar_degradation)
                pesticide_remaining = canopy_mass * np.exp(-rain[day] * washoff_coeff)
                pesticide_mass_soil[day] += canopy_mass - pesticide_remaining
                last_application = day  # JCH - sure?
    return pesticide_mass_soil


@njit
def pesticide_to_water(pesticide_mass_soil, runoff, erosion, leaching, bulk_density, soil_water, kd, deg_aqueous,
                       runoff_effic, delta_x, erosion_effic, soil_depth):
    # Initialize output arrays
    runoff_mass, erosion_mass = np.zeros(runoff.size), np.zeros(runoff.size)

    # Initialize running variables
    total_mass, degradation_rate = 0, 0

    # Initialize erosion intensity
    erosion_intensity = erosion_effic / soil_depth

    # Loop through days
    for day in range(runoff.size):
        daily_runoff = runoff[day] * runoff_effic
        total_mass = total_mass * degradation_rate + pesticide_mass_soil[day]
        retardation = (soil_water[day] / delta_x) + (bulk_density * kd)
        deg_total = deg_aqueous + ((daily_runoff + leaching[day]) / (delta_x * retardation))
        if leaching[day] > 0:
            degradation_rate = np.exp(-deg_total)
        else:
            degradation_rate = np.exp(-deg_aqueous)

        average_conc = ((total_mass / retardation / delta_x) / deg_total) * (1 - degradation_rate)

        if leaching[day] > 0:
            runoff_mass[day] = average_conc * daily_runoff
        if erosion[day] > 0:
            enrich = np.exp(2.0 - (0.2 * np.log10(erosion[day])))
            enriched_eroded_mass = erosion[day] * enrich * kd * erosion_intensity * 0.1
            erosion_mass[day] = average_conc * enriched_eroded_mass

    return runoff_mass, erosion_mass


if __name__ == "__main__":
    import cProfile
    from .pesticide_calculator import main

    if False:
        cProfile.run('main()')
    else:
        main()
