import os
import sys
import math

import numpy as np
import pandas as pd

from numba import njit, guvectorize
from tempfile import mkdtemp
from collections import namedtuple


class MemoryMatrix(object):
    def __init__(self, index, y_size, z_size=None, path=None, base=None, name=None, dtype=np.float32, overwrite=True,
                 input_only=False):
        self.name = name
        self.dtype = dtype
        self.index = index
        self.count = len(self.index)
        self.lookup = dict(zip(self.index, np.arange(self.count)))
        self.shape = (self.count, y_size) if not z_size else (self.count, y_size, z_size)
        overwrite = False if input_only else overwrite

        # Load from saved file if one is specified, else generate
        path = mkdtemp() if not path else path
        if not os.path.exists(path):
            os.makedirs(path)
        base = name if not base else base
        self.path = os.path.join(path, base + ".dat")
        assert not input_only or os.path.exists(self.path), "Matrix {} not found".format(self.path)
        if overwrite or not os.path.exists(self.path):
            np.memmap(self.path, dtype=dtype, mode='w+', shape=self.shape)
            self.new = True
        else:
            self.new = False

    def fetch(self, get_index, verbose=True, raw_index=False, dtype=None, copy=False):
        array = self.copy if copy else self.reader
        location = self.lookup.get(get_index) if not raw_index else get_index
        if location is not None:
            output = array[location]
        else:
            if verbose:
                print("{} not found".format(get_index))
            output = None
        del array
        return output

    def fetch_multiple(self, indices, copy=False, verbose=False, aliased=True, return_index=False, columns=None):

        mode = 'c' if copy else 'r'
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

        array = np.memmap(self.path, dtype='float32', mode=mode, shape=self.shape)
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
        return np.memmap(self.path, dtype=self.dtype, mode='r', shape=self.shape)

    @property
    def copy(self):
        return np.memmap(self.path, dtype=self.dtype, mode='c', shape=self.shape)

    @property
    def writer(self):
        mode = 'r+' if os.path.isfile(self.path) else 'w+'
        return np.memmap(self.path, dtype=self.dtype, mode=mode, shape=self.shape)


class FlowMatrix(MemoryMatrix):
    def __init__(self, file_path, region, dates):
        self.region = region
        self.dates = dates
        self.months = np.array(list(map(lambda x: x.month, self.dates))) - 1

        # Load key
        self.key_file = os.path.join(file_path, "region_{}_key.npy".format(self.region))
        key = np.load(self.key_file, mmap_mode='r')
        self.shape, self.reaches = tuple(key[:2]), key[2:]

        # Initialize matrix
        super(FlowMatrix, self).__init__(self.reaches, self.shape[1], name="region_{}".format(self.region),
                                         path=file_path, input_only=True)

    def fetch(self, reach_id, verbose=False):
        row = super(FlowMatrix, self).fetch(reach_id, verbose)
        if row is not None:
            length = row[1] * 1000.  # km -> m
            q, v = np.split(row[2:], [12])
            return namedtuple('Flow', ['q', 'v', 'l'])(q[self.months], v[self.months], length)
        else:
            pass


class Hydroregion(object):
    """
    Contains all datasets and functions related to the NHD Plus region, including all hydrological features and links
    """

    def __init__(self, i, map_path, flowfile_dir, upstream_dir, lakefile_dir):
        from .parameters import time_of_travel

        self.i = i
        self.id = i.region
        self.irf = None if not time_of_travel.gamma_convolve else ImpulseResponseMatrix(i.n_dates)

        # Read hydrological input files
        self.flow_file = FlowMatrix(flowfile_dir, self.id, i.dates)
        self.nav = Navigator(self.id, upstream_dir)
        self.lake_table = pd.read_csv(lakefile_dir.format(self.id), index_col="LakeID")
        self.recipe_map = RecipeMap(map_path)
        self.years = sorted({year for _, year in self.recipe_map.lookup.keys() if year})

        # Confine to available reaches and assess what's missing
        self.active_reaches = self.confine()
        self.run_reaches = set()

    def cascade(self):
        confined_lake_table = self.lake_table[self.lake_table.OutletID.isin(self.active_reaches)]
        for i, lake in confined_lake_table.iterrows():
            upstream_reaches = set(self.nav.upstream_watershed(lake.OutletID, mode='reach', return_times=False))
            reaches = upstream_reaches & self.active_reaches - {lake.OutletID}
            yield reaches - self.run_reaches, lake
            self.run_reaches |= reaches
        remaining_reaches = self.active_reaches - self.run_reaches
        yield remaining_reaches, None

    def confine(self):
        # Recipe/reaches that are (1) in the upstream file and (2) have a recipe file in at least 1 yr
        map_reaches = {reach_id for reach_id, _ in self.recipe_map.lookup.keys()}
        region_reaches = set(self.nav.reach_to_alias.keys())
        active_reaches = map_reaches & region_reaches
        active_reaches.discard(0)

        # Identify full watershed extent of reaches and get matching lakes
        # full_watershed = {us for r in active_reaches for us in self.upstream_watershed(r, return_times=False)}

        return active_reaches


class ImpulseResponseMatrix(MemoryMatrix):
    def __init__(self, n_dates):
        self.index = np.arange(1000)
        self.n_dates = n_dates
        super(ImpulseResponseMatrix, self).__init__(np.arange(50), n_dates, name="impulse response")

    def fetch(self, index):
        irf = super(ImpulseResponseMatrix, self).fetch(index, verbose=False)
        if irf is None:
            irf = impulse_response_function(index, 1, self.n_dates)
            self.update(index, irf)
        return irf


class InputFile(object):
    """ User-specified parameters and parameters derived from them """

    def __init__(self, input_data):

        # Adjust variables
        from .parameters import to_be_added_params

        self.format_inputs(input_data)

        # Crops
        self.crops = {application.crop for application in self.applications}

        # Convert half-lives to degradation rates
        self.deg_aqueous, self.deg_photolysis, self.deg_hydrolysis, self.deg_wc, self.deg_benthic = \
            map(lambda x: 0.693 / x if x else 0.,
                (
                    self.soil_hl, self.aq_photolysis_hl, self.hydrolysis_hl, self.wc_metabolism_hl,
                    self.ben_metabolism_hl))
        self.koc /= 1000.0  # now in m3/kg

        self.region = '07' if self.region == 'mtb' else self.region  # temporary, for utool pilot

        # Add in hardwired stuff that will eventually go in front end
        self.__dict__.update(to_be_added_params)

    def format_inputs(self, data):

        def application_matrix(applications):
            new_applications = []
            Application = namedtuple("Application", sorted({key for app in applications for key in app.keys()}))
            for application_dict in applications:
                application_dict['rate'] /= 10000.  # kg/ha -> kg/m2
                new_applications.append(Application(**application_dict))
            return new_applications

        def date(datestring):
            m, d, y = datestring.split("/")
            return np.datetime64("{}-{}-{}".format(y, m, d))

        input_format = \
            {"chemical_name": str,  # Atrazine
             "region": str,  # Ohio Valley
             "applications": application_matrix,
             "endpoints": dict,
             "soil_hl": float,  # Soil half life
             "wc_metabolism_hl": float,  # Water column metabolism half life
             "ben_metabolism_hl": float,  # Benthic metabolism half life
             "aq_photolysis_hl": float,  # Aqueous photolysis half life
             "hydrolysis_hl": float,  # Hydrolysis half life
             "kd_flag": int,  # 1
             "koc": float,  # 100
             "sim_date_start": date,  # 01/01/1984
             "sim_date_end": date,  # 12/31/2013
             }

        # Check if any required input data are missing or extraneous data are provided
        provided_fields = set(data.keys())
        required_fields = set(input_format.keys())
        unknown_fields = provided_fields - required_fields
        missing_fields = required_fields - provided_fields
        if unknown_fields:
            print("Input field(s) \"{}\" not understood".format(", ".join(unknown_fields)))

        assert not missing_fields, "Required input field(s) \"{}\" not provided".format(", ".join(missing_fields))

        input_data = {field: input_format[field](data[field]) for field in input_format.keys()}
        self.__dict__.update(input_data)


    @property
    def dates(self):
        return pd.date_range(self.sim_date_start, self.sim_date_end)

    @property
    def year(self):
        return self.dates.year

    @property
    def year_length(self):
        unique_years, days = np.unique(self.year, return_counts=True)
        return days

    @property
    def new_year(self):
        return [(np.datetime64("{}-01-01".format(year)) - self.sim_date_start).astype(int)
                for year in self.unique_years]

    @property
    def n_dates(self):
        return len(self.dates)


class Navigator(object):
    def __init__(self, region_id, upstream_path):
        self.file = upstream_path.format(region_id)
        self.paths, self.times, self.map, self.alias_to_reach, self.reach_to_alias = self.load()

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

        try:
            start_row, end_row, col = map(int, self.map[reach])
            start_col = list(self.paths[start_row]).index(reach)
        except TypeError:
            print("Reach {} not found in region".format(reach))
            return [] if not return_times else [], []
        except ValueError:
            print("{} not in upstream lookup".format(reach))
            return [] if not return_times else [], []
        else:
            # Fetch upstream reaches and times
            aliases = unpack(self.paths)
            reaches = aliases if mode == 'alias' else np.int32(self.alias_to_reach[aliases])
        if not return_times:
            return reaches
        else:
            times = unpack(self.times)
            adjusted_times = np.int32(times - self.times[start_row][start_col])
            return reaches, adjusted_times


class RecipeMap(MemoryMatrix):
    def __init__(self, in_path):
        self.dir, self.base = os.path.split(in_path)
        self.memmap_file = in_path + ".dat"
        self.key_file = in_path + "_key.npy"
        self.name = "recipe map"

        self.key = list(map(tuple, np.load(self.key_file)))
        self.n_cols = self.key.pop(-1)[0]

        super(RecipeMap, self).__init__(self.key, self.n_cols, 2, name=self.name, path=self.dir, base=self.base,
                                        dtype=np.int32, input_only=True)

    def fetch(self, get_index, verbose=False):
        results = super(RecipeMap, self).fetch(get_index, verbose)
        if results is not None:
            results = results[results[:, 0] > [0]]
            if not results.any():
                results = None
        return results


class RecipeMatrices(object):
    def __init__(self, i, year, region, scenario_matrix, output_path, write_list=None, save_file=None):
        self.i = i
        self.year = year
        self.region = region
        self.output_dir = os.path.join(output_path, i.chemical_name)
        self.scenario_matrix = scenario_matrix
        self.filter = write_list
        self.recipe_ids = sorted(region.active_reaches)
        self.outlets = set(self.recipe_ids) & set(self.region.lake_table.OutletID)

        # Initialize output matrix: matrix containing SAM results that don't require recall in model
        self.output_fields = ['benthic_conc', 'total_flow', 'total_runoff', 'total_mass', 'total_conc']
        self.output = MemoryMatrix(self.recipe_ids, len(self.output_fields), self.i.n_dates, path=self.output_dir,
                                   base=save_file, name="output")

        # Initialize contributions matrix: loading data broken down by crop and runoff v. erosion source
        self.contributions = np.zeros((2, 255))

        # Initialize local matrix: matrix of local runoff and mass, for rapid internal recall
        self.local_fields = ['local_mass', 'local_runoff']
        self.local = MemoryMatrix(self.recipe_ids, len(self.local_fields), self.i.n_dates, name="recipe")

    def burn_reservoir(self, lake, upstream_reaches):

        from .parameters import time_of_travel as time_of_travel_params

        if lake is not None and upstream_reaches:

            # Get the convolution function
            irf = impulse_response_function(1, lake.ResidenceTime, self.i.n_dates)

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
                self.local.update(lake.OutletID, np.array([new_mass, new_runoff]))

    def calculate_exceedances(self, concentration):
        if concentration is not None:
            categories = [("Human health DWLOC (ug/L)", 'human', (4, 21, 60)),
                          ("Freshwater Fish (Tox x LOC)", 'fw_fish', (4, 60,)),
                          ("Freshwater Invertebrate (Tox x LOC)", 'fw_inv', (4, 21,)),
                          ("Estuarine/Marine Fish (Tox x LOC)", 'em_fish', (4, 21,)),
                          ("Estuarine/Marine Invertebrate (Tox x LOC)", 'em_inv', (4, 21,)),
                          ("Aquatic nonvascular plant (Tox x LOC)", 'nonvasc_plant', (4, 21,)),
                          ("Aquatic vascular plant (Tox x LOC)", 'vasc_plant'), (4, 21,)]

            durations, endpoints, names = [], [], []
            for label, field, durations in categories:
                for level, duration in zip(('acute', 'chronic', 'human'), durations):
                    endpoint = self.i.endpoints(["{}_{}".format(level, field)])
                    if endpoint:
                        durations.append(duration)
                        endpoints.append(endpoint)
                        names.append(label)

            years = self.i.dates.year - self.i.dates.year.min()
            exceed = moving_window(concentration, *map(np.int16, (durations, endpoints, years, self.i.year_length)))

        return exceed

    def calculate_contributions(self, scenarios, masses):
        # Sum the total contribution by land cover class and add to running total
        scenarios = np.array(scenarios)
        scenario_names = self.scenario_matrix.scenarios[scenarios]
        classes = [int(name.split("cdl")[1]) for name in scenario_names]
        self.contributions[0] += np.bincount(classes, weights=masses[0], minlength=255)
        self.contributions[1] += np.bincount(classes, weights=masses[1], minlength=255)

    def compute_concentration(self, reach, transported_mass, runoff):
        """ Concentration function for time of travel """
        try:
            q = self.region.flow_file.fetch(reach).q
        except AttributeError:
            return None, (None, None)
        else:
            """
            JCH - the idea here is that baseflow is estimated as the difference between total predicted q (erom) and
            mean modeled runoff. Predicted baseflow isn't event specific and is not sensitive to runoff itself apart
            from the mean.  This balances the sum of the modeled Q with the sum of the predicted Q
            """
            mean_runoff = runoff.mean()  # m3/d
            baseflow = np.subtract(q, mean_runoff, out=np.zeros(self.i.n_dates), where=(q > mean_runoff))
            total_flow = runoff + baseflow
            concentration = np.divide(transported_mass, total_flow, out=np.zeros(self.i.n_dates),
                                      where=(total_flow != 0))
            runoff_concentration = np.divide(transported_mass, runoff, out=np.zeros(self.i.n_dates),
                                             where=(runoff != 0))

            return total_flow, map(lambda x: x * 1000000., (concentration, runoff_concentration))  # kg/m3 -> ug/L

    def fetch_scenarios(self, recipe_id):
        try:
            data = self.region.recipe_map.fetch((recipe_id, self.year)).T
            scenarios, areas = data
        except AttributeError:
            return None, None
        else:
            # Fetch all scenarios and multiply by area.  For erosion, area is adjusted.
            array = self.scenario_matrix.processed_matrix.fetch_multiple(scenarios, copy=True, aliased=False)
            array = np.swapaxes(array, 0, 2)
            array[:, :2] *= areas  # runoff, runoff_mass
            array[:, 2:] *= (areas / 10000.) ** .12  # erosion, erosion_mass
            return scenarios, array

    def local_loading(self, recipe_id, cumulative):
        """Pull all scenarios in recipe from scenario matrix and adjust for area"""

        # Sum time series from all scenarios
        runoff, runoff_mass, erosion, erosion_mass = cumulative.sum(axis=2).T

        # Run benthic/water column partitioning
        benthic_conc = self.partition_benthic(recipe_id, runoff, runoff_mass,
                                              erosion_mass) if self.i.process_benthic else None

        # Add together all masses and update the array
        transported_mass = runoff_mass + erosion_mass

        if recipe_id in self.outlets:
            transported_mass, runoff = np.array([transported_mass, runoff]) + self.local.fetch(recipe_id)

        return transported_mass, runoff, benthic_conc

    def partition_benthic(self, reach, runoff, runoff_mass, erosion_mass):
        from .parameters import soil, stream_channel, benthic

        try:
            reach = self.region.flow_file.fetch(reach)
            q, v, l = reach.q, reach.v, reach.l
        except AttributeError:
            return None, None, (None, None)

        mean_runoff = runoff.mean()  # m3/d
        baseflow = np.subtract(q, mean_runoff, out=np.zeros(self.i.n_dates), where=(q > mean_runoff))
        total_flow = runoff + baseflow
        mixing_cell = 40.  # meters
        cross_section = total_flow / v
        width = stream_channel.a * np.power(cross_section, stream_channel.b)
        depth = cross_section / width
        surface_area = width * l
        volume = np.array([(depth * surface_area),  # Water column
                           (benthic.depth * surface_area * benthic.porosity)])  # Benthic zone

        # Compute concentration in runoff of runoff mass and erosion mass
        runoff_conc = np.divide(runoff_mass, runoff, out=np.zeros(self.i.n_dates), where=(runoff != 0))
        daily_conc = np.divide(runoff_mass + erosion_mass, mixing_cell, out=np.zeros(self.i.n_dates),
                               where=(runoff_mass + erosion_mass > 0.0) & (mixing_cell > 0.0))

        # Divide mass loading between water column and benthic zones
        mass_input = np.vstack([runoff_mass + ((1. - soil.prben) * erosion_mass),  # Water Column
                                soil.prben * erosion_mass]).T  # Benthic
        # Partition concentration into benthic and water column concentrations
        # This needs to be cleaned up
        # Compute benthic solute holding capacity
        fw1, fw2, theta, sed_conv_factor, omega = solute_holding_capacity(depth, surface_area, self.i.koc)

        k_adj = np.array((total_flow / mixing_cell) + (self.i.deg_photolysis + self.i.deg_hydrolysis) * fw1 + \
                         (self.i.deg_wc * fw1) + self.i.deg_benthic * (1 - fw1))

        aqconc_avg_wb, daily_avg, daily_peak = \
            concentration_loop(self.i.n_dates, daily_conc, k_adj, volume,
                               mass_input, fw1, fw2, omega, theta, self.i.deg_aqueous)

        return map(lambda x: x * 1000000., (runoff_conc, aqconc_avg_wb, daily_avg, daily_peak))

    def process_recipes(self, recipe_ids):
        for _, recipe_id in enumerate(recipe_ids):
            print(recipe_id)
            # Fetch time series data from all scenarios in recipe
            scenarios, time_series = self.fetch_scenarios(recipe_id)

            if scenarios is not None:
                # Assess the contributions to the recipe from each source (runoff/erosion) and crop
                self.calculate_contributions(scenarios, time_series[:, [1, 3]].sum(axis=0))

                # Process local contributions
                local_mass, local_runoff, benthic_conc = self.local_loading(recipe_id, time_series)

                # Update local array with mass and runoff
                self.local.update(recipe_id, np.array([local_mass, local_runoff]))

                # Process upstream contributions
                total_flow, total_runoff, total_mass, total_conc = self.upstream_loading(recipe_id)

                # Calculate exceedances
                exceedances = self.calculate_exceedances(total_conc)

                print(total_flow.sum(), total_runoff.sum(), total_mass.sum(), total_conc.sum())

                # Store results in output array
                self.update_output(recipe_id, benthic_conc, total_flow, total_runoff, total_mass, total_conc)

        return exceedances

    def update_output(self, recipe_id, benthic_conc=None, total_flow=None, total_runoff=None, total_mass=None,
                      total_conc=None):

        writer = self.output.writer
        index = self.recipe_ids.index(recipe_id)

        # This must match self.fields as designated in __init__
        rows = [benthic_conc, total_flow, total_runoff, total_mass, total_conc]
        for i, row in enumerate(rows):
            if row is not None:
                writer[index, i] = row
        del writer

    def upstream_loading(self, reach):
        def confine_reaches():
            indices = np.int16([i for i, r in enumerate(reaches) if r not in self.region.run_reaches])
            return reaches[indices], times[indices]

        from .parameters import time_of_travel

        # Fetch all upstream reaches and corresponding travel times
        reaches, times = self.region.nav.upstream_watershed(reach)
        reaches, times = confine_reaches()
        if len(reaches) > 1:  # Don't need to do this if it's a headwater
            # Check here to confirm that all upstream_reaches have been processed?
            mnr, index = self.local.fetch_multiple(reaches, return_index=True)  # (reaches, mass/runoff, dates)
            if index is not None:
                reaches, times = reaches[index], times[index]
            totals = np.zeros((2, self.i.n_dates))  # (mass/runoff, dates)
            for tank in range(np.max(times) + 1):
                in_tank = mnr[times == tank].sum(axis=0)
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

        total_flow, (concentration, runoff_conc) = self.compute_concentration(reach, mass, runoff)

        return total_flow, runoff, mass, concentration

    def write_contributions(self):
        out_file = os.path.join(self.output_dir, "{}_contributions.csv".format(self.i.chemical_name))
        active_crops = np.where(self.contributions.sum(axis=0) > 0)[0]
        df = pd.DataFrame(data=self.contributions[:, active_crops].T, index=active_crops, columns=["Runoff", "Erosion"])
        df.to_csv(out_file)

    def write_time_series(self, recipe_id, fields='all'):
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
        out_data = self.output.fetch(recipe_id)[field_indices].T
        df = pd.DataFrame(data=out_data, index=self.i.dates, columns=headings)
        df.to_csv(os.path.join(self.output_dir, "time_series_{}.csv".format(recipe_id)))


class ScenarioMatrices(object):
    def __init__(self, i, input_memmap_path, retain=None):
        self.i = i
        self.path, self.base = os.path.split(input_memmap_path)
        self.keyfile_path = input_memmap_path + "_key.txt"

        # Load key file for interpreting input matrices
        self.arrays, self.variables, self.scenarios, \
        self.array_shape, self.variable_shape, self.start_date, self.n_dates = self.load_key()

        # Calculate date offsets
        self.end_date = self.start_date + np.timedelta64(self.n_dates, 'D')
        self.start_offset, self.end_offset = self.date_offsets()

        # Initialize input matrices
        self.array_matrix = MemoryMatrix(self.scenarios, len(self.arrays), self.n_dates, self.path,
                                         self.base + "_arrays", name="arrays", input_only=True)
        self.variable_matrix = MemoryMatrix(self.scenarios, len(self.variables), path=self.path,
                                            base=self.base + "_vars", name="variables", input_only=True)
        if retain:
            self.processed_matrix = MemoryMatrix(self.scenarios, 4, i.dates.size, path=self.path, base=retain,
                                                 name="scenario", overwrite=False)
        else:
            self.processed_matrix = MemoryMatrix(self.scenarios, 4, i.dates.size, name="scenario")

        # Populate scenario matrix if it doesn't exist
        if self.processed_matrix.new:
            self.process_scenarios()

        self.inspect()

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

    def extract_scenario(self, scenario_id):

        arrays = self.array_matrix.fetch(scenario_id)[:, self.start_offset:self.end_offset]
        variables = self.variable_matrix.fetch(scenario_id)
        leaching, runoff, erosion, soil_water, plant_factor, rain = arrays
        covmax, org_carbon, bulk_density, plant_beg, harvest_beg, emerg_beg, bloom_beg, mat_beg = variables
        events = {'plant': plant_beg, 'emergence': emerg_beg, 'maturity': mat_beg, 'harvest': harvest_beg}
        return events, runoff, erosion, leaching, org_carbon, bulk_density, soil_water, plant_factor, rain, covmax

    def inspect(self):
        """Make sure expected arrays and variables are present"""
        """Make sure array index and variable index are identical"""
        pass

    def load_key(self):
        with open(self.keyfile_path) as f:
            arrays, variables, scenarios = (next(f).strip().split(",") for _ in range(3))
            start_date = np.datetime64(next(f).strip())
            shape = np.array([int(val) for val in next(f).strip().split(",")])
        return arrays, variables, np.array(scenarios), shape[:3], shape[3:], start_date, int(shape[2])

    def pesticide_applications(self, active_crops, event_dates, plant_factor, rain, covmax):

        from .parameters import soil, plant

        scenario_applications = set(filter(lambda x: x.crop in active_crops, self.i.applications))
        application_mass = np.zeros((2, self.i.n_dates))
        canopy_applications = False

        # Determine how much pesticide is applied and when
        for app in scenario_applications:
            index = ['groud', 'foliar'].index(app.method)
            if index:
                canopy_applications = True
            start_dates = np.int16(self.i.new_year + event_dates[app.event] + app.offset)
            first_window = \
                np.repeat(start_dates, app.window1) + np.tile(np.arange(app.window1), len(start_dates))
            application_mass[index, first_window] += (app.rate * (app.pct1 / 100.)) / app.window1
            if app.refine == 'step':
                second_window = \
                    np.repeat(start_dates + app.window1, app.window2) + \
                    np.tile(np.arange(app.window2), len(start_dates))
                application_mass[index, second_window] += (app.rate * (app.pct2 / 100.)) / app.window2

        pesticide_mass_soil = application_mass[0] * soil.distrib_2cm
        if canopy_applications:
            pesticide_mass_soil += \
                canopy_loop(self.i.n_dates, application_mass[1], np.array(plant_factor), covmax,
                            soil.distrib_2cm, plant.foliar_degradation, np.array(rain), plant.washoff_coeff)

        return pesticide_mass_soil

    def pesticide_transport(self, pesticide_mass_soil, runoff, erosion, leaching, org_carbon, bulk_density, soil_water):
        """ Simulate transport of pesticide through runoff and erosion """
        from .parameters import soil

        # Initialize output
        runoff_mass = np.zeros(self.i.n_dates)
        erosion_mass = np.zeros(self.i.n_dates)
        runoff = runoff * soil.runoff_effic
        leach_dates = np.where(leaching > 0.0)[0]

        kd = self.i.koc * org_carbon if self.i.kd_flag else self.i.koc

        # Get retardation and deg_total
        retardation = (soil_water / soil.delta_x) + (bulk_density * kd)
        deg_total = self.i.deg_aqueous + ((runoff + leaching) / (soil.delta_x * retardation))

        # Get degradation rate for each day
        degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-self.i.deg_aqueous))  # Non-leaching days
        degradation_rate[leach_dates] = np.exp(-deg_total[leach_dates])  # Leaching days

        # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
        total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

        # Compute conc
        average_conc = ((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate)

        # Compute the mass of pesticide in runoff
        runoff_mass[leach_dates] = average_conc[leach_dates] * runoff[leach_dates]  # conc[kg/m3]*[m] = kg/m2

        # Compute the mass of pesticide from erosion
        erosion_dates = np.where((erosion > 0) & (leaching > 0))[0]
        erosion_intensity = soil.erosion_effic / soil.soil_depth  # Assume uniform extraction, no decline, MMF
        enrich = np.exp(2.0 - (0.2 * np.log10(erosion[erosion_dates])))
        enriched_eroded_mass = erosion[erosion_dates] * enrich * kd * erosion_intensity * 0.1
        erosion_mass[erosion_dates] = average_conc[erosion_dates] * enriched_eroded_mass

        return runoff_mass, erosion_mass

    def process_scenario(self, scenario_id, leaching, runoff, erosion, soil_water, plant_factor, rain, covmax,
                         org_carbon, bulk_density, event_dates):

        from .parameters import crop_groups

        crop = int(scenario_id.split("cdl")[1])
        all_crops = {crop} | crop_groups.get(crop, set())
        active_crops = self.i.crops & all_crops

        if active_crops:

            # Compute pesticide application that winds up in soil
            pesticide_mass_soil = \
                self.pesticide_applications(active_crops, event_dates, plant_factor, rain, covmax)

            # Determine the loading of pesticide into runoff and erosion
            runoff_mass, erosion_mass = \
                self.pesticide_transport(pesticide_mass_soil, runoff, erosion, leaching, org_carbon, bulk_density,
                                         soil_water)

        else:
            runoff_mass, erosion_mass = np.zeros((2, self.i.n_dates))

        return np.array([runoff, runoff_mass, erosion, erosion_mass])

    def process_scenarios(self, chunk=5000):
        array_reader = self.array_matrix.reader
        variable_reader = self.variable_matrix.reader
        writer = self.processed_matrix.writer
        for n, scenario_id in enumerate(self.scenarios):
            if not n % chunk:
                del array_reader, variable_reader, writer
                array_reader = self.array_matrix.reader
                variable_reader = self.variable_matrix.reader
                writer = self.processed_matrix.writer

            # Extract arrays
            leaching, runoff, erosion, soil_water, plant_factor, rain = \
                array_reader[n][:, self.start_offset:self.end_offset]

            # Extract variables
            covmax, org_carbon, bulk_density, plant_beg, harvest_beg, emerg_beg, bloom_beg, mat_beg = \
                variable_reader[n]

            events = {'plant': plant_beg, 'emergence': emerg_beg, 'maturity': mat_beg, 'harvest': harvest_beg}

            writer[n] = self.process_scenario(scenario_id, leaching, runoff, erosion, soil_water, plant_factor,
                                              rain, covmax, org_carbon, bulk_density, events)

        del array_reader, variable_reader, writer


@njit
def canopy_loop(n_dates, application_mass, plant_factor, covmax, soil_2cm, foliar_degradation, rain, washoff_coeff):
    canopy_mass = 0
    canopy_to_soil = np.zeros(n_dates)  # Cumulative
    last_application = 0
    for day in range(n_dates):
        if application_mass[day] > 0:
            canopy_pesticide_additions = application_mass[day] * plant_factor[day] * covmax
            canopy_to_soil[day] = (application_mass[day] - canopy_pesticide_additions) * soil_2cm
            canopy_mass = canopy_pesticide_additions + canopy_mass * np.exp(
                (day - last_application) * foliar_degradation)
            last_application = day
        if rain[day] > 0:
            canopy_mass *= np.exp((day - last_application) * foliar_degradation)
            pesticide_remaining = canopy_mass * np.exp(-rain[day] * washoff_coeff)
            canopy_to_soil[day] += canopy_mass - pesticide_remaining
            last_application = day
    return canopy_to_soil


@njit
def concentration_loop(n_dates, daily_concentration, k_adj, daily_volume, mass_input, fw1, fw2, omega, theta, deg_aq):
    # Beginning day aquatic concentrations, considered Peak Aqueous Daily Conc in Water Column
    daily_peak = np.zeros((2, n_dates))
    daily_avg = np.zeros((2, n_dates))
    aqconc_avg_wb = np.zeros(n_dates)

    # Reset starting values
    exp_k = np.exp(-k_adj)
    aqconc_wb = 0
    antecedent_mass = np.zeros(2)  # mn

    for day in range(daily_concentration.size):
        # Add mass input to antecedent mass
        daily_mass = antecedent_mass + mass_input[day]

        # Convert to aqueous concentrations (peak) at beginning of day
        # JCH - fw comes from solute_holding_capacity. Fraction going into each section. Should fw[0] + fw[1] = 1?
        daily_peak[0, day] = daily_mass[0] * fw1[day] / daily_volume[day, 0]
        daily_peak[1, day] = daily_mass[1] * fw2[day] / daily_volume[day, 1]

        # Compute daily average concentration in the water body - when no Benthic layer considered
        aqconc_wb += daily_concentration[day]  # initial water body concentration for current time step

        # Daily avg aq conc in water body, area under curve/t = Ci/k*(1-e^-k), NO benthic
        aqconc_avg_wb[day] = aqconc_wb / k_adj[day] * (1 - exp_k[day])

        # initial water body concentration for next time step
        aqconc_wb *= exp_k[day]

        # For simul diffeq soln: mn1,mn2,mavg1,mavg2 = new_aqconc1, new_aqconc2, aqconc_avg1[d], aqconc_avg2[d]
        # Note: aqconc_avg1 and aqconc_avg2 are outputted - Daily avg aq conc in WC and Benthic regions
        new_aqconc, wc_avg, benthic_avg = simultaneous_diffeq(k_adj[day], deg_aq, omega, theta[day], daily_peak[:, day])
        daily_avg[0, day] = wc_avg
        daily_avg[1, day] = benthic_avg

        # Masses m1 and m2 after time step, t_end
        antecedent_mass[0] = new_aqconc[0] / fw1[day] * daily_volume[day, 0]
        antecedent_mass[1] = new_aqconc[1] / fw2[day] * daily_volume[day, 1]

    return aqconc_avg_wb, daily_avg, daily_peak


@guvectorize(['void(float64[:], float64[:], float64[:])'], '(n),(o)->(n)', nopython=True)
def cumulative_multiply_and_add(a, b, res):
    res[0] = a[0]
    for i in range(1, a.size):
        res[i] = res[i - 1] * b[i - 1] + a[i]


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)

    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


@guvectorize(['void(float64[:], int16[:], int16[:], int16[:], int16[:], float64[:])'], '(p),(o),(o),(p),(n)->(o)',
             nopython=True)
def moving_window(time_series, window_sizes, endpoints, years_since_start, year_sizes, res):
    # Count the number of times the concentration exceeds the test threshold in each year
    counts = np.zeros((year_sizes.size, window_sizes.size))
    for test in range(window_sizes.size):
        window_size = window_sizes[test]
        threshold = endpoints[test]
        window_sum = np.sum(time_series[:window_size])
        for day in range(window_size, len(time_series)):
            year = years_since_start[day]
            window_sum += time_series[day] - time_series[day - window_size]
            avg = window_sum / window_size
            if avg > threshold:
                counts[year, test] += 1

    # Average the number of yearly exceedances for each test
    for test in range(window_sizes.size):
        exceedance = 0
        for year in range(year_sizes.size):
            exceedance += counts[year, test] / year_sizes[year]
        res[test] = exceedance / year_sizes.size


@njit
def simultaneous_diffeq(gamma1, gamma2, omega, theta, daily_aq_peak):
    """
    ANALYTICAL SOLUTION FOR THE TWO SIMULTANEOUS DIFFERENTIAL EQNS:
              dm1/dt = Am1 + Bm2
              dm2/dt = Em1 + Fm2
    WITH INITIAL VALUES m1 AND m2 FOR m1 AND m2
    mn1 IS OUTPUT VALUE FOR m1 AFTER TIME T
    mn2 IS OUTPUT VALUE FOR m2 AFTER TIME T
    mavg1 IS AVERAGE VALUE OF m1 OVER TIME T
    """

    t_end = 86400.  # seconds, time step of ONE DAY
    m1, m2 = daily_aq_peak

    # Calculate constants for simultaneous_diffeq: A,B,E,F
    # This reduces the model equivalent parameters to the coefficients needed for solving simultaneous_diffeq
    a = -gamma1 - omega * theta
    b = omega * theta
    e = omega
    f = -gamma2 - omega

    af = a + f
    dif = 4 * ((f * a) - (b * e))
    bbb = np.sqrt(af * af - dif)

    root1 = (af + bbb) / 2.
    root2 = (af - bbb) / 2.

    dd = (root1 - a) / b
    ee = (root2 - a) / b
    ff = ee - dd
    x1 = (ee * m1 - m2) / ff
    y1 = (m2 - dd * m1) / ff

    # Calculate new concentrations for next step
    rt1 = root1 * t_end
    rt2 = root2 * t_end
    exrt1 = np.exp(rt1)
    exrt2 = np.exp(rt2)
    ccc = x1 * exrt1
    ddd = y1 * exrt2

    # values for m1 and m2 after time step t_end
    mn = np.zeros(2)
    mn[0] = ccc + ddd  # Water column
    mn[1] = dd * ccc + ee * ddd  # Benthic

    # AVERAGE DAILY CONCENTRATION SOLUTION: set up for daily average, but can be changed by adjusting time step
    gx = x1 / root1
    hx = y1 / root2

    term1 = gx * exrt1  # term3 = -X1/root1*exp(root1*T1)
    term2 = hx * exrt2  # term4 = -Y1/root2*exp(root2*T1
    term3 = -gx
    term4 = -hx

    mavg_wc = (term1 + term2 + term3 + term4) / t_end  # Water column
    mavg_ben = (term1 * dd + term2 * ee + term3 * dd + term4 * ee) / t_end  # Benthic

    return mn, mavg_wc, mavg_ben


def solute_holding_capacity(depth, surface_area, koc):
    """Calculates Solute Holding capacities and mass transfer between water column and benthic regions"""

    from .parameters import benthic, water_column

    # Aqueous volumes in each region
    vol1 = depth * surface_area  # total volume in water column, approximately equal to water volume alone
    vol2a = benthic.depth * surface_area  # total benthic volume
    vol2 = vol2a * benthic.porosity  # total benthic pore water volume

    # Default EXAMS conditions for partitioning
    kow = koc / .35  # DEFAULT EXAMS CONDITION ON Kow  p.35
    kpdoc1 = kow * .074  # DEFAULT RELATION IN EXAMS (LITTORAL)
    kpdoc2 = koc  # DEFAULT RELATION IN EXAMS (BENTHIC) p.16 of EXAMS 2.98 (or is it Kow*.46 ?)
    xkpb = 0.436 * kow ** .907  # DEFAULT RELATION IN EXAMS

    # mass in littoral region
    vol1a = depth[0] * surface_area  # initial volume corresponding with suspended matter reference
    m_sed_1 = water_column.sused * vol1a * .001  # SEDIMENT MASS LITTORAL
    m_bio_1 = water_column.plmas * vol1a * .001  # BIOLOGICAL MASS LITTORAL
    m_doc_1 = water_column.doc * vol1a * .001  # DOC MASS LITTORAL

    # partitioning coefficients of individual media
    kd_sed_1 = koc * water_column.froc * .001  # Kd of sediment in littoral [m3/kg]
    kd_sed_2 = koc * benthic.froc * .001  # Kd of sediment in benthic
    kd_bio = xkpb / 1000.  # Kd of biological organisms
    kd_doc_1 = kpdoc1 / 1000.  # Kd of DOC in littoral region
    kd_doc_2 = kpdoc2 / 1000.  # Kd of DOC in benthic region

    # mass in benthic region
    m_sed_2 = benthic.bulk_density * vol2a * 1000.  # as defined by EXAMS parameters m_sed_2 = BULKD/PCTWA*VOL2*100000.
    m_bio_2 = benthic.bnmas * surface_area * .001
    m_doc_2 = benthic.doc * vol2 * .001

    # solute holding capacity in regions 1 and 2
    capacity_1 = kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + vol1
    capacity_2 = kd_sed_2 * m_sed_2 + kd_bio * m_bio_2 + kd_doc_2 * m_doc_2 + vol2

    # Fraction going to water column and benthic
    fw1 = vol1 / capacity_1  # fw1 is daily, vol1 is daily
    fw2 = vol2 / capacity_2

    theta = capacity_2 / capacity_1

    sed_conv_factor = vol2 / fw2 / m_sed_2  # converts pore water to [Total Conc normalized to sed mass]

    # Omega mass transfer - Calculates littoral to benthic mass transfer coefficient
    omega = benthic.d_over_dx / benthic.depth  # (m3/hr)/(3600 s/hr)

    return fw1, fw2, theta, sed_conv_factor, omega


if __name__ == "__main__":
    from pesticide_calculator import main
    main()